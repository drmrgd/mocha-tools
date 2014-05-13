#!/usr/bin/perl
# Script to collect variant call files from an Ion Torrent run into a central collectedVariants folder within the main
# results folder.  This script requires a sampleKey file indicating the barcode and the sample it represents. 
#
# 11/7/2013 - Updated for TSv4.0
#
# Created: 2/25/2013	Dave Sims
#########################################################################################################################
use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use Data::Dump;
use Cwd; 
use File::Copy;

( my $scriptname = $0 ) =~ s/^(.*\/)+//;
my $version = "v2.1.4";
my $description = <<"EOT";
Script to collect variant calls from an Ion Torrent run into a central 'collectedVariants' directory located
within the main results folder.  This script requires a sampleKey consisting of the barcode and sample it
represents, which can be generated by running 'sampleKeyGen'.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options]
    -s, --samplekey   Custom sample key file	
    -v, --version     Version Information
    -h  --help        Display the Help information
EOT

my $ver_info = 0;
my $help = 0;
my $custom_key;

GetOptions( "samplekey=s"   => \$custom_key,
            "version"       => \$ver_info,
            "help"          => \$help )
        or print $usage;

sub help {
	printf "%s - %s\n\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub version {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

help if $help;
version if $ver_info;

#########------------------------------ END ARG Parsing ---------------------------------#########

# Defaults vars
my $TVCout = "plugin_out/variantCaller_out/";

# Check the directory
my $resultsDir = getcwd;
if ( ! -d "$resultsDir/$TVCout" ) {
	die "ERROR: Either you are not running this script from a run results folder, or TVC has not yet been run on this sample\n";
}

#my ($runid) = $resultsDir =~ /([P|M]CC-\d+)/;
my ($runid) = $resultsDir =~ /((?:[P|M]CC|MC[12])-\d+)/;

# Make collectedVariants directory
my $colVarsDir = "$resultsDir/collectedVariants/";
mkdir( $colVarsDir ) unless ( -e $colVarsDir );

my $sampleKey;
# Allow for custom sample key in order to run standalone
if ( defined $custom_key ) {
	$sampleKey = $custom_key;
} else {
	$sampleKey = "$colVarsDir/sampleKey.txt";
}

# Check to make sure we have a valid sampleKey and load up a list
open( my $skey_fh, "<", "$sampleKey" ) || die "No sample key found: $!\n";
my %barcodes = map { split } <$skey_fh>;
close( $skey_fh );

# Read list of IonXpress dirs in TVC results dir
opendir( TVCRESULTS, "$TVCout" ) || die "Can't read TVC results dir: $!\n";
my @sample_dirs = grep { $barcodes{$_} } readdir( TVCRESULTS );
closedir( TVCRESULTS );

# Sanity check to make sure that the sampleKey is a good match for the data available
my $sampleNumber = 0;
$sampleNumber = scalar( @sample_dirs ) if @sample_dirs; 
my $numBC = keys %barcodes;
if ( $sampleNumber == 0  || ( $sampleNumber + 1 ) < $numBC ) { # Have to add 1 to sample number to account for LibNTC
	print "WARNING: Data not in good agreement with sample key:\n\n";
	print "\t\t'$sampleNumber' matches out of '$numBC' barcodes\n\n";
	print "Are you sure you want to continue [y/N]? ";
	chomp( my $ans = <STDIN> );
	if ( $ans =~ /^n|N$/ || $ans eq "" ) {
		die "Check for non-conforming BC index in Experiment Plan. Exiting....\n";
	} 
	elsif ($ans =~ /^y|Y$/ ) {
		print "Continuing on then. Check to be sure all samples are represented in the final output results\n";
	}
	else {
		die "Not a yes or no answer. Exiting to be safe\n";
	}
}

# Move to the variantCaller_out directory and grab copies of the variants files  Switch to the alleles.xls table as variants.xls is deprecated.
chdir( "$TVCout" );
foreach my $sample ( @sample_dirs ) {
    my $tab_file = "$sample/alleles.xls"; # replaces variants.xls
	my $vcf_file = "$sample/TSVC_variants.vcf";

	# Grab alleles.xls files
	if ( -f $tab_file ) {
		copy( $tab_file, "$colVarsDir/$sample.txt" ); #make copy of unfiltered data
		filter_vars( \$tab_file ); #filter and add to collectedVariants dir
   
	} 
	else {
		warn "The tabular variant call file: '$sample/alleles.xls' can not be found! $!\n";
	}

	# Grab vcf files
	if ( -f $vcf_file ) {
		copy( $vcf_file, "$colVarsDir/$barcodes{$sample}_${runid}.vcf" ); 
		filter_vars( \$vcf_file );
	} else {
		warn "No TSVC_variants.vcf file found.  Skipping...\n";
	}
}

# Create a summary table based on the data generated from above.
summary_table();

sub filter_vars {
	# Get rid of the 'Absent' and 'No Call' data to make more informative output.  Will keep raw data for troubleshooting later  
	
	my $varfile = shift;

	my @filtered_data;
	my ($barcode) = $$varfile =~ /^(IonXpress_\d+)/;
	my $sample_name = $barcodes{$barcode};

    # Setup the headers here for use downstream
    my $data_format = "%-7s %-12s %-8s %-16s %-6s %-12s %-12s %8.2f %12.0f %8.0f %8.0f %8.0f %5.0f     %-12s\n";
    my $theader_format = "%-7s %-12s %-8s %-16s %-6s %-12s %-12s   %-8s      %-6s    %-6s   %-7s %-8s %-5s %-12s\n";
    my @tab_fields = qw{ Chrom Position Gene AmpID Type Ref Alt Freq Qual Cov RefCov VarCov HP Hotspot };
    my $tab_header = sprintf( $theader_format, @tab_fields );

	if ( $$varfile =~ /\.vcf$/ ) {
        my $outfile = "$colVarsDir/${sample_name}_${runid}_filtered.vcf";
        return;
	}
    elsif ( $$varfile =~ /\.xls$/ ) {
		# Filter alleles.xls 
		open( my $xls_fh, "<", $$varfile ) || die "Can't open the variants.xls file for reading";
		@filtered_data =  grep { ! /(Absent|No Call)/ } <$xls_fh>; # Get only variant containing lines

        my $outfile = "$colVarsDir/${sample_name}_filtered.txt";

        open( my $out_fh, ">", $outfile ) || die "Can't create file '$outfile' for writing: $!";
        print $out_fh $tab_header;

		for ( @filtered_data ) {
			my @data = split;
            next if /Chrom/;

            my $refcov = $data[18]-$data[24];
            my @subset = ( @data[0,14,12,13,9,2,3,6,7,18], $refcov, @data[24,37,11] );
            my $formatted_data = sprintf( $data_format, @subset );

            print $out_fh $formatted_data;
		} 
    }
	return;
}

sub summary_table {
    # Create a summary table of all of the variants found in all of the samples.  Need to work with a filtered
	# dataset or this file will be a mess!
	
    my $theader_format = "%-7s %-12s %-8s %-16s %-6s %-12s %-12s   %-8s      %-6s    %-6s   %-7s %-8s %-5s %-12s\n";
    my @tab_fields = qw{ Chrom Position Gene AmpID Type Ref Alt Freq Qual Cov RefCov VarCov HP Hotspot };
    my $tab_header = sprintf( $theader_format, @tab_fields );

	my $outputFile = "$colVarsDir" .  $runid . "_allVariants.tsv";
	open ( my $summary_fh, ">", $outputFile ) || die "Output file: '$outputFile' can not be opened for writing! $!";

	opendir( VARSDIR, "$colVarsDir" ) || die "Can't read collectedVariants directory! $!";
    my @var_files = map { "$barcodes{$_}_filtered.txt" } map { $_ =~ /(IonXpress_\d+)\.txt/ } readdir( VARSDIR );

	# Read through all the variant call files, concatenate them, and add the sample name to the row.
	foreach my $varFile ( sort( @var_files ) ) {
        my ($sample_name) = $varFile =~ /^(.*)_filtered.*/;
		open ( VARFILE, "<", "$colVarsDir/$varFile" ) || die "Can't open this file for some reason: '$varFile'\n";
		print $summary_fh "::: Variant Calls for $sample_name :::\n\n";
        print $summary_fh "$tab_header";

        print $summary_fh "$sample_name\t$_" while (<VARFILE>);
        print $summary_fh "\n";
	}
    return;
}
