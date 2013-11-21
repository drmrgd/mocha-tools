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
my $version = "v0.6.1";
my $description = <<"EOT";
Script to collect variant calls from an Ion Torrent run into a central 'collectedVariants' directory located
within the main results folder.  This script requires a sampleKey consisting of the barcode and sample it
represents.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options]
	-s, --samplekey		Custom sample key file	
	-v, --version		Version Information
	-h  --help			Display the Help information
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
my ($runid) = $resultsDir =~ /([P|M]CC-\d+)/;

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
	if ( $ans =~ /^n|N/ || $ans eq "" ) {
		die "Check for non-conforming BC index in Experiment Plan. Exiting....\n";
	} 
	elsif ($ans =~ /^y|Y/ ) {
		print "Continuing on then. Check to be sure all samples are represented in the final output results\n";
	}
	else {
		die "Not a yes or no answer. Exiting to be safe\n";
	}
}

# XXX: Figure out what to keep and copy over.
# Move to the variantCaller_out directory and grab copies of the variants files
chdir( "$TVCout" );
foreach my $sample ( @sample_dirs ) {
	my $tab_file = "$sample/variants.xls";
	my $vcf_file = "$sample/TSVC_variants.vcf";
	# Grab variants.xls files
	if ( -f $tab_file ) {
		#copy( "$sample/variants.xls", "$colVarsDir/$barcodes{$sample}_variants.txt" );
		copy( $tab_file, "$colVarsDir/$sample.txt" );
		filter_vars( \$tab_file );
	} 
	else {
		warn "The variants .xls file: '$sample/variants.xls' can not be found! $!\n";
	}

	# Grab vcf files
	if ( -f $vcf_file ) {
		#copy( "$sample/TSVC_variants.vcf", "$colVarsDir/$barcodes{$sample}_variants.vcf" ); 
		copy( $vcf_file, "$colVarsDir/$sample.vcf" ); 
		filter_vars( \$vcf_file );
	} else {
		warn "No TSVC_variants.vcf file found.  Skipping...\n";
	}
}

sub filter_vars {
	# Get rid of all of the NoCall and Ref entries in the variant call files.
	
	my $varfile = shift;

	my @filtered_data;
	my ($barcode) = $$varfile =~ /^(IonXpress_\d+)/;
	my $sample_name = $barcodes{$barcode};

	if ( $$varfile =~ /\.vcf$/ ) {
		# TODO:
		# 	- Can filter VCF on 'NOCALL'
		# 	- Can also consider filtering VCF file on GT field (either ./. or 0/0 
		# 			perl -ane 'next if /^#/; print join( "\t", @F[0,1,2,3,4,5,6], ".", @F[8,9] ), "\n" if ( $F[9] !~ /^([.0]\/[.0]):/g )'

	    print "Looks like you got a VCF file on your hands buddy: $$varfile\n";
	}
    elsif ( $$varfile =~ /\.xls$/ ) {

		# Filter variants.xls file
		# TODO: 
		#     - Remove the Genotype column to shrink the table a bit.	#
		#     - Dynamically set field widths for Ref and Alt columns.
	    print "I spy with my little eye a tabular variants file: $$varfile\n";
		# TODO: Work on the output formatting.  Get dynamic field widths? Just print as is?
		#my $header = join( "\t", qw{ Chrom Position Gene AmpID Type Zygosity Genotype Ref Alt Freq Qual Coverage RefCov VarCov Hotspot } );
		my $header = qw{ Chrom Position Gene AmpID Type Zygosity Genotype Ref Alt Freq Qual Coverage RefCov VarCov Hotspot };

		open( my $xls_fh, "<", $$varfile ) || die "Can't open the variants.xls file for reading";
		@filtered_data =  grep { ! /(Ref|NC)/ } <$xls_fh>; # Get only variant containing lines
		my $format = "%-7s %-13s %-9s %-17s %-7s %-7s %-12s %-12s %-12s %-9.2d %-9.2d %-9.0d %-9.0d %-9.0d %-12s\n";

		open( my $out_fh, ">", "$colVarsDir/${sample_name}_filtered_variants.txt" ) || die "Can't create file for writing: $!";
		print $out_fh "$header\n";
		for my $line ( @filtered_data ) {
			my @data = split( /\s+/, $line );
			my $formatted_data = sprintf( "%-7s %-13s %-9s %-17s %-7s %-7s %-12s %-12s %-12s %-9.2d %-9.2d %-9.0d %-9.0d %-9.0d %-12s\n", @data );
			print $out_fh $formatted_data;
		} 
	} else {
		print "Something else: $$varfile\n";
	}
	return;
}


sub summary_table {
    # Create a summary table of all of the variants found in all of the samples.  Need to work with a filtered
	# dataset or this file will be a mess!
	
	my $outputFile = "$colVarsDir" .  $runid . "_allVariants.tsv";
	open ( OUTFILE, ">", $outputFile ) || die "Output file: '$outputFile' can not be opened for writing! $!\n";

	my @varFiles;
	opendir( VARSDIR, "$colVarsDir" ) || die "Can't read collectedVariants directory! #!\n";
	while ( my $files = readdir( VARSDIR ) ) {
		next unless ( $files =~ /IonXpress_\d+\.txt/ );
		push( @varFiles, $files );
	}

	# Read through all the variant call files, concatenate them, and add the sample name to the row.
	foreach my $varFile ( sort( @varFiles ) ) {
		( my $index = $varFile ) =~ s/\..*//;
		open ( VARFILE, "<", $varFile ) || die "Can't open this file for some reason: '$varFile'\n";
		print OUTFILE "::: Variant Calls for $barcodes{$index} :::\n\n";
		print OUTFILE "$header";

		while (<VARFILE>) {
			next if ( /Chrom/ );
			my @lines = split( /\t/, $_ );
			print OUTFILE join( "\t", $barcodes{$index}, @lines );
		}
		print OUTFILE "\n";
	}
}
