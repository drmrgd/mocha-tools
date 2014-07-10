#!/usr/bin/perl
# Using custom plasmid lookup table for version of CPSC used in experiment, query variants.xls
# file to check to see if the plasmid was seen in the sample and print out the data.  If the
# plasmid was not found, create a new list indicating which plasmids were not found. Print out
# results to a file.
#
# Created 2/24/2013 - Dave Sims
###################################################################################################		

use warnings;
use strict;
use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use Data::Dump;
use Cwd 'abs_path';
use Sort::Naturally;
use File::Basename;
use Term::ANSIColor;

my $scriptname = basename($0);
my $version = "v2.1.0_071014";
my $description = <<"EOT";
Using a plasmid lookup table for the version of the CPSC used in the experiment, query
a TVC 'alleles.xls' file to check to see if the plasmid was seen in the sample, and print
out the data.  If the plasmid was not found, create a new list indicating which plasmids
were not found.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] -l <3,4,4.1,5,mc>  <alleles.xls>

    -l, --lookup           CPSC lookup version to use (3, 4, 4.1  5, or mc)
    -V, --VCF              EXPERIMENTAL: run on vcfExtractor output instead of 'alleles.xls'
    -o, --output <file>    Output file name (default is 'cpscChecker_output.txt')
    -p, --preview          Do not write output to a file.  Print to STDOUT just for testing.
    -h, --help             Display Help information
    -v, --version          Version information
EOT

# Defaults
my $ver_info;
my $outfile = "cpscChecker_output.txt";
my $help;
my $ltable;
my $preview;
my $vcf;

GetOptions( "lookup|l=s"    => \$ltable,
            "VCF|V"         => \$vcf,
	        "preview|p"     => \$preview,
            "output|o=s"    => \$outfile,
            "version|v"     => \$ver_info,
            "help|h"        => \$help )
		or die $usage;

sub help {
	printf "%s - %s\n\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub version {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

# Set up some nice output coloring
my $warn = colored("WARN:", 'yellow on_black');
my $info = colored("INFO:", 'cyan on_black');
my $err = colored("ERROR:", 'red on_black');
my $pass = colored("PASS:", 'green on_black');

help if $help;
version if $ver_info;

if ( @ARGV != 1 ) {
		#print $err . "Invalid number of arguments\n\n";
		print "$err Invalid number of arguments\n\n";
        print $usage;
		exit 1;
}

if ( ! $ltable ) {
	print "$err You must supply a lookup table for this analysis.\n\n";
	print $usage;
	exit 1;
}

# Send output to either standard file or to STDOUT for testing purposes
my $out_fh;
if ( $preview ) {
	$out_fh = \*STDOUT;
} else {
	open( $out_fh, ">", $outfile ) || die "$err Can't open the output file '$outfile' for writing: $!";
}


#########------------------------------ END ARG Parsing ---------------------------------#########

# Get the CPSC lookup table
my $scriptdir = dirname(abs_path($0));
my $lookup;

my $lookup_path = "$scriptdir/lookup_tables";
( $ltable eq '3' )   ? $lookup = "${lookup_path}/cpsc_v3_lookupTable.txt"  :
( $ltable eq '4' )   ? $lookup = "${lookup_path}/cpsc_v4_lookupTable.txt"  :
( $ltable eq '4.1' ) ? $lookup = "${lookup_path}/cpsc_v41_lookupTable.txt"  :
#( $ltable eq '4' )   ? $lookup = "${lookup_path}/cpsc_v4_lookupTable.old"  :
( $ltable eq '5' )   ? $lookup = "${lookup_path}/cpsc_v5_lookupTable.txt"  :
( $ltable eq 'mc' )  ? $lookup = "${lookup_path}/cpsc_mc_lookupTable.txt"  :
die "$err '$ltable' is not a valid CPSC lookup table.  Valid options are: '3','4','4.1', '5', or 'mc'.\n";
print colored("The selected lookup table ['$ltable'] is: '$lookup'\n", 'green on_black');

# Load up hash lookup table
open( my $lookup_fh, "<", $lookup ) || die "Lookup file not found: '$!'";
my %plas_lookup; 

# Need to use 0 based position for XLS file and 1 based for VCF
if ($vcf) {
    %plas_lookup = map { my @fields = split; "$fields[1]:$fields[3]" => [@fields] } <$lookup_fh>;
} else {
    %plas_lookup = map { my @fields = split; "$fields[1]:$fields[2]" => [@fields] } <$lookup_fh>;
}
close( $lookup_fh );

#dd \%plas_lookup;
#exit;

my $query = shift;
my ( %results, %counter );

# XXX
# TODO: this sub works. Check that we're getting the input we expect.  Not seeing any indels from the vcfExtractor output 
($vcf) ? vcf_proc( \$query ) : alleles_proc( \$query );

#dd \%counter;
#dd \%results;
#exit;

# Format the fields
#my ( $rwidth, $awidth ) = field_format( \%results );
my $col_width = field_format( \%results, [4,5] );

# Print out the header for the detected plasmids section of the report
print $out_fh "::: Plasmid Variants Detected in Sample :::\n\n";
#printf $out_fh "%-9s %-9s %-13s %-17s %-${rwidth}s %-${awidth}s %-9s %-9s %-13s\n", qw( Gene Chr Position AmpID Ref Alt Freq Cov COSID );
printf $out_fh "%-9s %-9s %-13s %-17s %-$${col_width[0]}s %-$${col_width[1]}s %-9s %-9s %-13s\n", qw( Gene Chr Position AmpID Ref Alt Freq Cov COSID );

# Print out variant call info for found plasmids
#for ( sort { $results{$a}[0] cmp $results{$b}[0] } keys %results ) {
for ( sort { 
        ncmp($results{$a}[0], $results{$b}[0]) ||
        $results{$a}[2] <=> $results{$b}[2]
    } keys %results ) {
	#printf $out_fh "%-9s %-9s %-13s %-17s %-${rwidth}s %-${awidth}s %-9.2f %-9.0f %-13s\n", @{$results{$_}};
	printf $out_fh "%-9s %-9s %-13s %-17s %-$${col_width[0]}s %-$${col_width[1]}s %-9.2f %-9.0f %-13s\n", @{$results{$_}};
}
# Print header for list of plasmids that aren't in the sample
$col_width = field_format( \%plas_lookup, [5,6] );
#print "col width: $$col_width[0]\n";
#exit;

print $out_fh "\n::: Plasmids missing in sample :::\n\n";
#print $out_fh "Gene\tChr\tStart\t\tEnd\t\tCOSID\t\tCDS\t\tSequence\tAmpID\n";
printf $out_fh "%-9s %-9s %-13s %-13s %-13s %-$${col_width[0]}s %-$${col_width[1]}s %-20s\n", qw( Gene Chr Start End COSID CDS Sequence AmpID );

# Print out list of plasmids that were missed
for my $plasmid ( keys %plas_lookup ) {
	if ( ! exists  $counter{$plasmid} ) {
        printf $out_fh "%-9s %-9s %-13s %-13s %-13s %-$${col_width[0]}s %-$${col_width[1]}s %-20s\n", @{$plas_lookup{$plasmid}} ;
		#print $out_fh join( "\t", @{$plas_lookup{$plasmid}} ), "\n";
	}
}

print "\ncpscChecker results written to file '$outfile'\n" if ( ! $preview );
close( $out_fh );

sub field_format {
	# Try to get the largest field width to format the ouput nicely
	my $data = shift;
    my $cols = shift;

    my @format_widths;
    #my $rwidth = 0;
    #my $awidth = 0;

    for my $col (@$cols) {
        my $width = 0;
        for my $elem ( keys %$data ) {
            #print $$data{$elem}->[$field], "\n";
            if ( length $$data{$elem}->[$col] > $width ) {
                $width = ( length $$data{$elem}->[$col] ) + 4;
            }
        }
        push( @format_widths, $width );
    }

    print "formatted fields: \n";
    dd \@format_widths;
	#for my $var ( keys %$data ) {
		#if ( length $$data{$var}->[4] > $rwidth ) {
			#$rwidth = (length $$data{$var}->[4]) + 4;
		#}
		#if ( length $$data{$var}->[5] > $awidth ) {
			#$awidth = (length $$data{$var}->[5]) + 4;
		#}
	#}

	#return ($rwidth, $awidth);
    return( \@format_widths );
}

sub alleles_proc {
    # Process the alleles.xls file and store data in %counter and %results
    my $data_file = shift;

    open ( my $cpsc_fh, "<", $query ) || die "Can not open the file '$query' for reading: $!";
    my $header = <$cpsc_fh>;
    if ( $header !~ /Chrom/ ) {
        print "\n$err input file '$$data_file' does not appear to be an alleles.xls type file.  Check your input file!\n\n";
        print $usage;
        exit 1;
    }
    while (<$cpsc_fh>) {
        #next if ( /Chrom/ || /Absent/ || /No Call/ ); 
        my @line = split;
        my $variant_id = "$line[0]:$line[14]"; # chr:pos
        my $ref = $line[2];
        my $alt = $line[3];

        next if ( ! exists $plas_lookup{$variant_id} );

        my $cosid = $plas_lookup{$variant_id}[4]; # Get COSMIC ID from lookup table as it's more accurate than BED file.

        $results{$variant_id} = [@line[12,0,14,13],$ref,$alt,@line[6,18],$cosid]; 
        $counter{$variant_id} = 1;
    }
}

sub vcf_proc {
    # Process the results of running vcfExtractor -nN <vcf_file> instead of relying on alleles.xls

    my $data_file = shift;

    open ( my $cpsc_fh, "<", $$data_file ) || die "Query table '$$data_file' not found: '$!'\n";
    my $header = <$cpsc_fh>;
    
    # Check to be sure we really have a vcfExtractor output file
    if ( $header !~ /^CHROM:POS/ ) {
        print "\n$err input file '$$data_file' does not appear to be vcfExtractor output. Check your input file!\n\n";
        print $usage;
        exit 1;
    }

    while (<$cpsc_fh>) {
        my @fields = split;
        my $variant_id = $fields[0];
        my $ref = $fields[1];
        my $alt = $fields[2];
        my ($chr, $pos) = split( /:/, $variant_id );

        next if ( ! exists $plas_lookup{$variant_id} );

        my $cosid = $plas_lookup{$variant_id}[4]; # Get COSMIC ID from lookup table as it's more accurate than BED file.
        my $gene = $plas_lookup{$variant_id}[0]; # Have to use lookup for this too; no data in the VCF file
        my $ampid = $plas_lookup{$variant_id}[7]; # Ditto
        
        $results{$variant_id} = [$gene, $chr, $pos, $ampid, $ref, $alt, @fields[5,6], $cosid]; 
        $counter{$variant_id} = 1;
    }
}
