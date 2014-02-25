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
use Getopt::Long;
use Data::Dumper;
use Data::Dump;
use Cwd 'abs_path';
use File::Basename;

( my $scriptname = $0 ) =~ s/^(.*\/)+//;
my $version = "v2.0.2";
my $description = <<"EOT";

Using a plasmid lookup table for the version of the CPSC used in the experiment, query
a TVC 'alleles.xls' file to check to see if the plasmid was seen in the sample, and print
out the data.  If the plasmid was not found, create a new list indicating which plasmids
were not found.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] -l <3,4,5,mc>  <alleles.xls>

    -l, --lookup           CPSC lookup version to use (3, 4, 5, or mc)
    -o, --output <file>    Output file name (default is 'output.txt')
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

GetOptions( "lookup=s"    => \$ltable,
	        "preview"     => \$preview,
            "output=s"    => \$outfile,
            "version"     => \$ver_info,
            "help"        => \$help )
		or die $usage;

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

if ( @ARGV != 1 ) {
		printf "%s: Invalid number of arguments\n\n", $scriptname;
		printf "%s\n", $usage;
		exit 1;
}

if ( ! $ltable ) {
	print "ERROR: You must supply a lookup table for this analysis.\n\n";
	print $usage;
	exit 1;
}

# Send output to either standard file or to STDOUT for testing purposes
my $out_fh;
if ( $preview ) {
	$out_fh = \*STDOUT;
} else {
	open( $out_fh, ">", $outfile ) || die "Can't open the output file '$outfile' for writing: $!";
}

#########------------------------------ END ARG Parsing ---------------------------------#########

# Get the CPSC lookup table
my $scriptdir = dirname(abs_path($0));
my $lookup;

my $lookup_path = "$scriptdir/lookup_tables";
( $ltable eq '3' )   ? $lookup = "${lookup_path}/cpsc_v3_lookupTable.txt"  :
( $ltable eq '4' )   ? $lookup = "${lookup_path}/cpsc_v4_lookupTable.txt"  :
( $ltable eq '5' )   ? $lookup = "${lookup_path}/cpsc_v5_lookupTable.txt"  :
( $ltable eq 'mc' )  ? $lookup = "${lookup_path}/cpsc_mc_lookupTable.txt"  :
die "ERROR: '$ltable' is not a valid CPSC lookup table.  Valid options are: '3','4','5', or 'mc'.\n";
print "The selected lookup table ['$ltable'] is: '$lookup'\n";

# Load up hash lookup table
open( my $lookup_fh, "<", $lookup ) || die "Lookup file not found: '$!'";
my %plas_lookup = map { my @fields = split; "$fields[1]:$fields[2]" => [@fields] } <$lookup_fh>;
close( $lookup_fh );

# load up query table and print matching records; have to use alleles.xls files now
my $query = shift;
open ( my $cpsc_fh, "<", $query ) || die "Query table '$query' not found: '$!'\n";
my ( %results, %counter );
while (<$cpsc_fh>) {
	chomp;
	next if ( /Chrom/ || /Absent/ || /No Call/ ); 
	my @line = split;
	my $variant_id = "$line[0]:$line[14]";
    my $ref = $line[2];
    my $alt = $line[3];

	next if ( ! exists $plas_lookup{$variant_id} );

	#my ($cosid) = grep { $line[11] =~ /($_)/ } @{$plas_lookup{$variant_id}}; 
    my $cosid = $plas_lookup{$variant_id}[4]; # Get COSMIC ID from lookup table as it's more accurate than BED file.

	# Use Genotype field to create REF and ALT data
	#my ( $ref, $alt ) = $line[6] =~ /^(\w+)\/(\w+)$/;
    
	$results{$variant_id} = [@line[12,0,14,13],$ref,$alt,@line[6,18],$cosid]; 
	$counter{$variant_id} = 1;
}

# Format the fields
my ( $rwidth, $awidth ) = field_format( \%results );

# Print out the header for the detected plasmids section of the report
print $out_fh "::: Plasmid Variants Detected in Sample :::\n\n";
printf $out_fh "%-9s %-9s %-13s %-17s %-${rwidth}s %-${awidth}s %-9s %-9s %-13s\n", qw( Gene Chr Position AmpID Ref Alt Freq Cov COSID );

# Print out variant call info for found plasmids
for ( sort { $results{$a}[0] cmp $results{$b}[0] } keys  %results ) {
	printf $out_fh "%-9s %-9s %-13s %-17s %-${rwidth}s %-${awidth}s %-9.2f %-9.0f %-13s\n", @{$results{$_}};
}

# Print header for list of plasmids that aren't in the sample
print $out_fh "\n::: Plasmids missing in sample :::\n\n";
print $out_fh "Gene\tChr\tStart\t\tEnd\t\tCOSID\t\tCDS\t\tSequence\tAmpID\n";

# Print out list of plasmids that were missed
for my $plasmid ( keys %plas_lookup ) {
	if ( ! exists  $counter{$plasmid} ) {
		print $out_fh join( "\t", @{$plas_lookup{$plasmid}} ), "\n";
	}
}

print "\ncpscChecker results written to file '$outfile'\n" if ( ! $preview );
close( $out_fh );

sub field_format {
	# Try to get the largest field width to format the ouput nicely
	my $data = shift;
	my $rwidth = 0;
	my $awidth = 0;

	for my $var ( keys %$data ) {
		if ( length $$data{$var}->[4] > $rwidth ) {
			$rwidth = (length $$data{$var}->[4]) + 4;
		}
		if ( length $$data{$var}->[5] > $awidth ) {
			$awidth = (length $$data{$var}->[5]) + 4;
		}
	}

	return ($rwidth, $awidth);
}
