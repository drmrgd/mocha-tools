#!/usr/bin/perl
# Using custom plasmid lookup table for version of CPSC used in experiment, query variants.xls
# file to check to see if the plasmid was seen in the sample and print out the data.  If the
# plasmid was not found, create a new list indicating which plasmids were not found. Print out
# results to a file.
#
# 7/10/2013 
#     - Added capability for CPSCv5
#
# 8/15/2013 
#     - Added capability for MPACT CPSC sample.
#     - Cleaned up and reformatted the code a bit.
#
# 8/30/2013 
#     - Added a fix to recover a COSMIC ID from the original lookup file where there's
#       none present (e.g. '---') due to the HotSpots BED file not having an entry
#       for that particular plasmid.  This will help some downstream analysis.
#
#     - Fixed the lookup tables
#
# 11/15/2013
#     - Started Dev of v2.0 to be compatible with TSv4.0
#
# Created 2/24/2013 - Dave Sims
#
###################################################################################################		

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;

( my $scriptname = $0 ) =~ s/^(.*\/)+//;
my $version = "v2.0.0";
my $description = <<"EOT";

Using a plasmid lookup table for the version of the CPSC used in the experiment, query
a TVC variants.xls file to check to see if the plasmid was seen in the sample, and print
out the data.  If the plasmid was not found, create a new list indicating which plasmids
were not found.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] -l <3,4,5,mc>  <variant_call_file>

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
my $lookup;

# TODO: use relative path
my $lookup_path = "/opt/mocha/cpscChecker/lookup_tables";
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

# load up query table and print matching records
my $query = shift;
open ( QUERY, "<", $query ) || die "Query table '$query' not found: '$!'\n";
my ( %results, %counter );
while (<QUERY>) {
	chomp;
	next if ( /Gene Sym/ || /Ref/ || /NC/ || /---/ );
	my @line = split;
	my $variant_id = "$line[0]:$line[1]";
	
	next if ( ! exists $plas_lookup{$variant_id} );

	my ($cosid) = grep { $line[14] =~ /($_)/ } @{$plas_lookup{$variant_id}}; 

	# Use Genotype field to create REF and ALT data
	my ( $ref, $alt ) = $line[6] =~ /^(\w+)\/(\w+)$/;

	$results{$variant_id} = [@line[2,0,1,3],$ref,$alt,@line[9,11],$cosid]; 
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
