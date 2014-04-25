#!/usr/bin/perl
# Program to create a sampleKey.txt file for Ion Torrent runs based on the information entered into the Plan
# Template during setup.  The program reads the ion_params_00.json file for the experiment, parses the file 
# for the barcode index, and cross-references that list with the 'IonXpress_xxx' nomenclature. Requires a 
# barcode lookup table.
#
# Created: 3/3/2013 - Dave Sims
#
###############################################################################################################
use warnings;
use strict;
#use JSON; 
use JSON::XS;
use Data::Dumper;
use feature qw{ switch };

( my $scriptname = $0 ) =~ s/^(.*\/)+//;
my $version = "v0.9.1";
my $description = <<"EOT";
Program to read in the ion_params_00.json file setup during the planning for an Ion Torrent run, and create
a sampleKey.txt file used in the rest of the variant reporting downstream.  The JSON file describes the 
barcodes used during the setup, and this script it necessary to associate the two as a different index is 
used by the system than the 'IonXpress_xxx' nomenclature used elsewhere.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <bc_lookup_table> <ion_params_00.json>
	-o, --output		Non-standard output filename
	-v, --version		Version Information
	-h, --help		Display the help information
EOT

sub help {
	printf "%s - %s\n\n%s\n%s\n", $scriptname, $version, $description, $usage;
}

my $verInfo = 0;
my $help = 0;
my $outfile = "sampleKey.txt";

while ( scalar( @ARGV ) > 0 ) {
	last if ( $ARGV[0] !~ /^-/ );
	my $opt = shift;
	if ( $opt eq '-h' || $opt eq '--help' ) { $help = 1; } 
	elsif ( $opt eq '-v' || $opt eq '--version' ) { $verInfo = 1; }
	elsif ( $opt eq '-o' || $opt eq '--output' ) { $outfile = shift; }
	else {
		print "$scriptname: Invalid option: $opt\n\n";
		print "$usage\n";
		exit 1;
	}
}

if ( $verInfo ) {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}
elsif ( $help ) {
	help;
	exit;
}

if ( scalar( @ARGV ) != 2 ) {
	print "ERROR: Not enough arguments passed to script\n\n";
	print "$usage\n";
	exit;
}

my $bctable_in = shift;
my $jsonfile = shift;

##----------------------------------------- End Command Arg Parsing ------------------------------------##

# Slurp in the JSON file and parse it for the barcodedSamples line
open( my $input, "<", $jsonfile ) || die "Can't open file for reading. $!\n";
my $json_in = <$input>;
my $parsed_json = decode_json( $json_in );
my $barcodedSamples = $parsed_json->{'plan'}->{'barcodedSamples'};
my @json_lines = split( /,/, $barcodedSamples );
my $instrument = $parsed_json->{'pgmName'};
my ( $run_num ) = $parsed_json->{'resultsName'} =~ /.*user_([PM]CC-\d+).*/;

if ( $instrument eq 'uploads' ) {
    print "Detected imported run.  Using previous index information from $run_num\n";
}

close ( $input );

# Load in the BC Lookup Table
my %bclist;
open( my $bctable, "<", $bctable_in  ) || die "Can't open the lookup file. $!\n";
while ( <$bctable> ) {
	chomp;
	my @line = split( /\t/, $_ );

	#$bclist{$line[0]} = $line[1] if ( $instrument =~ /MCC-CLIA/ ); # Standard for CLIA TS
	#$bclist{$line[0]+16} = $line[1] if ( $instrument =~ /PCCPGM/ ); # Add offset of 16 to account R&D TS difference.
    given ( $instrument ) {
        when ( /MCC-CLIA/ ) { $bclist{$line[0]} = $line[1] }
        when ( /PCCPGM/ )    { $bclist{$line[0]+16} = $line[1] }
        when ( /MCC-DEV/ )   { $bclist{$line[0]} = $line[1] }
        when ( /uploads/ )   { ( $run_num =~ /PCC/ ) ? $bclist{$line[0]+16} = $line[1] : $bclist{$line[0]} = $line[1];}
        default { die 'No valid instrument found in barcode lookup list.' }
    }
}
close ( $bctable );

# Extract the bc index and sample information and load it into a new hash
my %bcSamples;
foreach ( @json_lines ) {
	( my $elem = $_ ) =~ s/\{|\}|\"//g;
	( my @sampleline ) = split( /:/, $elem );
	#print "$elem\n";
	$bcSamples{$sampleline[0]} = $sampleline[2];
}

# Cross-reference the bclist hash with the bcsamples hash to create a summaryKey.txt file
my %sampleKey = map { exists $bcSamples{$_} ? ( $bclist{$_} => $bcSamples{$_} ) : () } keys %bclist;
#print Dumper( \%sampleKey );
#exit;

# Create a sampleKey file from the final hash
open( my $fh_out, ">", $outfile ) || die "Can't write the file 'sampleKey.txt'. $1\n";
foreach my $bc ( sort  keys %sampleKey ) {
	
	my @line = join ("\t", $bc, $sampleKey{$bc} );
	print $fh_out "$_\n" for @line;
}
close( $fh_out );
