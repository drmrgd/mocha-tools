#!/usr/bin/perl
# Program to create a sampleKey.txt file for Ion Torrent runs based on the information entered into the Plan
# Template during setup.
#
# 11/7/2013 - updated for TSv4.0
#
# Created: 3/3/2013 - Dave Sims
###############################################################################################################
use warnings;
use strict;
use JSON;
use Getopt::Long;
use Data::Dump;

( my $scriptname = $0 ) =~ s/^(.*\/)+//;
my $version = "v2.0.0";
my $description = <<"EOT";
Program to read in the datasets_pipeline.json file setup during the planning for an Ion Torrent run, and 
create a tab delimited sampleKey.txt file used in the rest of the variant reporting downstream. 

If using the default JSON file, this script should be run from an experiment results directory.  If running 
from outside the results directory, a datasets_pipeline.json file can be passed to the script with the '-f'
option.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] 
    -f, --file        Custom JSON file if not running from default path
    -t, --test        Write output to STDOUT rather than 'sampleKey.txt
    -v, --version     Version Information
    -h, --help        Display the help information
EOT

my $ver_info;
my $help;
my $preview;
my $outfile = "sampleKey.txt";
my $jsonfile = "basecaller_results/datasets_pipeline.json";

GetOptions( "test"       => \$preview,
	        "file=s"     => \$jsonfile,
	        "version"    => \$ver_info,
			"help"       => \$help )
			or print $usage;

sub help {
	printf "%s - %s\n\n%s\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub version {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

help if $help;
version if $ver_info;

# Create an output fh; either STDOUT for testing or sampleKey.txt for production
my $out_fh;
if ( ! $preview ) {
	open( $out_fh, ">", $outfile ) || die "Can't write the file '$outfile': $!\n";
} else {
	$out_fh = \*STDOUT;
}

if ( ! -f $jsonfile ) {
	print "ERROR: The JSON file selected '$jsonfile' can not be found.  Check your path.\n";
	exit 1;
}

##----------------------------------------- End Command Arg Parsing ------------------------------------##

# Read in the JSON file and parse it for the sample and barcode info
my %bc_hash;
open( my $json_fh, "<", $jsonfile ) || die "Can't open file for reading. $!\n";
my $parsed_json;
{
	local $/;
	my $json = JSON->new;
	eval { $parsed_json = $json->decode(<$json_fh>) } 
}
close( $json_fh );

# Sanity check to make sure the JSON file is valid
if ( ! $$parsed_json{'barcode_config'} || $$parsed_json{'barcode_config'}->{'barcode_id'} ne 'IonXpress' ) {
	print "ERROR: The JSON file loaded '$jsonfile' does not contain information about the IonXpress barcodes.\n\n";
	print $usage;
	exit 1;
}


for my $read_bc ( keys %{$parsed_json->{'read_groups'}} ) {
	next if ( $read_bc =~ /nomatch/ );
	my $barcode = $$parsed_json{'read_groups'}->{$read_bc}->{'barcode_name'};
    my $sample = $$parsed_json{'read_groups'}->{$read_bc}->{'sample'};
	$bc_hash{$barcode} = $sample;
}

# Create a sampleKey file from the final hash
foreach my $bc ( sort  keys %bc_hash ) {
	print $out_fh "$bc\t$bc_hash{$bc}\n" unless ( $bc_hash{$bc} eq 'None' );
}
close( $out_fh );
