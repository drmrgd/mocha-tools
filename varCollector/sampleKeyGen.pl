#!/usr/bin/env perl
# Program to create a sampleKey.txt file for Ion Torrent runs based on the 
# information entered into the Plan Template during setup.
#
# Created: 3/3/2013 - Dave Sims
################################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long qw(:config ignore_case bundling no_auto_abbrev);
use Data::Dump;
use JSON;
use File::Basename;

my $scriptname = basename($0);
my $version = "v3.5.012616";
my $description = <<"EOT";
Program to read in the ion_params_00.json file setup during the planning for
an Ion Torrent run, and create a tab delimited sampleKey.txt file used in the 
rest of the variant reporting downstream. 

If using the default JSON file, this script should be run from an experiment 
results directory.  If running from outside the results directory, a 
datasets_pipeline.json file can be passed to the script with the '-f' option. 
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] 
    -f, --file        Custom JSON file if not running from default path
    -p, --preview     Write output to STDOUT rather than 'sampleKey.txt
    -o, --output      Write output to custom file (DEFAULT: 'sampleKey.txt')
    -v, --version     Version Information
    -h, --help        Display the help information
EOT

my $ver_info;
my $help;
my $preview;
my $outfile = "sampleKey.txt";
my $jsonfile = 'ion_params_00.json';
my $custom_json;

GetOptions( "preview|p"    => \$preview,
	        "file|f=s"     => \$custom_json,
            "output|o=s"   => \$outfile,
	        "version|v"    => \$ver_info,
			"help|h"       => \$help )
			or print $usage;

sub help {
	printf "%s - %s\n\n%s\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub version {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

help() if $help;
version() if $ver_info;

# Create an output fh; either STDOUT for testing or sampleKey.txt for production
my $out_fh;
if ( ! $preview ) {
	open( $out_fh, ">", $outfile ) || die "Can't write the file '$outfile': $!\n";
} else {
	$out_fh = \*STDOUT;
}

# Allow for custom JSON if we want.
$jsonfile = $custom_json if $custom_json; 
if ( ! -f $jsonfile) {
    die "ERROR: The JSON file '$jsonfile' can not be found. Check your path!\n";
}

##----------------------------------------- End Command Arg Parsing ------------------------------------##
my %bc_hash;

# Create hash map of data.
map_samples(\$jsonfile);

# Make sure that there is data to be printed, otherwise have to build a sample key manually
if ( ( grep { $bc_hash{$_} =~ /none/i } keys %bc_hash ) == scalar( keys %bc_hash ) ) {
    print "ERROR: No samples listed for any of the barcodes in the experiment ",
        "JSON file.  A run plan was possibly not uploaded prior to running.\n";
    print "You will need to manually create a sampleKey for this experiment\n";
    exit 1;
}

# Print out the final sampleKey data
foreach my $bc ( sort keys %bc_hash ) {
	print {$out_fh} "$bc\t$bc_hash{$bc}\n" unless ( $bc_hash{$bc} =~ /none/i );
}
close $out_fh;

sub read_json {
    # Read in the Ion JSON file
    my $json_file = shift;
    my $parsed_json;
    local $/;

    open( my $json_fh, "<", $$json_file );
    my $json = JSON->new;
    eval { $parsed_json = $json->decode(<$json_fh>) };
    close $json_fh;

    return $parsed_json;
}

sub map_samples {
    my $jsonfile = shift;
    my $parsed_json = read_json($jsonfile);
    my $sample_table = $parsed_json->{'experimentAnalysisSettings'}{'barcodedSamples'};

    while (my ($sample, $barcode_info) = each %$sample_table) {
        for my $barcode (keys %{$barcode_info->{'barcodeSampleInfo'}}) {
            my $na_type = $barcode_info->{'barcodeSampleInfo'}{$barcode}{'nucleotideType'};
            $bc_hash{$barcode} = "$sample-$na_type";
        }
    }
    return;
}
