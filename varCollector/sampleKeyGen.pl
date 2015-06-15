#!/usr/bin/perl
# Program to create a sampleKey.txt file for Ion Torrent runs based on the information entered into the Plan
# Template during setup.
#
# Created: 3/3/2013 - Dave Sims
###############################################################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long qw(:config ignore_case bundling no_auto_abbrev);
use Data::Dump;
use JSON;
use File::Basename;

my $scriptname = basename($0);
my $version = "v3.0.0_061515";
my $description = <<"EOT";
Program to read in the datasets_pipeline.json file setup during the planning for an Ion Torrent run, and 
create a tab delimited sampleKey.txt file used in the rest of the variant reporting downstream. 

If using the default JSON file, this script should be run from an experiment results directory.  If running 
from outside the results directory, a datasets_pipeline.json file can be passed to the script with the '-f'
option. This version will also process pre-v4.0 data by using the old ion_params_00.json file.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] 
    -f, --file        Custom JSON file if not running from default path
    -r, --retro       Use the 'ion_params_00.json' file instead for older analysis versions (<TSSv4.0) 
    -p, --preview     Write output to STDOUT rather than 'sampleKey.txt
    -o, --output      Write output to custom file (DEFAULT: 'sampleKey.txt')
    -v, --version     Version Information
    -h, --help        Display the help information
EOT

my $ver_info;
my $help;
my $preview;
my $outfile = "sampleKey.txt";
my $jsonfile;
my $custom_json;
my $retro_data;


GetOptions( "preview|p"    => \$preview,
            "retro|r"      => \$retro_data,
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

help if $help;
version if $ver_info;

# Create an output fh; either STDOUT for testing or sampleKey.txt for production
my $out_fh;
if ( ! $preview ) {
	open( $out_fh, ">", $outfile ) || die "Can't write the file '$outfile': $!\n";
} else {
	$out_fh = \*STDOUT;
}

# Allow for custom JSON file processing.
if ($retro_data) {
    if ($custom_json) {
        $jsonfile = $custom_json;
    } else {
        $jsonfile = 'ion_params_00.json';
    }
} else {
    if ($custom_json) {
        $jsonfile = $custom_json;
    } else {
        $jsonfile = "basecaller_results/datasets_pipeline.json";
    }
}

die "ERROR: The JSON file '$jsonfile' can not be found. Check your path!" unless (-f $jsonfile);

##----------------------------------------- End Command Arg Parsing ------------------------------------##
my %bc_hash;
($retro_data) ? (%bc_hash = map_samples_retro(\$jsonfile)) : (%bc_hash = map_samples(\$jsonfile));

# Make sure that there is data to be printed, otherwise have to build a sample key manually
if ( ( grep { $bc_hash{$_} =~ /none/i } keys %bc_hash ) == scalar( keys %bc_hash ) ) {
    print "ERROR: No samples listed for any of the barcodes in the experiment JSON file.  A run plan was possibly not uploaded prior to running.\n";
    print "You will need to manually create a sampleKey for this experiment\n";
    exit 1;
}

# Print out the final sampleKey data
foreach my $bc ( sort keys %bc_hash ) {
	print {$out_fh} "$bc\t$bc_hash{$bc}\n" unless ( $bc_hash{$bc} eq 'None' );
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

sub map_samples_retro {
    # Used for pre-v4.0 TSS versions, which mapped the ion_params_00.json file
    my $jsonfile = shift;
    
    #my $bc_lookup = ion_params_decoder();
    my $parsed_json = read_json($jsonfile);

    my %bc_lookup;
    my $shift;
    # For some reason the PCCPGM file is different than others.  Have to shift the index.
    ($parsed_json->{'pgmName'} eq 'PCCPGM') ? ($shift = 16) : ($shift = 0);
    for (my $index = 1; $index <= 96; $index++ ) {
        my $bc_number = sprintf("%03d", $index);
        $bc_lookup{$shift+$index+32} = "IonXpress_$bc_number";
    }

    # Extract the sample name part of the JSON and prepare for mapping.
    my @json_sample_info = split( /,/, $parsed_json->{'plan'}->{'barcodedSamples'});
    @json_sample_info = map { $_ =~ s/\{|\}|\"//g; $_ } @json_sample_info;

    # Map the samples to the correct barcode index
    foreach (@json_sample_info) {
        my ($index, undef, $sample_name) = split(/:/, $_);
        $bc_hash{$bc_lookup{$index}} = $sample_name;
    }
    return %bc_hash;
}

sub map_samples {
    # Use the datasets_pipeline.json file to map samples to the barcode used based on the run plan.
    my $jsonfile = shift;
    my $parsed_json = read_json($jsonfile);

    # Sanity check to make sure the JSON file is valid
    if ( ! $$parsed_json{'barcode_config'} || $$parsed_json{'barcode_config'}->{'barcode_id'} ne 'IonXpress' ) {
        print "ERROR: The JSON file loaded '$jsonfile' does not contain information about IonXpress barcodes.\n\n";
        print $usage;
        exit 1;
    }

    for my $read_bc ( keys %{$parsed_json->{'read_groups'}} ) {
        next if ( $read_bc =~ /nomatch/ );
        my $barcode = $$parsed_json{'read_groups'}->{$read_bc}->{'barcode_name'};
        my $sample = $$parsed_json{'read_groups'}->{$read_bc}->{'sample'};
        $bc_hash{$barcode} = $sample;
    }
    return %bc_hash;
}
