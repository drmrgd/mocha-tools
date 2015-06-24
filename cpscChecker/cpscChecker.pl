#!/usr/bin/perl
# Generate a CPSC Checker report starting with a VCF file. This replaces the original CPSC Checker
# script which required a variants.xls file to run.  
#
# 5/28/2015 - D Sims
######################################################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Cwd qw(abs_path getcwd);
use Term::ANSIColor;
use Data::Dump;
use Sort::Versions;

# Remove when in prod.
print "\n";
print colored("*" x 50, 'bold yellow on_black');
print colored("\n    DEVELOPMENT VERSION OF CPSC CHECKER\n", 'bold yellow on_black');
print colored("*" x 50, 'bold yellow on_black');
print "\n\n";

my $scriptname = basename($0);
my $version = "v3.6.0_062415-dev";
my $description = <<"EOT";
Using a plasmid lookup table for the version of the CPSC used in the experiment, query a TVC VCF
file to check to see if the plasmids were seen in the sample, and print out the data.  If the 
plasmid was not found, print out a list of missing plasmids.  
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] -l <lookup_file> -V <vcf_file> | -T <tabular file>
    -l, --lookup           CPSC lookup file to use.  Use '?' to see list of available files.
    -c, --custom_lookup    Use a custom lookup file rather than one of the built-in ones.
    -V, --VCF              Input is from a TVC VCF File
    -T, --Tab              Input is from a TVC tab delimited file (i.e. alleles.xls).
    -o, --output           Send output to custom file.  Default is STDOUT.
    -v, --version          Version information
    -h, --help             Print this help information
EOT

my $help;
my $ver_info;
my $outfile;
my $lookup;
my $custom_lookup;
my $vcf;
my $tab;

GetOptions( "lookup|l=s"         => \$lookup,
            "custom_lookup|c=s"  => \$custom_lookup,
            "VCF|V=s"            => \$vcf,
            "Tab|T=s"            => \$tab,
            "output|o=s"         => \$outfile,
            "version|v"          => \$ver_info,
            "help|h"             => \$help )
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

# Make sure enough args passed to script
unless ($tab || $vcf) {
    print "ERROR: Not enough arguments passed to script!\n\n";
    print "$usage\n";
    exit 1;
}

#if ( scalar( @ARGV ) < 1 ) {
    #print "ERROR: Not enough arguments passed to script!\n\n";
    #print "$usage\n";
    #exit 1;
#}

# Set up some formatted and colored output
my $err = colored( 'ERROR:', 'bold red on_black');
my $warn = colored( 'WARNING:', 'bold yellow on_black');
my $info = colored( 'INFO:', 'bold cyan on_black');

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
	open( $out_fh, ">", $outfile ) || die "Can't open the output file '$outfile' for writing: $!";
} else {
	$out_fh = \*STDOUT;
}

#########------------------------------ END ARG Parsing ---------------------------------#########
my $cpsc_lookup;
if ($lookup) {
    $cpsc_lookup = validate_lookup(\$lookup);
}
elsif ($custom_lookup) {
    if ( ! -e $custom_lookup ) {
        print "$err Selected lookup file '$custom_lookup' does not exist!\n";
        validate_lookup(\'?');
        exit 1;
    }
    $cpsc_lookup = validate_lookup(\$custom_lookup);
} else {
    print "$err You must load a lookup file with either the '--lookup' or '--custom_lookup' options\n";
    exit 1;
}

# Generate lookup dataset
my %cpsc_lookup_data = load_lookup(\$cpsc_lookup);

# read variant file and store data
my %plas_results;
my %raw_data;
if ($vcf) {
    if ( ! $vcf && ! -e $vcf ) {
        die "$err Can't read VCF file '$vcf': $!";
    } else {
        %raw_data = read_vcf(\$vcf);
    }
} else {
    my %tab_data;
    if ( ! $tab && ! -e $tab ) {
        die "$err Can't read alleles tab file: '$tab': $!";
    } else {
        %raw_data = read_tab(\$tab);
    }
}

# Check VCF dataset against lookup hash for final results.
%plas_results = proc_plas_data(\%cpsc_lookup_data,\%raw_data);

#dd \%plas_results;
#exit;

# Print out the results.
print_results(\%plas_results, \%cpsc_lookup_data);

sub read_tab {
    my $input_file = shift;
    my %tab_data;

    print "Checking tab delimited variants file '$$input_file' for CPSC data...\n";

    open( my $tab_fh, "<", $$input_file );
    while (<$tab_fh>) {
        chomp;
        next unless /^chr/;
        my @fields = split(/\t/);
        my $var_coord = join(':', @fields[1,2]);
        my $ref_cov = $fields[18]-$fields[24];
        $tab_data{join(':', @fields[0..3])} = [@fields[0..6],$fields[18],$ref_cov,$fields[24],$fields[11]];
    }
    close $tab_fh;

    return %tab_data;
}

sub read_vcf {
    my $vcf_file = shift;
    my %vcf_data;

    print "Checking VCF file '$$vcf_file' for CPSC data...\n";

    chomp(my $path = qx( which vcfExtractor.pl ));
    my $vcf_extractor_cmd = qq{ $path -Nn $$vcf_file };
    open( my $parsed_vcf, "-|", $vcf_extractor_cmd ); 
    while (<$parsed_vcf>) {
        next unless /^chr/;
        my @fields = split;
        #my $varid = join( ':', @fields[0..2] );
        $vcf_data{join(':', @fields[0..2])} = [@fields];
    } 
    return %vcf_data;
}

sub validate_lookup {
    my $file_ref = shift;
    my $lfile_path = dirname(abs_path($0));
    my %lfiles = (
        '3'    => 'cpsc_v3_lookupTable.txt',
        '4'    => 'cpsc_v4_lookupTable.txt',
        '4.1'  => 'cpsc_v41_lookupTable.txt',
        '4ocp' => 'cpsc_v4ocp_lookupTable.txt',
        'mc'   => 'cpsc_mc_lookupTable.txt',
        '5'    => 'cpsc_v5_lookupTable.txt',
        'sc1'  => 'seracare_misc_v1.0_022515.txt',
        'sc2'  => 'seracare_v2.0_060115.txt',
    );
    my $lookup_file;

    if ( $$file_ref eq '?' ) {
        print "Valid lookup files are: \n";
        printf "\t%-8s=>  %s\n", ($_, $lfiles{$_}) for sort keys %lfiles;
        exit 0;
    } 
    elsif ( ! exists $lfiles{$$file_ref} ) {
        print "$err Lookup table '$$file_ref' is not valid!\n";
        validate_lookup(\'?');
        exit 1;
    } else {
        $lookup_file = "$lfile_path/lookup_tables/$lfiles{$$file_ref}";
    }

    print "$info Lookup file '$lookup_file' selected\n";
    return $lookup_file;
}

sub load_lookup {
    my $file = shift;
    my %cpsc_data;

    open( my $lookup_fh, "<", $$file );
    while (<$lookup_fh>) {
        next if /^\s*$/;
        my @fields = split;
        $cpsc_data{join(':', @fields[1,2,5,6])} = [@fields];
    }
    close $lookup_fh;
    return %cpsc_data;
}

sub proc_plas_data {
    my ($cpsc_lookup,$vcf_data) = @_;
    my %results;
    my %missing;

    for my $var (keys %$vcf_data) {
        next unless (exists $$cpsc_lookup{$var});
        @{$results{$var}} = @{$$vcf_data{$var}};

        # Do a little post processing to get the gene name and fix variant ID
        push( @{$results{$var}}, $$cpsc_lookup{$var}->[0] );
        if ($results{$var}->[9] eq '.') {
            $results{$var}->[9] = $$cpsc_lookup{$var}->[3];
        }
    }
    #dd \%results;
    #exit;

    return %results;
}

sub print_results {
    my $found_plasmids = shift;
    my $lookup = shift;

    select $out_fh;
    my @header_elems = qw( Position Gene REF ALT VAF TotCov RefCov AltCov VARID );
    my $col_widths = field_width($found_plasmids, [1,2]);
    my $format = "%-20s%-14s%-$${col_widths[0]}s%-$${col_widths[1]}s%-10s%-8s%-8s%-8s%-10s\n";

    print ":::  Plasmids Detected in Sample  :::\n";
    printf $format, @header_elems;
    for my $var ( sort { versioncmp( $a, $b ) } keys %$found_plasmids ) {
        printf $format, @{$$found_plasmids{$var}}[0,10,1,2,5,6,7,8,9];
    }

    # Print out a list of thos missing only if we have some missing.
    if ( (scalar keys %$lookup) > (scalar keys %$found_plasmids) ) {
        print "\n:::  Plasmids Missing in Sample  :::\n";
        $col_widths = field_width( $lookup, [4,5] );
        $format = "%-8s%-9s%-14s%-15s%-$${col_widths[0]}s%-$${col_widths[1]}s\n";
        printf $format, qw(Gene Chr Position VARID CDS Sequence);
        for my $plas (keys %$lookup) {
            printf $format, @{$$lookup{$plas}} if (! exists $$found_plasmids{$plas});
        }
    }

    return;
}

sub field_width {
    my $data = shift;
    my $cols = shift;

    my @widths;

    for my $col (@$cols) {
        my $width = 0;
        for my $elem (keys %$data) {
            my $length = length $$data{$elem}->[$col];
            if ( $length > $width ) {
                $width = ($length);
            }
        }
        push( @widths, $width );
    }
    # Add some extra padding in there.
    @widths = map { $_ + 4 } @widths;
    return \@widths;
}
