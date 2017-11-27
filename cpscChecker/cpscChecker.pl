#!/usr/bin/perl
# Generate a CPSC Checker report starting with a VCF file. This replaces the original CPSC Checker
# script which required a variants.xls file to run.  
#
# 5/28/2015 - D Sims
######################################################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long qw(:config bundling auto_abbrev no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Term::ANSIColor;
use Data::Dump;
use Sort::Versions;

my $version = "v4.2.112717";

print "\n";
print colored("*" x 75, 'bold yellow on_black');
print colored("\n    DEVELOPMENT VERSION OF CPSC CHECKER (version: $version)\n", 'bold yellow on_black');
print colored("*" x 75, 'bold yellow on_black');
print "\n\n";

my $scriptname = basename($0);
my $description = <<"EOT";
Using a plasmid lookup table for the version of the CPSC used in the experiment, query a TVC VCF
file to check to see if the plasmids were seen in the sample, and print out the data.  If the 
plasmid was not found, print out a list of missing plasmids.  
EOT

my $help;
my $ver_info;
my $outfile;
my $lookup;
my $custom_lookup;
my $assay_type = 'solid';

my $usage = <<"EOT";
USAGE: $scriptname [options] -l <lookup_file> <vcf_file>
    -l, --lookup           CPSC lookup file to use.  Use '?' to see list of available files.
    -c, --custom_lookup    Use a custom lookup file rather than one of the built-in ones.
    -a, --assay            Assay type. Choices are 'cfdna' or 'solid' (DEFAULT: $assay_type).
    -o, --output           Send output to custom file.  Default is STDOUT.
    -v, --version          Version information
    -h, --help             Print this help information
EOT

GetOptions( 
    "lookup|l=s"         => \$lookup,
    "custom_lookup|c=s"  => \$custom_lookup,
    "assay|a=s"          => \$assay_type,
    "output|o=s"         => \$outfile,
    "version|v"          => \$ver_info,
    "help|h"             => \$help 
) or die $usage;

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

# Set up some formatted and colored output
my $err = colored( 'ERROR:', 'bold red on_black');
my $warn = colored( 'WARNING:', 'bold yellow on_black');
my $info = colored( 'INFO:', 'bold cyan on_black');

# Check to make sure that we have vcfExtractor in our path.
my $vcf_extractor_cmd = __check_prog($assay_type);

# Make sure enough args passed to script
my $vcf = shift;
if ( defined $lookup && $lookup ne '?' ) {
    unless ($vcf) {
        print "$err You must load a VCF file!\n\n";
        print "$usage\n";
        exit 1;
    }
}

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
	open( $out_fh, ">", $outfile ) || die "Can't open the output file '$outfile' for writing: $!";
} else {
	$out_fh = \*STDOUT;
}

#########------------------------------ END ARG Parsing ---------------------------------#########

# Generate lookup dataset
my %cpsc_lookup_data;
if ($lookup) {
    my $cpsc_lookup = validate_lookup(\$lookup);
    %cpsc_lookup_data = load_lookup(\$cpsc_lookup);
}
elsif ($custom_lookup) {
    if ( ! -e $custom_lookup ) {
        print "$err Selected lookup file '$custom_lookup' does not exist!\n";
        validate_lookup(\'?');
        exit 1;
    }
    %cpsc_lookup_data = load_lookup(\$custom_lookup);
} else {
    print "$err You must load a lookup file with either the '--lookup' or '--custom_lookup' options\n";
    exit 1;
}

# read variant file and store data
my %raw_data = read_vcf(\$vcf, $vcf_extractor_cmd);

dd \%raw_data;
exit;

# Check VCF dataset against lookup hash for final results.
my %plas_results = proc_plas_data(\%cpsc_lookup_data, \%raw_data);

# Print out the results.
print_results(\%plas_results, \%cpsc_lookup_data);

sub read_vcf {
    my ($vcf_file, $cmd) = @_;
    my %vcf_data;

    print "Checking VCF file '$$vcf_file' for CPSC data...\n";

    $cmd .= " $$vcf_file";

    my ($parsed_vcf, $child_pid, $child_rc);
    # Better error handling to capture failed VCF Extractor output.
    unless($child_pid = open($parsed_vcf, "-|")) {
        open(STDERR, ">&STDOUT");
        exec($cmd);
        die "ERROR: could not execute program: $!";
    }
    waitpid($child_pid,0);
    $child_rc = $? >> 8;
    if ($child_rc) {
        my ($err_msg) = grep { $_ =~ /ERROR/ } <$parsed_vcf>;
        print "$err Can not parse VCF file '$$vcf_file'! Message from calling prog:\n";
        print "  ->  " . $err_msg;
        exit 1;
    }

    while (<$parsed_vcf>) {
        next unless /^chr/;
        my @fields = split;
        my ($ref, $alt) = @fields[1,2];
        my ($rev_ref, $rev_alt) = rev_and_trim( \$ref, \$alt );
        my ($norm_ref, $norm_alt) = rev_and_trim( \$rev_ref, \$rev_alt );
        
        # Reassign REF and ALT alleles.
        $fields[1] = $norm_ref;
        $fields[2] = $norm_alt;

        $vcf_data{join(':', @fields[0..2])} = [@fields];
    } 
    return %vcf_data;
}

sub rev_and_trim {
    # reverse and trim common sequence from both ends of ref and alt seqs, so that we can remove the anchor bases.
    no warnings; # Turn off so we don't have to deal with empty var warnings.
    my ($ref, $alt) = @_;
    
    my @rev_ref = split(//, reverse($$ref));
    my @rev_alt = split(//, reverse($$alt));

    while ($rev_ref[0] eq $rev_alt[0]) {
        shift @rev_ref;
        shift @rev_alt;
    }

    # If there's nothing left, replace with a hyphen.
    $rev_ref[0] //= '-';
    $rev_alt[0] //= '-';

    return (join('', @rev_ref), join('', @rev_alt));
}

sub validate_lookup {
    my $file_ref = shift;
    my $lfile_path = dirname(abs_path($0));

    # Instead of hard coding tables, use a file with tables to load.  Easier when we add
    # new ones downstream.
    open (my $index_fh, "<", "$lfile_path/lookup_tables/lookup_file_index.txt");
    my %lfiles = map{ chomp; split(/,/) } <$index_fh>;
    close $index_fh;
    
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

    print "$info Lookup file '" . basename($lookup_file) . "' selected\n";
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
    # Do a little post processing on variant call data to add in COSMIC IDs and whatnot
    my ($cpsc_lookup, $variant_data) = @_;
    my %results;
    my %missing;

    for my $var (keys %$variant_data) {
        next unless (exists $$cpsc_lookup{$var});
        @{$results{$var}} = @{$$variant_data{$var}};

        # Do a little post processing to get the gene name and fix variant ID
        push( @{$results{$var}}, $$cpsc_lookup{$var}->[0] );

        # Have to adjust the fields since we removed all NOCALL data from final output.
        if ($results{$var}->[7] =~ /^[-.]/) {
            $results{$var}->[7] = $$cpsc_lookup{$var}->[3];
        }
    }
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
        #printf $format, @{$$found_plasmids{$var}}[0,10,1,2,5,6,7,8,9];
        printf $format, @{$$found_plasmids{$var}}[0,8,1..7];
    }

    # Print out a list of thos missing only if we have some missing.
    if ( (scalar keys %$lookup) > (scalar keys %$found_plasmids) ) {
        print "\n:::  Plasmids Missing in Sample  :::\n";
        $col_widths = field_width( $lookup, [4,5] );
        $format = "%-8s%-9s%-14s%-15s%-$${col_widths[0]}s%-$${col_widths[1]}s\n";
        printf $format, qw(Gene Chr Position VARID CDS Sequence);
        for my $plas ( sort{ versioncmp( $a, $b ) } keys %$lookup) {
            printf $format, @{$$lookup{$plas}}[0..4,7] if (! exists $$found_plasmids{$plas});
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

sub __check_prog {
    my $type = shift;
    my %prog_map = (
        'solid' => 'vcfExtractor.pl',
        'cfdna' => 'cfDNA_snv_indel_report.pl',
    );

    # Check assay type params.
    if (! grep { $type eq $_ } keys %prog_map) {
        print "$err You must choose 'cfdna' or 'solid' for the assay type only!\n";
        exit 1;
    }
    my $prog = $prog_map{$type};
        
    unless (qx(which $prog)) {
        print "$err You must have $prog installed in your \$PATH to run this script. ",
            "This script can be downloaded from:\n\thttps://github.com/drmrgd.\n\n";
        exit 1;
    }
    ($type eq 'solid') ? return $prog .= ' -Nn' : return $prog;
}

