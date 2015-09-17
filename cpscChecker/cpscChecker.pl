#!/usr/bin/perl
# Generate a CPSC Checker report starting with a VCF file. This replaces the original CPSC Checker
# script which required a variants.xls file to run.  
# 
# TODO:
#     - Fix issue with differential reporting betwen the VCF and TSV files.  VCF is not capturing the 
#       indel variants for some reason!
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

# Remove when in prod.
print "\n";
print colored("*" x 50, 'bold yellow on_black');
print colored("\n    DEVELOPMENT VERSION OF CPSC CHECKER\n", 'bold yellow on_black');
print colored("*" x 50, 'bold yellow on_black');
print "\n\n";

my $scriptname = basename($0);
my $version = "v3.9.10_091715-dev";
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

# Set up some formatted and colored output
my $err = colored( 'ERROR:', 'bold red on_black');
my $warn = colored( 'WARNING:', 'bold yellow on_black');
my $info = colored( 'INFO:', 'bold cyan on_black');

# Make sure enough args passed to script
if ( defined $lookup && $lookup ne '?' ) {
    unless ($tab || $vcf) {
        print "$err You must define either a VCF file with '-V' or a TSV file with '-T'!\n\n";
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
my %raw_data;
if ($vcf) {
    (-e $vcf) ? (%raw_data = read_vcf(\$vcf)) : die "$err Can't read VCF file '$vcf': $!";
} else {
    my %tab_data;
    (-e $tab) ? (%raw_data = read_tab(\$tab)) : die "$err Can't read alleles tab file '$tab': $!";
}

# Check VCF dataset against lookup hash for final results.
# XXX
my %plas_results = proc_plas_data(\%cpsc_lookup_data,\%raw_data);

# TODO
print "Got to line 134\n";
dd \%plas_results;
exit;

# Print out the results.
print_results(\%plas_results, \%cpsc_lookup_data);

sub read_tab {
    # TODO: Issue here with REF / ALT and OREF / OALT not always agreeing (see the ATM indel for example).  Need to find 
    # a better way to map this.
    my $input_file = shift;
    my %tab_data;

    print "Checking tab delimited variants file '$$input_file' for CPSC data...\n";

    open( my $tab_fh, "<", $$input_file );
    my $sample_line = <$tab_fh>;
    die "$err '$$input_file' does not appear to be a TSV file!\n" if $sample_line =~ /^#/;
    while (<$tab_fh>) {
        chomp;
        next unless /^chr/;
        my @fields = split(/\t/);
        my $var_coord = join(':', @fields[0,1]);

        #next unless $var_coord eq "chr3:178952150"; # PIK3CA insA
        next unless $var_coord eq 'chr11:108117847'; # ATM delGT
        #next unless $var_coord eq "chr7:55242466"; # EGFR del15

        next if $fields[6] == 0; # get rid of ref calls.
        my $ref_cov = $fields[18]-$fields[24];
        $tab_data{join(':', @fields[0..3])} = [$var_coord, @fields[2,3,6,18], $ref_cov, $fields[24], $fields[11]];
    }
    close $tab_fh;
    return %tab_data;
}

sub read_vcf {
    # TODO: Issue here with REF / ALT and OREF / OALT not always agreeing (see the ATM indel for example).  Need to find 
    # a better way to map this.
    my $vcf_file = shift;
    my %vcf_data;

    print "Checking VCF file '$$vcf_file' for CPSC data...\n";

    chomp(my $path = qx( which vcfExtractor.pl ));
    my $vcf_extractor_cmd = qq{ $path -Nn $$vcf_file };

    my ($parsed_vcf,$child_pid, $child_rc);
    # Better error handling to capture failed VCF Extractor output.
    unless($child_pid = open($parsed_vcf, "-|")) {
        open(STDERR, ">&STDOUT");
        exec($vcf_extractor_cmd);
        die "ERROR: could not execute program:$!";
    }
    waitpid($child_pid,0);
    $child_rc = $? >> 8;
    die "$err Can not parse VCF file '$$vcf_file'\n" if $child_rc;

    while (<$parsed_vcf>) {
        next unless /^chr/;
        my @fields = split;
        #next unless $fields[0] eq "chr3:178952150"; # PIK3CA insA
        #next unless $fields[0] eq 'chr11:108117847'; # ATM delGT
        #next unless $fields[0] eq "chr7:55242466"; # EGFR del15
        dd @fields;
        my ($ref, $alt) = @fields[1,2];
        print "original ref: $ref\n";
        print "original alt: $alt\n";
        my ($rev_ref, $rev_alt) = rev_and_trim( \$ref, \$alt );
        #print "reversed ref: $rev_ref\n";
        #print "reversed alt: $rev_alt\n";
        my ($norm_ref, $norm_alt) = rev_and_trim( \$rev_ref, \$rev_alt );
        
        #print "normalized ref: $norm_ref\n";
        #print "normalized alt: $norm_alt\n";

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

    #while (@rev_ref > 1 && @rev_alt > 1 && $rev_ref[0] eq $rev_alt[0]) {
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

    # XXX
    #dd $variant_data;
    #exit;

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
