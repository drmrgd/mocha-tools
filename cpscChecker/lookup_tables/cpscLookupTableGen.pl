#!/usr/bin/perl
# Script to generate a new cpsc lookup table for the cpscChecker utility based on a list 
# of COSMIC IDs that will be loaded via a flat file and a master cpsc variants table containing
# all the cpsc species that we've created so far. 
#
# 8/15/13 - D Sims
##################################################################################################
use warnings;
use strict;
use autodie;

use File::Basename;
use Data::Dump;
use Term::ANSIColor;
use Getopt::Long qw(:config bundling auto_abbrev no_ignore_case);

my $scriptname = basename($0);
my $version = "v3.0.0_052815";
my $description = <<"EOT";
Program to generate new lookup tables used by the cpscChecker utility based on a master file of all 
variants and a file containing a list of COSMIC IDs to pull from that table.  
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <master_cpsc_file> <lookup_file>
    -o, --output    Send output to custom output file.  Default is STDOUT.
    -v, --version   Version information
    -h, --help      Print this help information
EOT

my $help;
my $ver_info;
my $outfile;

GetOptions( "outfile|o=s"  => \$outfile,
            "version|v"    => \$ver_info,
            "help|h"       => \$help,
    ) 
        or die $usage;

sub help {
	printf "%s - %s\n\n%s\n%s", $scriptname, $version, $description, $usage;
	exit;
}

sub version_info {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

help if $help;
version_info if $ver_info;


if ( scalar( @ARGV ) < 2 ) {
	print "ERROR: Not enough arguments passed to script!\n\n";
	print "$usage\n";
	exit 1;
}

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
	open( $out_fh, ">", $outfile ) || die "Can't open the output file '$outfile' for writing: $!";
} else {
	$out_fh = \*STDOUT;
}

#########------------------------------ END ARG Parsing ---------------------------------#########

# Load in the master cpsc variants table and check to make sure it's OK
my $master_cpsc_file = shift;
open( my $cpsc_fh, "<", $master_cpsc_file ) || die "Can not open the master CPSC file for reading: $!";
my %plas_data = map { /(COSM\d+)/; $1 => $_ } <$cpsc_fh>;
close $cpsc_fh;

# Load in the cosmic variants list to extract and make sure it's formatted correctly.
my $cos_lookup_file = shift;
my @cosids;
open( my $lookup_fh, "<", $cos_lookup_file ) || die "Can't open the COSMIC Lookup file for reading: $!";
while (<$lookup_fh>) {
	chomp;
	unless( /^COSM\d+/ ) {
		print "ERROR: The COSMIC lookup file does not appear to be the correct format (COSM####).\n";
		exit 1;
	} else {
		push( @cosids, $_ );
	}
}

my (@err, @found_list);
for my $cosid ( @cosids ) {
    if ( exists $plas_data{$cosid} ) {
        push( @found_list, $plas_data{$cosid} ); # Create array that I can sort later
    } else {
       push( @err, "WARNING: '$cosid' does not exist in the master lookup table!\n" );
    }
}

# Add sort routine to make output cleaner
my @sorted_found = map  { $_->[0] }
                   sort { $a->[1] cmp $b->[1] } # Sort by gene name
                   sort { $a->[3] <=> $b->[3] } # Then sort by position
                   map  { [$_, split(/\t/, $_)] } @found_list;

# Output errors at the bottom so we can see them with long output
print {$out_fh} $_ for @sorted_found;
#print color('bold red');
#print "\n", join( "\n", @err ) if @err;
if (@err) {
    print "\n";
    print colored($_, 'bold red on_black') for @err;
}

