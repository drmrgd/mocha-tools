#!/usr/bin/perl
# Script to generate a new cpsc lookup table for the cpscChecker utility based on a list 
# of COSMIC IDs that will be loaded via a flat file and a master cpsc variants table containing
# all the cpsc species that we've created so far. 
#
# 8/15/13 - D Sims
##################################################################################################

use warnings;
use strict;
use File::Basename;
use Data::Dump;

#( my $scriptname = $0 ) =~ s/^(.*\/)+//;
my $scriptname = basename($0);
my $version = "v2.0.0_062714";
my $description = <<"EOT";
Program to generate new lookup tables used by the cpscChecker utility based on a master file of all 
variants and a file containing a list of COSMIC IDs to pull from that table.  
EOT

my $usage = <<"EOT";
USAGE: $0 [options] <master_cpsc_file> <lookup_file>
    -o, --output    Send output to custom output file.  Default is STDOUT.
    -v, --version   Version information
    -h, --help      Print this help information
EOT

sub help {
	printf "%s - %s\n\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub version_info {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

my $help = 0;
my $verInfo = 0;
my $outfile;

while ( scalar( @ARGV ) > 0 ) {
	last if ( $ARGV[0] !~ /^-/ );
	my $opt = shift;
	if    ( $opt eq '-v' || $opt eq '--version' )  { $verInfo = 1; }
	elsif ( $opt eq '-h' || $opt eq '--help' )     { $help = 1; }
	elsif ( $opt eq '-o' || $opt eq '--output' )   { $outfile = shift; }
	else {
		print "Invalid option: '$opt'\n";
		print "$usage\n";
		exit 1;
	}
}

&help if ( $help );
&version_info if $verInfo;


if ( scalar( @ARGV ) != 2 ) {
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
#chomp( my @cpsc_vars = <$cpsc_fh> );
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

my @err;
for my $cosid ( @cosids ) {
    if ( exists $plas_data{$cosid} ) {
        print $out_fh $plas_data{$cosid};
    } else {
       push( @err, "WARNING: '$cosid' does not exist in the master lookup table!\n" );
    }
}

# Output errors at the bottom so we can see them with long output
print "\n", join( "\n", @err ) if @err;
