#!/usr/bin/perl
# Script to pull together the cpsc lookup files to generate a larger one
# unfinished!
# 8/15/13 - D Sims

use warnings;
use strict;
use Data::Dumper;

my $cpsc_v4 = shift;
my $cpsc_v5 = shift;

open( V4, "<", $cpsc_v4 );
chomp( my @v4_vars = <V4> );
close(V4);

open( V5, "<", $cpsc_v5 );
chomp( my @v5_vars = <V5> );
close(V5);

my @results;
foreach my $line ( @v5_vars ) {
	push( @results, $line ) if ( ! grep { $line ~~ $_ } @v4_vars );
}

print Dumper( \@results );

