#!/usr/bin/perl
# Quickie to convert the CPSG text files to BED files for easier comparison.
use strict;
use warnings;
use autodie;

use Data::Dump;
use Sort::Versions;

my $text_file = shift or die "ERROR: You must input a text file!\n";
my @parsed_lines;

open(my $fh, "<", $text_file);
while (<$fh>) {
    my ($gene, $chr, $pos, $id, $cds, $ref, $alt, $seq) = split;
    my ($parsed_ref, $parsed_alt, $stop);
    # map an anchor style ref / alt from sequence field.
    if ($seq =~ /\[/) {
        # We have an SNV / Indel
        $parsed_ref = $ref;
        $parsed_alt = $alt;
    } 
    elsif ($seq =~ /\|/) {
        #print "We have an insertion\n";
        $parsed_ref = get_anchor($seq);
        $parsed_alt = "$parsed_ref$alt";
    }
    elsif ($seq =~ /</) {
        #print "We have a deletion.\n";
        my $anchor = get_anchor($seq);
        $parsed_ref = "$anchor$ref";
        $parsed_alt = $anchor;
    }
    else {
        #should not get here
        next;
    }
    #print "($pos) $parsed_ref -> $parsed_alt: ";
    $stop = get_start_stop($parsed_ref, $parsed_alt, $pos);
    my $info = "GENE=$gene;CDS=$cds;SEQ=$seq";
    push(@parsed_lines, [$chr, $pos, $stop, $parsed_ref, $parsed_alt, $id, $info]);
    #print join("\t", $chr, $pos, $stop, $parsed_ref, $parsed_alt, $id, $info), "\n";
}

#dd \@parsed_lines;
for my $entry (sort {versioncmp($a->[0], $b->[0]) || versioncmp($a->[1], $b->[1]) } @parsed_lines) {
    print join("\t", @$entry), "\n";
}

sub get_anchor {
    my $string = shift;
    my ($anchor) = $string =~ /.*?([ACTG])[\|\<][ACTG]+[\|\>].*/i;
    return $anchor;
}

sub get_start_stop {
    my ($ref, $alt, $pos) = @_;
    if (length($ref) <= length($alt)) {
        return $pos + 1;
    }
    elsif (length($ref) > length($alt)) {
        return $pos + (length($ref) - 1);
    }
}
