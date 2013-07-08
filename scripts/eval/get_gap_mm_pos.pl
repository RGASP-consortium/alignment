#!/bin/env perl

# Determine positional distribution of indels, introns and mismatches

use warnings;
use strict;

use constant DEBUG => 0;

# The script should be run without arguments
# Print usage message if any arguments provided
if(@ARGV) {

    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Determine positional distribution of indels, introns and mismatches.

This script takes no arguments and expects SAM on stdin. The
distribution will be printed to stdout.

EndOfUsage
    exit(1);

}

# Initialize arrays to store counts
my @intron_counts = (0) x 76;
my @del_counts = (0) x 76;
my @ins_counts = (0) x 76;
my @mismatch_counts = (0) x 76;
my $terminal_introns = 0;

# Process SAM input, one line per iteration
while(my $line = <STDIN>) {

    chomp $line;
    my @sam = split /\t/, $line;

    # Ignore mitochondrial chromosome as due to difference in human chrM seq between UCSC and Ensembl
    next if($sam[2] eq 'chrM');

    # Get MD tag
    my $md_tag;
    foreach my $tag (@sam[11..$#sam]) {
	if($tag =~ /^MD:Z:(.+)/) {
	    $md_tag = $1;
	    last;
	}
    }
    defined($md_tag) or die "No MD tag for alignment: $line";

    # Check for commmon case for which there is nothing to count
    next if($sam[5] eq '76M' and $md_tag eq '76');

    # Get strand of alignment
    my $ori = $sam[1] & 0x10 ? -1 : 1;

    print STDERR "\n$sam[0] $ori $sam[5] $md_tag\n" if(DEBUG);

    # Get indel and intron pos from CIGAR
    # Also create an array @read_pos_map that has the read coordinate for each M operation
    my @cig_fields = split(/([A-Z])/, $sam[5]);
    my @read_pos_map;
    my $read_pos = $ori == 1 ? 0 : 75;
    while(@cig_fields) {
	my $l = shift @cig_fields;
	my $op = shift @cig_fields;
	if($op eq 'M') {
	    while($l) {
		push @read_pos_map, $read_pos;
		$read_pos += $ori;
		$l--;
	    }
	}
	elsif($op eq 'N') {
	    my $gap_pos = $ori == 1 ? $read_pos-1 : $read_pos; # Count at read pos before gap
	    if($gap_pos < 0 or $gap_pos > 74) {
		if($terminal_introns == 0) {
		    warn "Warning: terminal intron on line: $line";
		}
		$terminal_introns++;
	    }
	    else {
		$intron_counts[$gap_pos]++;
		print STDERR "Intron $gap_pos\n" if(DEBUG);
	    }
	}
	elsif($op eq 'D') {
	    my $gap_pos = $ori == 1 ? $read_pos-1 : $read_pos; # Count at read pos before gap
	    if($gap_pos < 0 or $gap_pos > 74) {
		die "Deletion at position $gap_pos for alignment: $line\n";
	    }
	    $del_counts[$gap_pos]++;
	    print STDERR "Deletion $gap_pos\n" if(DEBUG);
	}
	elsif($op eq 'I') {
	    while($l) {
		$ins_counts[$read_pos]++;
		print STDERR "Insertion $read_pos\n" if(DEBUG);
		$read_pos += $ori;
		$l--;
	    }
	}
	elsif($op eq 'S') {
	    $read_pos += ($ori == 1 ? $l : -$l);
	}
	else {
	    die "I don't know what to do with CIGAR operation $op";
	}
    }

    # Check that we counted 76 bases
    unless(($ori == 1 and $read_pos == 76) or ($ori == -1 and $read_pos == -1)) {
	die "Read pos counting error for line: $line\n"
    }

    # Split the MD tag into components
    # Walk through to find location of each mismatch
    my @md_fields = split /(\D+)/, $md_tag;
    my $md_tag_pos = 0;
    while(@md_fields > 1) {
	$md_tag_pos += shift @md_fields;
	my $op = shift @md_fields;
	if($op =~ /^[CATGcatg]$/) {
	    my $mm_pos = $read_pos_map[$md_tag_pos];
	    my $read_base = $ori == 1 ? substr($sam[9], $mm_pos, 1) : substr($sam[9], 75-$mm_pos, 1);
	    if($read_base ne 'N' and $read_base ne 'n') {
		print STDERR "Mismatch $mm_pos $read_base:$op\n" if(DEBUG);
		$mismatch_counts[ $mm_pos ]++;
	    }
	    $md_tag_pos++;
	}
	elsif($op eq 'N' or $op eq 'n') {
	    $md_tag_pos++;
	}
	elsif($op !~ /^\^[A-Za-z]+$/) {
	    die "Error parsing MD tag: $md_tag";
	}
    }
}

if($terminal_introns) {
    warn "Warning: $terminal_introns terminal introns ignored\n";
}

print join("\t", "intron", "mm", "del", "ins"), "\n";
for my $i (0..75) {
    print join("\t", $intron_counts[$i], $mismatch_counts[$i], $del_counts[$i], $ins_counts[$i]), "\n";
}
