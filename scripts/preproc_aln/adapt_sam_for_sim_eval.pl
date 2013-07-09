#!/bin/env perl

# Adapt alignment files for evaluation by calc_sim_accuracy.pl or the scripts in the BEERS toolkit.

use warnings;
use strict;

# Do the job and exit
main();
exit(0);


sub usage {
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Adapt alignment files for evaluation by calc_sim_accuracy.pl or the
scripts in the BEERS toolkit.

Usage: perl $cmd

The script expects SAM on stdin and writes SAM to stdout.  Aligmments
should be sorted by fragment ID (column 1). The input should not
contain unmapped reads.

To use BAM files as input and output, the script can be run as:

samtools view -hF 0x4 \$IN_FILE | $cmd | samtools view -bSo \$OUT_FILE -

EndOfUsage
    exit(1);
}


sub main {
    (@ARGV == 0) or usage();

    my $line;

    # Pass through header lines
    while($line = <STDIN>) {
	last if(substr($line, 0, 1) ne '@');
	print $line;
    }

    my $i = 0;
    my ($pri1, $sec1, $pri2, $sec2);

    while(defined($line)) {

	# Split line.
	chomp $line;
	my @sam = split "\t", $line;

	# Get read number
	my ($j) = $sam[0] =~ /^seq.(\d+)$/;
	defined($j) or die "ERROR: Failed to parse read id $sam[0]";

	# A new ID?
	if($i != $j) {
	    process_aln($i, $pri1, $sec1, $pri2, $sec2) if($i);
	    print_unmapped_pairs($i+1, $j-1);
	    $i = $j;
	    undef $pri1;
	    undef $sec1;
	    undef $pri2;
	    undef $sec2;
	}   

	# Store this alignment
	if($sam[1] & 0x40) {
	    if($sam[1] & 0x100) {
		$sec1 = \@sam unless(defined $sec1);
	    }
	    else {
		$pri1 = \@sam;
	    }
	}
	elsif($sam[1] & 0x80) {
	    if($sam[1] & 0x100) {
		$sec2 = \@sam unless(defined $sec2);
	    }
	    else {
		$pri2 = \@sam;
	    }
	}
	else {
	    die "ERROR: No mate number in flags: $sam[1]";
	}
   
	# Read next line
	$line = <STDIN>;
    }

    process_aln($i, $pri1, $sec1, $pri2, $sec2) if($i);
    print_unmapped_pairs($i+1, 40000000);
}


sub process_aln {
    my ($seq_nr, $pri1, $sec1, $pri2, $sec2) = @_;

    if($sec1 or $sec2) {
	print_aln($pri1, $seq_nr, 1, "HI:i:1", "IH:i:2");
	print_aln($pri2, $seq_nr, 2, "HI:i:1", "IH:i:2");
	print_aln($sec1, $seq_nr, 1, "HI:i:2", "IH:i:2");
	print_aln($sec2, $seq_nr, 2, "HI:i:2", "IH:i:2");
    }
    else {
	print_aln($pri1, $seq_nr, 1, "HI:i:1", "IH:i:1");
	print_aln($pri2, $seq_nr, 2, "HI:i:1", "IH:i:1");
    }
}


sub print_unmapped_pairs {
    my ($i, $j) = @_;
    for my $seq_nr ($i..$j) {
	print_aln(undef, $seq_nr, 1, "HI:i:1", "IH:i:1");
	print_aln(undef, $seq_nr, 2, "HI:i:1", "IH:i:1");
    }
}


sub print_aln {
    my ($aln, $seq_nr, $mate_nr, @tags) = @_;
    if($aln) {
	$#$aln = 10;
	$aln->[0] .= chr(96 + $mate_nr);
	$aln->[9] = "*";
	$aln->[10] = "*";
    }
    else {
	$aln = ["seq.".$seq_nr.chr(96 + $mate_nr), $mate_nr * 0x40 + 0x4 + 0x1, "*", 0, 0, "*", "*", 0, 0, "*", "*"];
    }
    print join("\t", @$aln, @tags), "\n";
}
