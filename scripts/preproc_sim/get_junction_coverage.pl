#!/bin/env perl
# Compute the read coverage for each splice junction in a set of alignments
use warnings;
use strict;

@ARGV and usage();

my $introns = read_sam();

foreach my $intron_id (keys %$introns) {
    print join ("\t", $intron_id, $introns->{$intron_id}), "\n";
}

exit(0);


# Function to report usage information
sub usage {
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Compute the read coverage for each splice junction in a set of alignments

Usage: perl $cmd

The script expects SAM on stdin.

EndOfUsage
    exit(1);
}


# Read SAM from STDIN and collect stats
sub read_sam {

    my %introns;

    # Process one alignment per iteration
    while(my $line = <STDIN>) {

	# Parse SAM
	my (undef, $flags, $chr, $pos, undef, $cigar) = split "\t", $line;
	next if($flags & 0x4); # ignore non-alignments

	# Split up CIGAR string
	my @cig_fields = split(/([A-Z])/, $cigar);
	my $i = 0;

	# Process one CIGAR operation per iteration
	while(@cig_fields) {
	    my $l = shift @cig_fields;
	    my $op = shift @cig_fields;
	    if($op eq 'M' or $op eq 'D') {
		$pos += $l;
	    }
	    elsif($op eq 'N') {
		if($i > 0 and @cig_fields > 0) { # Check that this is an internal operation
		    my $intron_key = $chr . ":" . ($pos - 1) . "-" . ($pos + $l);
		    $introns{$intron_key}++;
		    $pos += $l;
		}
		else {
		    warn "Warning: terminal N operation in CIGAR string: $cigar"; 
		}
	    }
	    elsif($op ne 'I' and $op ne 'S') {
		die "I don't know what to do with CIGAR operation $op";
	    }
	    $i++;
	}

    }

    # Return introns
    return \%introns;
}

