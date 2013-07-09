#!/bin/env perl

# Compute rates for incorporation of false junctions (from alignments) into transcript models

# Pragmas and includes
use warnings;
use strict;
use Getopt::Std;

my $MAX_ALN_COUNT = 10; # Max alignments to count for splice stats
 # ^ keep this as a variable in case we make it a commandline option at some point
my $REPORT_FALSE;

# Do the job and exit
main();
exit(0);


# Function to report usage information
sub usage {
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Compute rates for incorporation of false junctions (from alignments)
into transcript models

Usage: perl $cmd [options] truth.gtf[.gz] pred.gtf[.gz]

The script expects SAM on stdin.

Options:

-r   Report false junctions on STDERR

EndOfUsage
    exit(1);
}


# Main function
sub main {
    
    my %args;
    getopts("r", \%args);
    $REPORT_FALSE = $args{r};
    @ARGV == 2 or usage();
    my ($truth_fn, $pred_fn) = @ARGV;

    my $truth_index = read_splices($truth_fn);

    my $pred_index = read_splices($pred_fn);

    my $introns = read_sam();

    print_report($introns, $truth_index, $pred_index);
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



# Make report on splices
sub print_report {
    my ($introns, $truth_index, $pred_index) = @_;

    my %class_counts;
    while(my ($intron_id, $nr_aln) = each %$introns) {
	my $class =
	    ($truth_index->{$intron_id} ? "true" : "false").".".
	    ($pred_index->{$intron_id} ? "pred" : "skipped");
	$nr_aln = $MAX_ALN_COUNT if($nr_aln > $MAX_ALN_COUNT);
	$class_counts{$class}[$nr_aln-1]++;
	if($REPORT_FALSE and ($class eq "false.pred" or $class eq "false.skipped")) {
	    print STDERR $class, "\t", $intron_id, "\n";
	}
    }

    print join("\t", "class", (map { "n$_" } (1..$MAX_ALN_COUNT))), "\n"; 

    for my $class ("true.pred", "true.skipped", "false.pred", "false.skipped") {
	print $class;
	for my $i (0..($MAX_ALN_COUNT-1)) {
	    print "\t", ($class_counts{$class}[$i] || 0);
	}
	print "\n";
    }

}


# This function reads a GTF file to create annotation index:
# %splices{chr:left-right}

sub read_splices {
    my $fn = shift;

    # Read file
    my (%tx_exons, %tx_chr);
    if($fn =~ /\.gz/i) {
	open(IN, "zcat $fn |") or die "ERROR: could not open file $fn with zcat\n";
    }
    else {
	open(IN, $fn) or die "ERROR: could not open file $fn for reading\n";
    }
    while(my $line = <IN>) {
	chomp $line;
	my ($chr, undef, $ft_type, $start, $end, undef, undef, undef, $attrib) = split "\t", $line;
	# Check that the feature is relevant
	next unless($ft_type eq 'exon'); # Only consider exons
	# Get tx id
	my ($tx_id) = $attrib =~ /transcript_id "([\w\.]+)"/;
	defined($tx_id) or die "Failed to parse attribute string: $attrib";
	# Store exon coords in index
	if(exists $tx_chr{$tx_id}) {
	    $tx_chr{$tx_id} eq $chr or die "Tx $tx_id annotated on multiple chromosomes"; 
	}
	else {
	    $tx_chr{$tx_id} = $chr;
	}
	push @{$tx_exons{$tx_id}}, ($start, $end);
    }
    close IN;

    # Get exon junction coords
    my %splices;
    foreach my $tx_id (keys %tx_exons) {
	my $chr = $tx_chr{$tx_id};
	my @sorted_coords = sort { $a <=> $b } @{$tx_exons{$tx_id}};
	shift @sorted_coords; # Ignore start of first exon
	pop @sorted_coords;   # Ignore end of last exon
	while(@sorted_coords) {
	    my $left = shift @sorted_coords;
	    my $right = shift @sorted_coords;
	    $splices{"$chr:$left-$right"} = 1;
	}
    }

    # Return index
    return \%splices;
}


