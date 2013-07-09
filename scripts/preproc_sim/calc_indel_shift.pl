#!/bin/env perl

# Determine placement ambiguity for each indel in a simulated transcriptome

use warnings;
use strict;
use FileHandle;
use Getopt::Std;
use Genoman::File::TwoBit;

# The evaluation script (calc_sim_accuracy.pl) uses the results from
# this script to allow flexibility in indel placement.

# We only allow shifts _within_ the same exon because changes outside are not functionally equivalent
# Also, it will often be possible to determine in which exon an indel should be if computing over reads,
# as there are alternative isoforms and retained introns in the data

my $GENOME_FN = "/nfs/nobackup/research/bertone/common/ucsc/gbdb/hg19/hg19.2bit";
my $VERBOSE;
use constant READ_LEN => 76;
use constant MAX_DEL => 20;

main();
exit(0);


# Function to report usage information
sub usage {
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Determine placement ambiguity for each indel in a simulated transcriptome

Usage: perl $cmd [options] indel-file tx-file

Arguments:
1. indel file produced by the simulator
2. transcript file produced by the simulator

The files can be gzipped.

Options:
-v    Verbose output

The script checks both indels and short introns (<= 20 nt).

The output is tab-delimited and written to stdout. For each gap that
can be shifted, the following data is given:
1. Chromosome
2. Genomic coordinate of indel
3. Indel size (positive if insertion, negative if deletion)
4. Left shift
5. Right shift 

EndOfUsage
    exit(1);
}


# Main function
sub main {   
    # Parse args
    my %args;
    getopts('v', \%args);
    $VERBOSE = $args{'v'};
    @ARGV == 2 or usage();
    my ($indel_fn, $tx_fn) = @ARGV;

    my $genome = Genoman::File::TwoBit->new(-file => $GENOME_FN);

    process_indel_file($indel_fn, $genome);
    process_tx_file($tx_fn, $genome);
}


sub process_indel_file {
    my ($indel_fn, $genome) = @_;

    my $fh = must_open_infile($indel_fn);

    my $prev_exon_loc = "";
    my ($chr, $start, $end, $exon_seq);
    my $i = 0;

    while(my $line = <$fh>) {
	chomp $line;

	# Parse indel line
        # Positive size means insertion, negative means deletion
        # Base & coords always refer to + strand
	my ($exon_loc, $indel_pos, $indel_size, $indel_seq) = split "\t", $line;

	# Get exon sequence (if not same as for previous record)
	if($exon_loc ne $prev_exon_loc) {
	    ($chr, $start, $end) = $exon_loc =~ /^(\w+):(\d+)-(\d+)$/;
	    defined($end) or die "Failed to parse exon loc string $exon_loc";
	    $exon_seq = uc $genome->get_seq_string($chr, $start, $end);
	    $prev_exon_loc = $exon_loc;
	}

	# Compute shift
	my ($left_shift, $right_shift) = compute_indel_shift($exon_loc, $exon_seq, $indel_pos, $indel_size, $indel_seq);

	# Output shift if any
	if($left_shift or $right_shift) {
	    print join("\t", $chr, $start + $indel_pos, $indel_size, $left_shift, $right_shift), "\n";
	}

	# Progress report
	$i++;
	if($VERBOSE) {
	    print STDERR "Processed $i indels\n" if($i % 1000 == 0);
	}
    }

    print STDERR "Processed $i indels (all)\n";

    close $fh;
}


sub process_tx_file {
    my ($tx_fn, $genome) = @_;

    my $junctions = get_short_introns($tx_fn);
    
    print STDERR "Found ", scalar(@$junctions), " short introns\n";

    my $i = 0;

    foreach my $junc (@$junctions) {

	# Get coords of the intron and its two associated exons
	my ($chr, $exon1_start, $exon1_end, $exon2_start, $exon2_end) = split ':', $junc;
	my $intron_pos = $exon1_end - $exon1_start + 1;  # check
	my $intron_size = $exon2_start - $exon1_end - 1;

	# Get sequence of exon1-intron-exon2 region, and of just the intron
	my $full_seq = uc $genome->get_seq_string($chr, $exon1_start, $exon2_end);
	my $intron_seq = substr($full_seq, $intron_pos, $intron_size);

	# Compute shift
	my ($left_shift, $right_shift) = compute_indel_shift($junc, $full_seq, $intron_pos, -$intron_size, $intron_seq);

	#print STDERR "$junc\n$intron_seq\n$left_shift $right_shift\n";

	# Output shift if any
	if($left_shift or $right_shift) {
	    print join("\t", $chr, $exon1_end + 1, -$intron_size, $left_shift, $right_shift), "\n";
	}

	# Progress report
	$i++;
	if($VERBOSE) {
	    print STDERR "Processed $i introns\n" if($i % 1000 == 0);
	}
    }

    print STDERR "Processed $i introns (all)\n";

}


sub get_short_introns {
    my $tx_fn = shift;

    my %junctions;

    my $fh = must_open_infile($tx_fn);

    ## Input should begin with ---------
    my $line = <$fh>;
    unless($line =~ /^---/) {
	die "ERROR: illegal input format";
    }

    ## Process one transcript per iteration
    while(defined $line) {

	# Get transcript id
	$line = <$fh>;
	my ($id) = $line =~ /^genes\[(\d+)\] =/;
	unless(defined($id)) {
	    die "ERROR: failed to parse transcript ID from line: $line";
	}
	
	# Parse exon starts and ends
	$line = <$fh>;
	($line) = $line =~ /^starts = ([\d,]+),$/;
	my @starts = map { $_ + 1 } split ",", $line;
	$line = <$fh>;
	($line) = $line =~ /^ends = ([\d,]+),$/;
	my @ends = split ",", $line;
	
	# Ignore strand
	$line = <$fh>;
	
	# Parse chr
	$line = <$fh>;
	chomp $line;
	my ($chr) = $line =~ /^chr = (\w+)$/;
	defined($chr) or die "ERROR: failed to parse chr from line: $line";
	
	# Skip exon coords (redundant with starts and ends read above)
	while($line = <$fh>) {
	    last if($line =~ /^---/);
	}

	# Store short introns
	# We store the coords of each of the bordering exons
	for my $i (0..@starts-2) {
	    my $intron_size = $starts[$i+1] - $ends[$i] - 1;
	    if($intron_size <= MAX_DEL) {
		($intron_size > 0) or die "Intron of size $intron_size found";
		my $key = join(":", $chr, $starts[$i], $ends[$i], $starts[$i+1], $ends[$i+1]);
		$junctions{$key} = 1;
	    }
	}
    }

    $fh->close();

    return [ keys %junctions ];
}


sub compute_indel_shift {
    my ($exon_id, $exon_seq, $indel_pos, $indel_size, $indel_seq) = @_;

    # Extract sequences to the left and right of indel ($c1, $c2)
    my ($c1, $c2);
    if($indel_size < 0) { # Deletion
	$indel_size = -$indel_size;
	my $genomic_indel_seq = substr($exon_seq, $indel_pos, $indel_size);
	unless($genomic_indel_seq eq $indel_seq) {
	    die "Deletion does not match genomic sequence ($indel_seq / $genomic_indel_seq) for region $exon_id pos $indel_pos";
	}
	$c1 = substr($exon_seq, 0, $indel_pos);
	$c2 = substr($exon_seq, $indel_pos + $indel_size);
    }
    else {  # Insertion
	$c1 = substr($exon_seq, 0, $indel_pos);
	$c2 = substr($exon_seq, $indel_pos);
    }

    # Check that provided indel sequence and size match
    unless(length($indel_seq) == $indel_size) {
	die "Indel size mismatch ($indel_seq / $indel_size) for region $exon_id pos $indel_pos";
    }

    # Determine left shift
    my $left_shift = 0;
    my $i = $indel_size - 1;
    while(length($c1)) {
	last if(substr($c1, -1, 1, "") ne substr($indel_seq, $i, 1));
	$left_shift++;
	$i--;
	$i = $indel_size - 1 if($i == -1);
    }

    # Determine right shift
    my $right_shift = 0;
    $i = 0; 
    while(length($c2)) {
	last if(substr($c2, 0, 1, "") ne substr($indel_seq, $i, 1));
	$right_shift++;
	$i++;
	$i = 0 if($i == $indel_size);
    }

    # Return left and right shift
    return ($left_shift, $right_shift);
}


sub must_open_infile {
    my $fn = shift;
    my $fh;
    if($fn =~ /\.gz$/i) {
	$fh = FileHandle->new("zcat $fn |") or die "ERROR: could not open file $fn with zcat\n";
    }
    else {
	$fh = FileHandle->new($fn) or die "ERROR: could not open file $fn for reading\n";
    }
    return $fh;
}
