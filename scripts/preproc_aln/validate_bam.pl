#!/bin/env perl

# Check that BAM files comply with format requirements

# TODO:
#  - generalize this script to run on other data sets
#  - eliminate hard-coded paths
#  - write better usage message

use warnings;
use strict;
use Getopt::Std;
use File::Basename;
use FileHandle;

# Definitions
use constant MAX_COUNT => 100;
use constant HUMAN_TARGET_FN => "/nfs/research2/bertone/external/RGASP3/targets/hg19.sizes";
use constant MOUSE_TARGET_FN => "/nfs/research2/bertone/external/RGASP3/targets/mm9.sizes";
use constant FASTQ_DIR => "/nfs/research2/bertone/external/RGASP3/reads";
# ^ note: we should not read from the FTP directory as that filesystem is not suited for many parallell reads
my %FASTQ_FILES = ('mouse' => FASTQ_DIR."/AdamsLab_mouse_brain/Mm_C57B6_brain_*.txt*",
		   'LID8465' => FASTQ_DIR."/GingerasLab_K562/cytosol/LID8465/LID8465_*.txt*",
		   'LID8466' => FASTQ_DIR."/GingerasLab_K562/cytosol/LID8466/LID8466_*.txt*",
		   'LID8556' => FASTQ_DIR."/GingerasLab_K562/nucleus/LID8556/LID8556_*.txt*",
		   'LID8557' => FASTQ_DIR."/GingerasLab_K562/nucleus/LID8557/LID8557_*.txt*",
		   'LID16627' => FASTQ_DIR."/GingerasLab_K562/whole_cell/LID16627/LID16627_*.txt*",
		   'LID16628' => FASTQ_DIR."/GingerasLab_K562/whole_cell/LID16628/LID16628_*.txt*",
		   'sim1' => FASTQ_DIR."/GrantLab_human_simulated/Human_simulated_reads_1.fq*",
		   'sim2' => FASTQ_DIR."/GrantLab_human_simulated/Human_simulated_reads_2.fq*");
my %VALID_CIGAR_OPS = map { $_ => 0 } ("M","I","D","N","S","H","P","=","X");

# Verbose mode
my $VERBOSE = 0;

# Global indexing valid target sequence IDs
my (%VALID_TARGETS);

# Globals to store counts
my ($N_UNKNOWN_READ_IDS, $N_MATE_NR_ERROR, $N_MUL_PRIMARY, $N_UNMAPPED, $N_SHORT_SEQ, $N_LONG_SEQ) = (0, 0, 0, 0, 0, 0);
my ($N_SPLICED, $N_XS_SPLICED, $N_XS_UNSPLICED) = (0, 0, 0);
my (%READS, %TARGETS, %CIGAR_OPS);

main();
exit(0);

sub usage {
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Check that BAM files comply with format requirements:

 - read IDs (match to input files)
 - flags (valid?)
 - valid counts of primary and secondary alignments per read
 - cigar strings (what operations used)
 - chromosome names
 - has XS flags?
 ... and more?

This script is currently hardcoded to work with RGASP data on the EBI
filesystem.

Usage: perl $cmd [options] in.bam outfile-prefix

Options:
-v   Verbose

EndOfUsage
    exit(1);
}


sub main {

    my %opts;
    getopts('v',\%opts);
    $VERBOSE ||= $opts{'v'};

    @ARGV == 2 or usage();

    my ($bam_fn, $report_fn_prefix) = @ARGV;
    my $dataset = get_dataset_name($bam_fn);

    # Index the target sequence IDs in the global hash %VALID_TARGETS
    progress_msg("Reading target IDs...");
    get_target_ids($dataset);

    # Index the read IDs in the global hash %READS
    progress_msg("Reading read IDs...");
    get_read_ids($dataset);

    progress_msg("Processing $bam_fn...");
    validate_bam($bam_fn);

    progress_msg("Creating report...");
    create_report($report_fn_prefix, $dataset);

    progress_msg("Done!");
}


sub validate_bam {
    my ($bam_fn) = @_;

    # Check each alignment
    open(IN, "samtools view $bam_fn |") or die "Failed to open file $bam_fn with samtools";
    while(my $line = <IN>) {

	# Parse line
	chomp $line;
	my ($query, $flags, $target, undef, undef, $cigar, undef, undef, undef, $seq, undef, @tags) = split "\t", $line;

	# Check if we know about the read ID
	unless(exists $READS{$query}) {
	    $N_UNKNOWN_READ_IDS++;
	    if($N_UNKNOWN_READ_IDS < 11) {
		print STDERR "Unknown read ID: $query\n";
		print STDERR "(will not report further unknown read IDs)\n" if($N_UNKNOWN_READ_IDS == 10);
	    }
	    next;
	}

	# Check that a mate nr is present
	unless(($flags & 0xc0) == 0x40 or ($flags & 0xc0) == 0x80) {
	    $N_MATE_NR_ERROR++;
	    next;
	}

	# Check if an alignment
	if($flags & 4) {
	    $N_UNMAPPED++;
	    next;
	}

	# Count CIGAR operations and compute sequence length
	my $seq_len = 0;
	my $is_spliced;
	my @cig_fields = split(/([A-Z])/, $cigar);
	while(@cig_fields) {
	    my $op_len = shift @cig_fields;
	    my $op = shift @cig_fields;
	    $CIGAR_OPS{$op}++;
	    if($op eq 'M' or $op eq 'I' or $op eq 'S' or $op eq '=' or $op eq 'X') {
		$seq_len += $op_len;
	    }
	    elsif($op eq 'N') {
		$is_spliced = 1;
	    }
	}

	# Count as spliced alignment?
	$N_SPLICED++ if($is_spliced);

	# Check computed sequence length
	if($seq ne '*' and $seq_len != length($seq)) {
	    die "ERROR: sequence length $seq_len computed from CIGAR ($cigar) does not match sequence ($seq) for read $query"; 
	}
	if($seq_len < 76) {
	    $N_SHORT_SEQ++;
	    if($N_SHORT_SEQ < 11) {
		print STDERR "Short sequence: $line";
		print STDERR "(will not report further short sequences)\n" if($N_SHORT_SEQ == 10);
	    }
	}
	elsif($seq_len > 76) {
	    $N_LONG_SEQ++;
	    if($N_LONG_SEQ < 11) {
		print STDERR "Long sequence: $line";
		print STDERR "(will not report further long sequences)\n" if($N_LONG_SEQ == 10);
	    }
	}

	# Check for XS tag
	foreach my $tag (@tags) {
	    if($tag eq 'XS:A:+' or $tag eq 'XS:A:-') {
		if($is_spliced) {
		    $N_XS_SPLICED++;
		}
		else {
		    $N_XS_UNSPLICED++;
		}
		last;
	    }
	}

	# Count the target ID
	$TARGETS{$target}++;

	# Count the read ID
	# We encode counts for both mates in a single value to save some memory
	# Bit 1: set if primary alignment seen for mate 1
	# Bits 2-16: counts for mate 1
	# Bit 17: set if primary alignment seen for mate 2
	# Bits 18-32: counts for mate 2
	if($flags & 0x40) {
	    $READS{$query} += 1 << 1 unless(($READS{$query} & 0xfffe) == 0xfffe);
	    unless($flags & 0x100) {
		if($READS{$query} & 1) {
		    $N_MUL_PRIMARY++;
		    if($N_MUL_PRIMARY < 11) {
			print STDERR "Multiple primary alignments for read $query mate 1\n";
			print STDERR "(will not report further additional primary alignment)\n" if($N_MUL_PRIMARY == 10);
		    }
		}
		else {
		    $READS{$query} |= 1;
		}
	    }
	}
	else {
	    $READS{$query} += 1 << 17 unless(($READS{$query} & 0xfffe0000) == 0xfffe0000);
	    unless($flags & 0x100) {
		if($READS{$query} & 0x10000) {
		    $N_MUL_PRIMARY++;
		    if($N_MUL_PRIMARY < 11) {
			print STDERR "Multiple primary alignments for read $query mate 2\n";
			print STDERR "(will not report further additional primary alignment)\n" if($N_MUL_PRIMARY == 10);
		    }
		}
		else {
		    $READS{$query} |= 0x10000;
		}
	    }
	}
    }

    close(IN) or die "ERROR: failed to read $bam_fn with samtools (exit code $?)\n";
}


sub create_report {
    my ($out_fn_prefix, $dataset) = @_;
    my $fh;

    # Write counts for target seq IDs
    write_counts("$out_fn_prefix.targets.txt", \%TARGETS, \%VALID_TARGETS);

    # Write counts for CIGAR operations
    write_counts("$out_fn_prefix.cigar_ops.txt", \%CIGAR_OPS, \%VALID_CIGAR_OPS);

    # Count counts of alignment per read
    my (@aln_counts, @aln_counts_pri, %lane_counts);
    my ($n_read_only_secondary, $n_frag_only_secondary, $max_aln) = (0, 0, 0);
    while (my ($read_id, $combined_count) = each %READS) {

	# Get lane ID
	my $lane_id;
	if($dataset eq "sim1" or $dataset eq "sim2") {
	    $lane_id = $dataset;
	}
	else {
	    ($lane_id) = $read_id =~ /^(\w+:\d+):/;
	    defined($lane_id) or die "ERROR: failed to parse lane ID from read ID: $read_id\n";    
	}

	# Determine number of alignments for the read
	my $count1 = ($combined_count & 0xfffe) >> 1;
	my $count2 = ($combined_count & 0xfffe0000) >> 17;

	# Increase fragment count for the lane
	$lane_counts{$lane_id}[0]++;

	# Check if each mate was aligned
	if($count1) {
	    $max_aln = $count1 if($max_aln < $count1); # Record the maximal number alignments seen for any read
	    $count1 = MAX_COUNT if($count1 > MAX_COUNT); # Cap alignment count to avoid enormous data structure
	    $lane_counts{$lane_id}[1]++;  # Increase mate 1 alignment count for lane
	}
	if($count2) {
	    $max_aln = $count2 if($max_aln < $count2); # Record the maximal number alignments seen for any read
	    $count2 = MAX_COUNT if($count2 > MAX_COUNT); # Cap alignment count to avoid enormous data structure
	    $lane_counts{$lane_id}[2]++;  # Increase mate 2 alignment count for lane
	}

	# Check for lack of primary alignments
	if(($count1 and !($combined_count & 1)) or ($count2 and (!$combined_count & 0x10000))) {
	    if($combined_count & 0x10001) { # Primary alignment for one of the mates?
		$n_read_only_secondary++;
		if($n_read_only_secondary < 11) {
		    print STDERR "Only secondary alignment for read $read_id\n";
		    print STDERR "(will not report further lack of primary alignments for reads)\n" if($n_read_only_secondary == 10);
		}
	    }
	    else {  # No primary alignment for either mate
		$n_frag_only_secondary++;
		if($n_frag_only_secondary < 11) {
		    print STDERR "Only secondary alignment for fragment $read_id\n";
		    print STDERR "(will not report further lack of primary alignments for fragments)\n" if($n_frag_only_secondary == 10);
		}
	    }
	}

	# Record number of alignments for the read pair
	$aln_counts[$count1][$count2]++;    # Alignment count irrespective of whether any primary alignment was reported
	$count1 = 0 unless($combined_count & 1);
	$count2 = 0 unless($combined_count & 0x10000);
	$aln_counts_pri[$count1][$count2]++;  # Alignment count requiring primary alignment
    }

    # Write the alignment counts
    write_pair_counts("$out_fn_prefix.aln_counts.txt", \@aln_counts_pri, \@aln_counts, MAX_COUNT, MAX_COUNT);

    # Write counts per lane of total sequenced and aligned reads
    $fh = must_open_outfile("$out_fn_prefix.lane_counts.txt");
    print $fh "lane\ttotal\tmate1\tmate2\n";
    foreach my $lane_id (sort keys %lane_counts) {
	my $c = $lane_counts{$lane_id};
	print $fh join("\t", $lane_id, $c->[0], ($c->[1]||0), ($c->[2]||0)), "\n";
    }
    close $fh;

    # Write misc stats
    my $n_alignments = sum_counts(\%TARGETS);
    $fh = must_open_outfile("$out_fn_prefix.misc.txt");
    print $fh "Unknown read ID\t", $N_UNKNOWN_READ_IDS, "\n";
    print $fh "Mate number error\t", $N_MATE_NR_ERROR, "\n";
    print $fh "Non-alignment\t", $N_UNMAPPED, "\n";
    print $fh "Extra primary alignment\t", $N_MUL_PRIMARY, "\n";
    print $fh "Only secondary alignment for mate\t", $n_read_only_secondary, "\n";
    print $fh "Only secondary alignment for fragment\t", $n_frag_only_secondary, "\n";
    print $fh "Max alignments for any read\t", $max_aln, "\n";
    print $fh "Short sequence\t", $N_SHORT_SEQ, "\n";
    print $fh "Long sequence\t", $N_LONG_SEQ, "\n";
    print $fh "Spliced\t", $N_SPLICED, "\n";
    print $fh "Spliced XS\t", $N_XS_SPLICED, "\n";
    print $fh "Unspliced\t", $n_alignments - $N_SPLICED, "\n";
    print $fh "Unspliced XS\t", $N_XS_UNSPLICED, "\n";
    close $fh;

}


sub sum_counts {
    my $counts = shift;
    my $sum = 0;
    foreach my $count (values %$counts) {
	$sum += $count;
    }
    return $sum;
}


sub write_counts {
    my ($fn, $counts, $valid_ids) = @_;

    my $fh = must_open_outfile($fn);
    print $fh "id\tvalid\tcount\n";
    foreach my $id (sort keys %$valid_ids) {
	print $fh join("\t", $id, "TRUE", $counts->{$id}||0), "\n";
    }
    foreach my $id (sort keys %$counts) {
	print $fh join("\t", $id, "FALSE", $counts->{$id}), "\n" unless(exists $valid_ids->{$id});
    }
    close $fh;
}


sub write_pair_counts {
    my ($fn, $counts_pri, $counts_all, $max_i, $max_j) = @_;

    # unless(defined $max_i) {
    # 	$max_i = $#$counts;
    # }
    # unless(defined $max_j) {
    # 	$max_j = 0;
    # 	for my $i (0..$max_i) {
    # 	    my $n = $#{$counts->[$i]};
    # 	    $max_j = $n if($max_j < $n);
    # 	}
    # }

    my $fh = must_open_outfile($fn);
    print $fh join("\t", "mate1", "mate2", "count.pri", "count.all"), "\n";
    for my $i (0..$max_i) {
	for my $j (0..$max_j) {
	    print $fh join("\t", $i, $j, ($counts_pri->[$i][$j]||0), ($counts_all->[$i][$j]||0)), "\n";
	}
    }
    close $fh;
}


sub get_target_ids {
    my $dataset = shift;

    my $fn = $dataset eq 'mouse' ? MOUSE_TARGET_FN : HUMAN_TARGET_FN;

    open IN, $fn or die "Failed to open file $fn for input\n";
    while(my $line = <IN>) {
	chomp $line;
	my ($id) = split "\t", $line;
	defined($id) or die "Failed to parse file $fn\n";
	$VALID_TARGETS{$id} = 0;
    }
    close IN;
}


sub get_read_ids {
    my $dataset = shift;

    # Get list of files to read
    my $fn_expr = $FASTQ_FILES{$dataset} or die "Unknown dataset: $dataset\n";
    my @files = glob($fn_expr);
    scalar(@files) or die "No files found: $fn_expr";

    # Read sequence IDs from files
    foreach my $fn (@files) {
	progress_msg("file: $fn");
	if($fn =~ /.gz/i) {
	    open(FASTQ, "zcat $fn |") or die "ERROR: could not open file $fn with zcat\n";
	}
	else {
	    open(FASTQ, $fn) or die "ERROR: could not open file $fn for reading\n";
	}
	if($dataset eq 'sim1' or $dataset eq 'sim2') {
	    while(my $id = <FASTQ>) {  # Read ID line
		##last if(scalar(keys %READS) == 100);
		chomp $id;
		($id) = $id =~ /^@(seq\.\d+)([ab])$/;
		defined($id) or die "Illegal fastq header";
		$READS{$id} = 0; # Initialize count to 0
		<FASTQ>; <FASTQ>; <FASTQ>; # Skip sequence and quality lines
	    }
	}
	else {
	    while(my $id = <FASTQ>) {  # Read ID line
		##last if(scalar(keys %READS) == 100);
		chomp $id;
		($id) = $id =~ /^@([^\/]+)(\/[12])$/;
		defined($id) or die "Illegal fastq header";
		$READS{$id} = 0; # Initialize count to 0
		<FASTQ>; <FASTQ>; <FASTQ>; # Skip sequence and quality lines
	    }
	}
	close(FASTQ) or die "ERROR: failed to read file $fn\n";
    }
}


sub get_dataset_name {
    my $fn = shift;
    return basename($fn, "_1.bam", "_2.bam" ,"_3.bam", "_4.bam", ".bam");
}


sub must_open_outfile {
    my $fn = shift;
    return \*STDOUT if($fn eq '' or $fn eq '-');
    my $fh = FileHandle->new(">$fn") or die "ERROR: could not open file $fn for writing.\n";
    return $fh;
}


sub progress_msg
{
    $VERBOSE or return;
    my $str = shift;
    print STDERR time_string()." $str\n";
}

sub time_string
{
    my ($s, $m, $h) = localtime(time);
    return sprintf("[%02d:%02d:%02d]", $h, $m, $s);
}
