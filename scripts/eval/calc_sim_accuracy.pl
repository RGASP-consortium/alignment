#!/bin/env perl
# Calc accuracy for alignments of simulated reads
use warnings;
use strict;
use FileHandle;
use Getopt::Long;
use List::Util qw(max);

use constant MAX_DEL => 18;
use constant MAX_JUNC_COV => 20;

# Indices for stats arrays
use constant {         # Basewise and read-level stats
    READS_TOTAL => 0,  # Total number of sequenced reads 
    BASE_RIGHT => 1,   # Basewise stats
    BASE_WRONG => 2,
    ALN_TOTAL => 3,    # Read level stats
    ALN_OVERLAP => 4,
    ALN_1_RIGHT => 5,
    ALN_PERFECT => 6
};
use constant {                  # Gap stats
    MI_READ_ACTUAL => 0,        # Multi-intron accuracy at the read level
    MI_READ_REPORTED => 1,
    MI_READ_CORRECT => 2,
    MI_PAIR_ACTUAL => 3,        # Multi-intron accuracy at the read pair level
    MI_PAIR_REPORTED => 4,
    MI_PAIR_CORRECT => 5,
    DEL_CORRECT => 6,           # Accuracy for individual indels
    DEL_MISSED => 7,
    DEL_FALSE => 8,
    INS_CORRECT => 9,
    INS_MISSED => 10,
    INS_FALSE => 11,
    INTRON_CORRECT => 12,       # Accuracy for individual introns
    INTRON_MISSED => 13,
    INTRON_FALSE => 14
};


# Variables to store global options
my ($AMBIG_READ, $NO_MITO, $VERBOSE);
my $REPORT_FALSE_INDELS = 1e9;
my $READ_LEN = 76;

## Call main function and exit
main();
exit(0);


# Function to report usage information
sub usage {
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Calc accuracy for alignments of simulated reads

Usage: perl $cmd [options] indels.txt junc.txt truth.bam result.bam outName

The indel file should contain indel shifts, as computed by calc_indel_shift.pl.
The file can be gzipped.

The file junc.txt should contain true coverage of splice junctions, as computed
by calc_junc_coverage.pl. The file can be gzipped.

The following 9 output files will be created in the current working directory:
sa_general.<outName>.txt        General stats
sa_indel.<outName>.txt          Indel accuracy
sa_intron_multi.<outName>.txt   Accuracy for multi-intron alignments
sa_intron_fn_amb.<outName>.txt  False negative splices (ambiguous mappings)
sa_intron_fn_uni.<outName>.txt  False negative splices (unique mappings)
sa_intron_fp_amb.<outName>.txt  False positive splices (ambiguous mappings)
sa_intron_fp_uni.<outName>.txt  False positive splices (unique mappings)
sa_intron_tp_amb.<outName>.txt  True positive splices (ambiguous mappings)
sa_intron_tp_uni.<outName>.txt  True positive splices (unique mappings)

General options:
-a      Determine ambiguous mapping per read instead of per fragment
-l int  Read length. Default: 76.
-m      Ignore mitochondrial chromosome

Diagnostic options:
-v   Verbose mode (not implemented)
--report-false-indels <min-size>   Report false indels on STDERR
--report-aln-info    Write info about each primary alignment to STDERR

EndOfUsage
    exit(1);
}


# Main function
sub main {   
    # Parse args
    GetOptions('a' => \$AMBIG_READ,
	       'l=i' => \$READ_LEN,
	       'm' => \$NO_MITO,
	       'v' => \$VERBOSE,
	       'report-false-indels=i' => \$REPORT_FALSE_INDELS);
    @ARGV == 5 or usage();
    my ($indel_fn, $junc_fn, $truth_fn, $result_fn, $out_name) = @ARGV;

    # Read indel and junction data
    my ($del_shifts, $ins_shifts) = read_indel_shifts($indel_fn);
    my $junc_coverage = read_junction_coverage($junc_fn);

    # Compute accuracy
    my $stats = calc_accuracy_for_file($truth_fn, $result_fn, $del_shifts, $ins_shifts, $junc_coverage);

    # Output stats
    write_stats($stats, $out_name);

    # Done.
}


sub write_stats {
    my ($stats, $out_name) = @_;
    
    # Report basewise and read-level stats
    my $fn = "sa_general.$out_name.txt";
    open OUT, ">$fn" or die "ERROR: failed to open $fn for output";
    print OUT join("\t", qw(cat bases.right bases.wrong aln overlap some.right perfect  total.reads total.bases)), "\n";
    foreach my $i ('uni', 'amb') {
	foreach my $j ('contig', 'spliced') {
	    my $x = $stats->{$i}{$j};
	    print OUT join("\t", "${i}_${j}", map { $_ || 0 }
			   @{$x}[ BASE_RIGHT, BASE_WRONG, ALN_TOTAL, ALN_OVERLAP, ALN_1_RIGHT, ALN_PERFECT, READS_TOTAL ],
			   ($x->[ READS_TOTAL ] || 0) * $READ_LEN), "\n";
	}
    }
    close OUT;

    # Report indel stats, stratified by indel size
    $fn = "sa_indel.$out_name.txt";
    my @col_index = (DEL_CORRECT, DEL_MISSED, DEL_FALSE, INS_CORRECT, INS_MISSED, INS_FALSE);
    my @col_names = ('del.correct', 'del.missed', 'del.false', 'ins.correct', 'ins.missed', 'ins.false');
    write_gap_stats_stratif($fn, $stats, \@col_index, \@col_names);

    # Report multi-intron stats, stratified by the number of introns per read
    $fn = "sa_intron_multi.$out_name.txt";
    @col_index = (MI_READ_ACTUAL, MI_READ_REPORTED, MI_READ_CORRECT,
		  MI_PAIR_ACTUAL, MI_PAIR_REPORTED, MI_PAIR_CORRECT);
    @col_names = ('mi.read.actual', 'mi.read.reported', 'mi.read.correct',
		  'mi.pair.actual', 'mi.pair.reported', 'mi.pair.correct');
    write_gap_stats_stratif($fn, $stats, \@col_index, \@col_names);

    # Report intron accuracy, stratified by read position and true coverage
    foreach my $cat ('uni', 'amb') {
	write_intron_accuracy("sa_intron_tp_$cat.$out_name.txt", $stats->{$cat}{gap}[INTRON_CORRECT]);
	write_intron_accuracy("sa_intron_fn_$cat.$out_name.txt", $stats->{$cat}{gap}[INTRON_MISSED]);
	write_intron_accuracy("sa_intron_fp_$cat.$out_name.txt", $stats->{$cat}{gap}[INTRON_FALSE]);
    }
}



sub write_gap_stats_stratif {
    my ($fn, $stats, $col_index_ref, $col_names_ref) = @_;

    open OUT, ">$fn" or die "ERROR: failed to open $fn for output";

    my @col_names = ((map { "uni.$_" } @$col_names_ref), (map { "amb.$_" } @$col_names_ref));
    my $i_max = max_array_size( @{$stats->{uni}{gap}}[ @$col_index_ref ],
				@{$stats->{amb}{gap}}[ @$col_index_ref ] );

    print OUT join("\t", "i", @col_names), "\n";

    for my $i (0..$i_max) {
	print OUT "$i";
	for my $cat ('uni', 'amb') {
	    for my $j (@$col_index_ref) {
		print OUT "\t", ($stats->{$cat}{gap}[$j][$i] || 0);
	    }
	}
	print OUT "\n";
    }

    close OUT;
}


sub write_intron_accuracy {
    my ($fn, $stats) = @_;

    open OUT, ">$fn" or die "ERROR: failed to open $fn for output";

    print OUT join("\t", "", 0..MAX_JUNC_COV), "\n";

    for my $i (0..($READ_LEN-2)) {
	print OUT $i+1;
	for my $j (0..MAX_JUNC_COV) {
	    print OUT "\t", ($stats->[$i][$j] || 0);
	}
	print OUT "\n";
    }

    close OUT;
}


sub calc_accuracy_for_file {
    my ($truth_fn, $result_fn, $del_shifts, $ins_shifts, $junc_coverage) = @_;
    
    # Stats to collect
    my $stats_unique = { contig => [], spliced => [], gap => [] };
    my $stats_ambig = { contig => [], spliced => [], gap => [] };

    # Open alignment files
    my $truth_fh = must_open_aln_file($truth_fn);
    my $result_fh = must_open_aln_file($result_fn);
    
    # Process each fragment
    while(my $truth_line1 = <$truth_fh>) {
	my $truth_line2 = <$truth_fh>;
	my $result_line1 = <$result_fh>;
	my $result_line2 = <$result_fh>;
	my ($stats1, $stats2, $stats_pair);

	# Check if ambig-mapped
	my ($ih) = $result_line1 =~ /\tIH:i:([12])$/;
	defined($ih) or die "Failed to parse sam line: $result_line1\n";
	if($ih == 1) {
	    $stats1 = $stats2 = $stats_pair = $stats_unique;
	}
	elsif($ih == 2) {
	    my $result_line3 = <$result_fh>;
	    my $result_line4 = <$result_fh>;
	    if($AMBIG_READ) { # Determine ambiguous mapping per read
		$stats1 = is_aligned($result_line3) ? $stats_ambig : $stats_unique;
		$stats2 = is_aligned($result_line4) ? $stats_ambig : $stats_unique;
	    }
	    else { # Determine ambiguous mapping per fragment
		$stats1 = $stats2 = $stats_ambig;
	    }
	    $stats_pair = $stats_ambig;
	}

	# Compute read level accuracy. Pass $stats1 and $stats2 for recording result.
	my @introns1 = calc_accuracy_for_aln($truth_line1, $result_line1, $stats1, $del_shifts, $ins_shifts, $junc_coverage);
	my @introns2 = calc_accuracy_for_aln($truth_line2, $result_line2, $stats2, $del_shifts, $ins_shifts, $junc_coverage);

	# Compute read pair level multi-intron stats 
	if(@introns1 and @introns2) {
	    calc_multi_intron_read_pair_stats(@introns1, @introns2, $stats_pair->{gap});
	}
	elsif(@introns1 or @introns2) {
	    die "ERROR: cannot compute pair statistics when only one read is to be counted.\n";
	}
    }

    # Close input files
    close($truth_fh) or die "Failed to close truth BAM file $truth_fn: $!";
    close($result_fh) or die "Failed to close result BAM file $result_fn: $!";

    # Return stats
    return { 'uni' => $stats_unique, 'amb' => $stats_ambig };
}


# This is the core function for determining accuracy
# Args:
# 1. SAM line from simulator, describing the true alignment
# 2. SAM line from aligner output
# 3. Reference to array where result is stored
# 4. Shift information for deletions
# 5. Shift information for insertions
# 6. True coverage for splice junctions

sub calc_accuracy_for_aln {
    my ($truth_line, $result_line, $stats, $del_shifts, $ins_shifts, $junc_coverage) = @_;
    
    # Parse SAM lines
    my ($read_id1, $flags1, $chr1, $start1, undef, $cigar1) = split "\t", $truth_line;
    my ($read_id2, $flags2, $chr2, $start2, undef, $cigar2) = split "\t", $result_line;

    # Confirm that read IDs match
    if($flags1 & 0x40) {
	$read_id1 .= "a";
    }
    elsif($flags1 & 0x80) {
	$read_id1 .= "b";
    }
    else {
	die "Both mate number flags unset for truth SAM line: $truth_line";
    }
    if($read_id1 ne $read_id2) {
	die "ERROR: read ID mismatch: $read_id1 $read_id2";
    }

    # Ignore mitochondrial reads?
    return () if($NO_MITO and $chr1 eq 'chrM');

    # For unaligned reads, make sure chr, pos and cigar are unset
    if($flags2 & 4) {
	($chr2, $start2, $cigar2) = ("", 0, "");
    }

    # Extract strand flag for each alignment
    my $strand1 = $flags1 & 0x10;
    my $strand2 = $flags2 & 0x10;

    # Get slot for gap stats (indel and intron)
    my $gap_stats = $stats->{ 'gap' };

    # Declare variables for storing alignment data 
    my ($coords1, $coords2, $introns1, $introns2); # Coordinates for individual read bases and introns
    my ($base_right_count, $base_wrong_count); # Nr of correct and misplaced bases

    # Compare true and reported alignments
    
    # First check if chromosome and strand is the same
    my $same_chr_strand = ($chr1 eq $chr2 and $strand1 eq $strand2) ? 1 : 0;

    # Check for common and easy case: the same alignment in both true case and test case
    if($same_chr_strand and $start1 == $start2 and $cigar1 eq $cigar2) {
	$base_right_count = $READ_LEN;
	$base_wrong_count = 0;
	# The most common case: ungapped alignment
	if($cigar1 eq $READ_LEN.'M') {
	    $introns1 = $introns2 = [];
	}
	# Gapped alignment
	else {
	    # Get the genomic coordinate for each position in the read
	    # We only need to compute this for one of the alignments, since they are identical
	    my ($del, $ins);
	    ($coords1, $introns1, $del, $ins) = get_coords($chr1, $start1, $cigar1, $del_shifts, $ins_shifts);
	    $introns2 = $introns1;
	    $coords2 = $coords1;
	    # Compare indels
	    compare_del($read_id1, $del, $del, $same_chr_strand, $gap_stats);
	    compare_ins($read_id1, $ins, $ins, $same_chr_strand, $gap_stats);
	}
    }
    # Otherwise we need to do some more work
    else {
	# Get the genomic coordinate for each position in the read
	my ($del1, $del2, $ins1, $ins2);
	($coords1, $introns1, $del1, $ins1) = get_coords($chr1, $start1, $cigar1, $del_shifts, $ins_shifts);
	($coords2, $introns2, $del2, $ins2) = get_coords($chr2, $start2, $cigar2);
	# Adjust true coords for indel ambguity (only needed if the alignments coincide)
	$coords1 = shift_coords($coords1, $del1, $ins1) if($same_chr_strand);
	# Compare base placement and indels
	($base_right_count, $base_wrong_count) = compare_basewise($coords1, $coords2, $same_chr_strand);
	compare_del($read_id1, $del1, $del2, $same_chr_strand, $gap_stats);
	compare_ins($read_id1, $ins1, $ins2, $same_chr_strand, $gap_stats);
    }

    # Get slot for read-level and base-level stats 
    # We count spliced and contiguous reads separately
    my $read_stats = $stats->{ @$introns1 == 0 ? 'contig' : 'spliced' };

    # Increment total reads and alignments seen
    $read_stats->[ READS_TOTAL ]++;
    $read_stats->[ ALN_TOTAL ]++ if($chr2 ne "");

    # Record read placement stats
    # (we could also make this quantitative, recording the distribution of number of correctly placed bases per alignment)
    if($base_right_count == $READ_LEN) {
	$read_stats->[ ALN_PERFECT ]++;  # Perfect alignment
    }
    elsif($base_right_count > 0) {
	$read_stats->[ ALN_1_RIGHT ]++;  # At least 1 base correct
    }
    elsif($same_chr_strand and compare_overlap($coords1, $coords2)) {
	$read_stats->[ ALN_OVERLAP ]++;  # Overlaps true alignment but no correctly placed base (unusual case)
    }

    # Record number of correct and wrong bases
    $read_stats->[ BASE_RIGHT ] += $base_right_count;
    $read_stats->[ BASE_WRONG ] += $base_wrong_count;
     
    # Intron comparison
    compare_introns($introns1, $introns2, $chr1, $chr2, $strand1, $strand2, $junc_coverage, $gap_stats);

    # Report info about each alignment?
    # if($REPORT_ALN) {
    # 	print STDERR join("\t", "@$read_id1",
    # 			  $chr1, $start1, $cigar1, ($strand1?'-':'+'),
    # 			  $chr2, $start2, $cigar2, ($strand2?'-':'+'),
    # 			  $base_right_count,
    # 			  $base_wrong_count), "\n";
    # }

    # In verbose mode, report some results
    # if($VERBOSE) {
    # 	if($base_wrong_count > 0) {
    # 	    print STDERR "Read $read_id1 has $base_wrong_count bases wrong (result=$chr2:$start2:$cigar2, truth=$chr1:$start1:$cigar1)\n";
    # 	}
    # }

    # Return intron strings (to be used for intron evaluation at the read pair level)
    return(get_intron_strings($chr1, $strand1, $introns1),
	   get_intron_strings($chr2, $strand2, $introns2));
}


# Compare coords between two alignments
sub compare_basewise {
    my ($coords1, $coords2, $same_chr_strand) = @_;

    my ($right, $wrong) = (0, 0);

    if($same_chr_strand) {
	for my $i (0..($READ_LEN-1)) {
	    if(defined($coords2->[$i])) { # Only count aligned parts of the read
		# Count this base in stats element 0 (correct) or 1 (wrong)
		if(find_in_array($coords2->[$i], $coords1->[$i])) {
		    $right++;
		}
		else {
		    $wrong++;
		}
	    }
	}
    }
    else {
	for my $i (0..($READ_LEN-1)) {
	    if(defined($coords2->[$i])) { # Only count aligned parts of the read
		$wrong++;
	    }
	}
    }

    return ($right, $wrong);
}


sub compare_overlap {
    my ($coords1_nested, $coords2) = @_;

    # Linearize and sort coords1
    # (coords2 should already be linear and sorted)
    my @coords1 = sort { $a <=> $b } map { @$_ } @$coords1_nested;

    # Walk through arrays checking for overlap
    my $c1 = shift @coords1;
    foreach my $c2 (@$coords2) {
	next unless(defined $c2);
	while($c1 < $c2) {  # Advance coords1 until we reach current position in coords2
	    $c1 = shift @coords1;
	    return 0 unless(defined $c1);
	}
	if($c1 == $c2) {
	    return 1; # overlap found!
	} 
    }

    return 0; # no overlap
}


# Compare introns between two alignments
sub compare_introns {
    my ($introns1, $introns2, $chr1, $chr2, $strand1, $strand2, $junc_coverage, $stats) = @_;

    my %correct;

    foreach my $i1 (@$introns1) {
	# $i1 contains the true intron data: [ query pos, target pos, length ]
	# To keep things simple we don't require the query pos to match.
	# Look for an identical intron in the aligner results
	my $match;
	if($chr1 eq $chr2 and $strand1 eq $strand2) {
	    foreach my $i2 (@$introns2) {
		if($i1->[1] == $i2->[1] and $i1->[2] == $i2->[2]) {
		    $correct{$i2} = 1;
		    $match = 1;
		    last;
		}
	    }
	}
	# Count as correct or missed depending on whether a match was found
	my $read_pos = $i1->[0];
	$read_pos = $READ_LEN - $read_pos - 2 if($strand1); # read position is zero-based, so subtract 2
	my $cov = $junc_coverage->{"$chr1:$i1->[1]:$i1->[2]"};
	unless($cov) { die "ERROR: no coverage for junction $chr1:$i1->[1]:$i1->[2]\n" }
	if($match) {
	    $stats->[ INTRON_CORRECT ][ $read_pos ][ $cov ]++;
	} else {
	    $stats->[ INTRON_MISSED ][ $read_pos ][ $cov ]++;
	}
    }

    # Count false positive introns
    foreach my $i2 (@$introns2) {
	next if($correct{$i2});
	my $read_pos = $i2->[0];
	$read_pos = $READ_LEN - $read_pos - 2 if($strand2); # read position is zero-based, so subtract 2
	my $cov = $junc_coverage->{"$chr2:$i2->[1]:$i2->[2]"} || 0;
	$stats->[ INTRON_FALSE ][ $read_pos ][ $cov ]++;
    }

    # Multi-intron counts
    $stats->[ MI_READ_ACTUAL ][ scalar @$introns1 ]++;
    $stats->[ MI_READ_REPORTED ][ scalar @$introns2 ]++;
    $stats->[ MI_READ_CORRECT ][ scalar(keys %correct) ]++;
}


# Gather multi-intron stats for read pairs
# Input:
#  1. True introns for read 1 (strings) 
#  2. Reported introns for read 1 (strings)
#  3. True introns for read 2 (strings) 
#  4. Reported introns for read 2 (strings)
#  5. Array ref for storing result
sub calc_multi_intron_read_pair_stats {
    my ($truth1, $pred1, $truth2, $pred2, $stats) = @_;

    my %index_truth = map { $_ => 1 } (@$truth1, @$truth2);
    my %index_pred = map { $_ => 1 } (@$pred1, @$pred2);

    my $correct = 0;
    foreach my $intron (keys %index_truth) {
	$correct++ if($index_pred{$intron});
    }

    $stats->[ MI_PAIR_ACTUAL ][ scalar(keys %index_truth) ]++;
    $stats->[ MI_PAIR_REPORTED ][ scalar(keys %index_pred) ]++;
    $stats->[ MI_PAIR_CORRECT ][ $correct ]++;
}


sub get_intron_strings {
    my ($chr, $strand, $introns) = @_;
    my @strings;
    foreach my $intron (@$introns) {
	push @strings, "$chr:$intron->[1]:$intron->[2]";
    }
    return \@strings;
}


# Compare deletions between two alignments
sub compare_del {
    my ($read_id, $del1, $del2, $same_chr_strand, $stats) = @_;

    my %correct;

    foreach my $d1 (@$del1) {
	# Get position of the true deletion, accounting for ambiguity
	my $query_pos_min = $d1->[0] - $d1->[3]; # query pos - left shift
	my $query_pos_max = $d1->[0] + $d1->[4]; # query pos + right shift
	my $target_pos_min = $d1->[1] - $d1->[3]; # target pos - left shift
	my $target_pos_max = $d1->[1] + $d1->[4]; # target pos + right shift
	# Look for an identical deletion in the aligner results
	my $match;
	if($same_chr_strand) {
	    foreach my $d2 (@$del2) {
		if(!$correct{$d2} and # make sure we only count each del once
		   $d2->[0] >= $query_pos_min and $d2->[0] <= $query_pos_max and
		   $d2->[1] >= $target_pos_min and $d2->[1] <= $target_pos_max and
		   $d1->[2] == $d2->[2]) {
		    $correct{$d2} = 1;
		    $match = 1;
		    last;
		}
	    }
	}
	# Count as correct or missed depending on whether a match was found
	# But don't count deletions at starts/ends of reads as they are subject
	# to formatting ambiguity
	if($query_pos_min >= 0 and $query_pos_max < $READ_LEN-1) {
	    if($match) {
		$stats->[ DEL_CORRECT ][ $d1->[2] ]++;
	    } else {
		$stats->[ DEL_MISSED ][ $d1->[2] ]++;
	    }
	}
    }

    # Count false positive deletions
    # Again, don't count deletions at starts or ends of reads
    foreach my $d2 (@$del2) {
	next if($correct{$d2} or $d2->[0] == -1 or $d2->[0] == $READ_LEN-1);
	$stats->[ DEL_FALSE ][ $d2->[2] ]++;
	if($d2->[2] >= $REPORT_FALSE_INDELS) {
	    print STDERR "False deletion of size $d2->[2] in read $read_id\n";
	}
    }
}


# Compare insertions between two alignments
sub compare_ins {
    my ($read_id, $ins1, $ins2, $same_chr_strand, $stats) = @_;

    my %correct;

    foreach my $i1 (@$ins1) {
	# Get position of the true insertion, accounting for ambiguity
	my $query_pos_min = $i1->[0] - $i1->[3]; # query pos - left shift
	my $query_pos_max = $i1->[0] + $i1->[4]; # query pos + right shift
	my $target_pos_min = $i1->[1] - $i1->[3]; # target pos - left shift
	my $target_pos_max = $i1->[1] + $i1->[4]; # target pos + right shift
	# Look for an identical insertion in the aligner results
	my $match;
	if($same_chr_strand) {
	    foreach my $i2 (@$ins2) {
		if(!$correct{$i2} and # make sure we only count each insertion once
		   $i2->[0] >= $query_pos_min and $i2->[0] <= $query_pos_max and
		   $i2->[1] >= $target_pos_min and $i2->[1] <= $target_pos_max and
		   $i1->[2] == $i2->[2]) {
		    $correct{$i2} = 1;
		    $match = 1;
		    last;
		}
	    }
	}
	# Count as correct or missed depending on whether a match was found
	# But don't count insertions at starts/ends of reads as they are subject
	# to formatting ambiguity
	if($query_pos_min > 0 and $query_pos_max + $i1->[2] < $READ_LEN) {
	    if($match) {
		$stats->[ INS_CORRECT ][ $i1->[2] ]++;
	    } else {
		$stats->[ INS_MISSED ][ $i1->[2] ]++;
	    }
	}
    }

    # Count false positive insertions
    foreach my $i2 (@$ins2) {
	next if($correct{$i2} or $i2->[0] <= 0 or $i2->[0] + $i2->[2] >= $READ_LEN); # It should be sufficient to use == instead of <= and >= here
	$stats->[ INS_FALSE ][ $i2->[2] ]++;
	if($i2->[2] >= $REPORT_FALSE_INDELS) {
	    print STDERR "False insertion of size $i2->[2] in read $read_id\n";
	}
    }
}


sub get_coords {
    my ($chr, $pos, $cigar, $del_shifts, $ins_shifts) = @_;

    my (@coords, @introns, @deletions, @insertions);

    # Parse CIGAR string to get genomic coord for each nt in the read.
    # Also collect genomic positions of indels.
    my @cig_fields = split(/([A-Z])/, $cigar);
    while(@cig_fields) {
	my $l = shift @cig_fields;
	my $op = shift @cig_fields;
	if($op eq 'M') {
	    while($l) {
		push @coords, $pos;
		$pos++;
		$l--;
	    }
	}
	elsif($op eq 'N' or $op eq 'D') {
	    if($l > MAX_DEL) {  # Count as intron or deletion depending on size
		# Store read pos before intron, target pos where intron starts, and intron size
		# We ignore terminal introns as they are subject to formatting ambiguity
		push @introns, [$#coords, $pos, $l] unless(@coords == 0 or @coords == $READ_LEN);
		$pos += $l;
	    } else {
		# Store read pos before deletion, target pos where deletion starts, deletion size and shifts
		my @shifts = get_del_shifts($del_shifts, $chr, $pos, $l) if($del_shifts);
		push @deletions, [$#coords, $pos, $l, @shifts];
		$pos += $l;
	    }
	}
	elsif($op eq 'I') {
	    # Store read pos where insertion starts, target pos after insertion, insertion size and shifts
	    my @shifts;
	    if($ins_shifts) {
		my $is_terminal = (@coords == 0 or @cig_fields == 0) ? 1 : 0; 
		@shifts = get_ins_shifts($ins_shifts, $chr, $pos, $l, $is_terminal)
	    }
	    push @insertions, [$#coords+1, $pos, $l, @shifts];
	    while($l) {
		push @coords, $pos - .5;
		$l--;
	    }    
	}
	elsif($op eq 'S') {
	    while($l) {
		push @coords, undef;
		$l--;
	    }
	}
	else {
	    die "I don't know what to do with CIGAR operation $op";
	}
    }

    return (\@coords, \@introns, \@deletions, \@insertions);
}


# shift_coords()

# This function takes care of the coordinate shifting magic. The aim
# is to take indel ambiguity into account, by adding additional genomic
# coordinates for each read position accordingly.

# Following the function code are some examples to illustrate the
# calculations.

sub shift_coords {
    my ($coords_in, $deletions, $insertions) = @_;

    # Make an array reference for each read position.  This allows us
    # to store multiple alternative genomic coordinates per read
    # position. Initialize the arrays to contain the input
    # coordinates.
    my @coords = map { [ $_ ] } @$coords_in;
   
    # Consider deletion shifts
    foreach my $del (@$deletions) { # For each deletion in the alignment
	my ($del_query_pos, $del_target_pos, $del_size, $left_shift, $right_shift) = @$del;
	# Compute left shifted coords
	for my $i (1 .. $left_shift) {
	    my $j = $del_query_pos - $i +  1;
	    last if($j < 0);
	    my $new_target_pos = $del_target_pos + $del_size - $i;
	    push @{$coords[$j]}, $new_target_pos;
	}
	# Compute right shifted coords
	for my $i (1 .. $right_shift) {
	    my $j = $del_query_pos + $i;
	    last if($j > $#coords);
	    my $new_target_pos = $del_target_pos + $i - 1;
	    push @{$coords[$j]}, $new_target_pos;
	}
    }	

    # Consider insertion shifts
    foreach my $ins (@$insertions) { # For each insertion in the alignment
	my ($ins_query_pos, $ins_target_pos, $ins_size, $left_shift, $right_shift) = @$ins;
	# Compute left shifted coords
	for my $i (1 .. $left_shift) {
	    my $j_min = $ins_query_pos - $i;
	    my $j = $j_min + $ins_size;
	    last if($j < 0);
	    $j_min = 0 if($j_min < 0);
	    push @{$coords[$j]}, $ins_target_pos - $i;
	    while($j > $j_min) {
		$j--;
		push @{$coords[$j]}, $ins_target_pos - $i - .5;
	    }
	}
	# Compute right shifted coords
	$ins_query_pos--; # Query position before the insertion
	$ins_target_pos--; # Genomic position before the insertion
	for my $i (1 .. $right_shift) {
	    my $j = $ins_query_pos + $i;
	    my $j_max = $j + $ins_size;
	    last if($j > $#coords);
	    $j_max = $#coords if($j_max > $#coords);
	    push @{$coords[$j]}, $ins_target_pos + $i;
	    while($j < $j_max) {
		$j++;
		push @{$coords[$j]}, $ins_target_pos + $i + .5;
	    }
	}
    }
    
    return \@coords;
}

# Examples to illustrate the calculations

# DELETION SHIFT
# --------------

# 01   2345  
# AA...TTAA  query (read)
# AATTTTTAA  target (genome)
# 123456789

# Deletion: TTT
# del_query_pos = 1
# del_target_pos = 3
# Can be shifted right by 2.

# Shift 1:

# 012   345
# AAT...TAA
# AATTTTTAA
# 123456789

# Shift 2:

# 0123   45
# AATT...AA
# AATTTTTAA
# 123456789

# For left shift, consider the reverse scenario.
# We start from the last example above.
# del_query_pos = 3
# del_target_pos = 5
# Can be shifted left by 2.

# INSERTION SHIFT
# ---------------

# 012345678
# AATTTTTAA  query (read)
# AA...TTAA  target (genome)
# 122223456  The insertion between pos 2 and 3 is nubmered 2.5
#   555

# Insertion: TTT
# ins_query_pos = 2
# ins_target_pos = 3
# Can be shifted right by 2

# Shift 1:

# 012345678
# AATTTTTAA
# AAT...TAA
# 123333456
#    555 

# Shift 2:

# 012345678
# AATTTTTAA
# AATT...AA
# 123444456
#     555

# For left shift, consider the reverse scenario.
# We start from the last example above.
# ins_query_pos = 4
# ins_target_pos = 5
# Can be shifted left by 2.


sub is_aligned {
    my $sam = shift;
    my (undef, $flags) = split /\t/, $sam;
    return $flags & 0x4 ? 0 : 1;
}


sub read_indel_shifts {
    my $fn = shift;

    my $in = must_open_infile($fn);

    my (%del_shifts, %ins_shifts);

    while(my $line = <$in>) {
	chomp $line;
	my ($chr, $pos, $size, $left_shift, $right_shift) = split "\t", $line;
	if($size < 0) {
	    $size = -$size;
	    push @{$del_shifts{"$chr:$pos:$size"}}, [ $left_shift, $right_shift ];
	}
	else {
	    push @{$ins_shifts{"$chr:$pos"}}, [ $left_shift, $right_shift, $size ];
	}
    }

    close($in) or die "Failed to close file $fn: $!";

    return(\%del_shifts, \%ins_shifts);
}


# Get the left and right shifts allowed for a deletion.
# There may be some (very rare) cases with multiple shifts
# entries with the same pos; use the max.
sub get_del_shifts {
    my ($del_shifts, $chr, $del_target_pos, $del_size) = @_;
    my $left_shift = 0;
    my $right_shift = 0;
    my $shifts_array = $del_shifts->{"$chr:$del_target_pos:$del_size"};
    if($shifts_array) {
	foreach my $shifts (@$shifts_array) { 
	    $left_shift = $shifts->[0] if($left_shift < $shifts->[0]);
	    $right_shift = $shifts->[1] if($right_shift < $shifts->[1]);
	}
    }
    return ($left_shift, $right_shift);
}


# Get the left and right shifts allowed for an insertion.
# There may be some (very rare) cases with multiple shifts
# entries with the same pos; use the max.
# We require the insertion size to match excatly except in
# cases where the insertion is at the start or end of the
# alignment, and may therefore be partial.
sub get_ins_shifts {
    my ($ins_shifts, $chr, $ins_target_pos, $ins_size, $is_terminal) = @_;
    my $left_shift = 0;
    my $right_shift = 0;
    my $shifts_array = $ins_shifts->{"$chr:$ins_target_pos"};
    if($shifts_array) {
	foreach my $shifts (@$shifts_array) { 
	    if($shifts->[2] == $ins_size or # Exact size match (ideal situtation)
	       ($shifts->[2] > $ins_size and $is_terminal)) { # Potentially partial insertion (OK if at start or end of alignment)
		$left_shift = $shifts->[0] if($left_shift < $shifts->[0]);
		$right_shift = $shifts->[1] if($right_shift < $shifts->[1]);
	    }
	}
    }
    return ($left_shift, $right_shift);
}


sub read_junction_coverage {
    my ($fn) = @_;

    my %index;
    my $in = must_open_infile($fn);

    while(my $line = <$in>) {
	chomp $line;
	my ($intron_id, $coverage) = split "\t", $line;
	$coverage = MAX_JUNC_COV if($coverage > MAX_JUNC_COV);
	my ($chr, $start, $end) = split /[:\-]/, $intron_id;
	$start++;
	my $size = $end - $start;
	$index{ "$chr:$start:$size" } = $coverage;
    }

    close $in;

    return \%index;
}


sub find_in_array {
    my ($value, $array) = @_;
    foreach my $array_value (@$array) {
	return 1 if($value == $array_value);
    }
    return 0;
}


sub max_array_size {
    return max( map { $#$_ } @_ );
}


sub must_open_aln_file {
    my $fn = shift;
    if($fn =~ /\.bam$/i) {
	return must_open_infile("samtools view $fn |");
    }
    elsif($fn =~ /\.sam$/i) {
	return must_open_infile($fn);
    }
    else {
	die "Alignment file extension must be .sam or .bam: $fn\n";
    }
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

