#!/bin/env perl

# Evaluate transcript reconstruction performance

use warnings;
use strict;
use List::Util qw(sum);
use Getopt::Std;
use Genoman::DataStore::Binner;
use Genoman::BlockSet;

# Globals

my $DEF_UNSPLICED_COVERAGE = 0.6;
my ($DEBUG, $UNSPLICED_STRANDED, $SPLICED_STRANDED, $UNSPLICED_COVERAGE);

my $BINNER = Genoman::DataStore::Binner->new;

# Definitions

use Class::Struct PredStats =>
    [ ref_tx => '@',
      ref_exons => '@',
      pred_tx_count => '$',
      pred_tx_corr => '$',
      pred_tx_part => '$',
      pred_tx_wrong => '$',
      pred_exon_count => '$',
      pred_exon_corr => '$',
      pred_exon_corr_known => '$'
    ];

use constant {
    TX_HIT_WRONG => 1,
    TX_HIT_PART => 2,
    TX_HIT_CORR => 3,
    EXON_SINGLE => 0,
    EXON_FIRST => 1,
    EXON_INTERNAL => 2,
    EXON_LAST => 3
};


# Do the job and exit
main();
exit(0);


# Function to report usage information
sub usage {
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Evaluate transcript reconstruction performance

Usage: perl $cmd [options] truth.gtf[.gz] outPrefix inputDir1 [inputDir2 ...]

truth.gtf   Correct transcript models
outPrefix   Prefix for output files
inputDir    Cufflinks result directory containing transcripts.gtf

Options:
-a file    GTF file with annotation provided to aligners (can be gzipped)
-d         Debug mode
-s         Treat single-exon transcripts as unstranded.
-S         Treat multi-exon transcripts and their constituent exons as unstranded.
-t string  Only consider reference transcripts of this type (from truth.gtf)
-u float   Coverage of single exons required. Default: $DEF_UNSPLICED_COVERAGE.

EndOfUsage
    exit(1);
}


# Main function
sub main {

    # Parse args
    my %args;
    getopts('a:dsSt:u:', \%args);
    my $known_fn = $args{'a'};
    $DEBUG = $args{'d'};
    $UNSPLICED_STRANDED = $args{'s'} ? 0 : 1;
    $SPLICED_STRANDED = $args{'S'} ? 0 : 1;
    my $ref_tx_type = $args{'t'};
    $UNSPLICED_COVERAGE = defined($args{'u'}) ? $args{'u'} : $DEF_UNSPLICED_COVERAGE;
    @ARGV >= 3 or usage();
    my ($truth_fn, $out_prefix, @in_dir_list) = @ARGV;

    # Read annotation
    my $ref_tx = read_gtf($truth_fn, $ref_tx_type);

    # Index annotation
    my ($tx_index, $splice_index, $ref_tx_groups) = index_tx($ref_tx);
    my ($exon_index, $ref_exons) = index_exons($ref_tx);

    # Classify reference transcripts as known or novel
    if($known_fn) {
	my $known_tx = read_gtf($known_fn);
	report_ref_tx_novelty($out_prefix, $tx_index, $known_tx);
	set_exon_novelty($exon_index, $ref_exons, $known_tx);
    }

    # Process each prediction file
    my @stats;
    foreach my $in_dir (@in_dir_list) {
	print STDERR "Processing results from $in_dir\n";
	my $pred_fn = "$in_dir/transcripts.gtf";
	my $pred_tx;
	if(-e $pred_fn) {
	    $pred_tx = read_gtf($pred_fn);
	} else {
	    warn "Warning: input file $pred_fn does not exist\n";
	    $pred_tx = [];
	}
	push @stats, cmp_pred_to_ref($pred_tx, $ref_exons, $tx_index, $splice_index, $exon_index);
    }

    # Output result
    print STDERR "Reporting\n";
    write_report($out_prefix, \@in_dir_list, $ref_tx_groups, $ref_exons, \@stats);
}



sub write_report {
    my ($out_prefix, $pred_id_array, $ref_tx_groups, $ref_exons, $stats) = @_;

    # Print mapping of splice patterns to transcript IDs
    my $fn = "$out_prefix.tx_groups.txt";
    open OUT, ">$fn" or die "Failed to open $fn for output";
    for my $i (0..$#$ref_tx_groups) {
	foreach my $tx (@{$ref_tx_groups->[$i]}) {
	    print OUT $i+1, "\t", $tx->[0], "\n";
	}
    }
    close OUT;

    # Print recall matrix for reference splice patterns
    # A pattern is an intron combination and a strand
    # Row: pattern. Column: protocol. Value: identification status (integer)
    $fn = "$out_prefix.tx_match.txt";
    open OUT, ">$fn" or die "Failed to open $fn for output";
    print OUT join("\t", "submission", @$pred_id_array), "\n";
    for my $i (0..$#$ref_tx_groups) {
	print OUT $i+1;
	for my $j (0..$#$pred_id_array) {
	    print OUT "\t", $stats->[$j]->ref_tx($i);
	}
	if($DEBUG) {
	    my $tx = $ref_tx_groups->[$i][0];
	    print OUT "\t# $tx->[0] $tx->[2]:$tx->[4][0]-$tx->[4][-1]";
	}
	print OUT "\n";
    }
    close OUT;

    # Write stats file
    $fn = "$out_prefix.stats.txt";
    open OUT, ">$fn" or die "Failed to open $fn for output";
    print OUT join("\t", "submission",
		   "ref.tx", "ref.tx.corr", "ref.tx.part", "ref.tx.wrong",
		   "pred.tx", "pred.tx.corr", "pred.tx.part", "pred.tx.wrong",
		   "ref.exons", "ref.exons.known", "ref.exons.corr", "ref.exons.corr.known",
		   "pred.exons", "pred.exons.corr", "pred.exons.corr.known"), "\n";
    for my $i (0..$#$pred_id_array) {
	my $s = $stats->[$i];
	my @ref_tx_counts = count_integers( $s->ref_tx(), TX_HIT_CORR );
	my ($correct_ref_exons, $correct_ref_exons_known, $ref_exons_known) = count_correct_exons( $s->ref_exons(), $ref_exons );
	print OUT join("\t", $pred_id_array->[$i],
		       scalar(@$ref_tx_groups), $ref_tx_counts[TX_HIT_CORR], $ref_tx_counts[TX_HIT_PART], $ref_tx_counts[TX_HIT_WRONG],
		       $s->pred_tx_count(), $s->pred_tx_corr(), $s->pred_tx_part(), $s->pred_tx_wrong(),
		       scalar(@$ref_exons), $ref_exons_known, $correct_ref_exons, $correct_ref_exons_known,
		       $s->pred_exon_count(), $s->pred_exon_corr(), $s->pred_exon_corr_known()), "\n";
    }
    close OUT;

}


sub count_correct_exons {
    my ($is_correct, $ref_exons) = @_;

    my ($correct, $correct_known, $known) = (0, 0, 0);
    
    unless(@$is_correct == @$ref_exons) { die "ERROR: array size mismatch" }

    for my $i (0..$#$is_correct) {
	if($is_correct->[$i]) {
	    $correct++;
	    $correct_known++ if($ref_exons->[$i][5]);
	}
	$known++ if($ref_exons->[$i][5]);
	if($DEBUG) {
	    my ($type, $chr, $strand, $start, $end, $is_known) = @{$ref_exons->[$i]};
	    print STDERR "known=$is_known correct=$is_correct->[$i] type=$type $chr:$start-$end $strand\n";
	}
   }

    return ($correct, $correct_known, $known);
}


sub report_ref_tx_novelty {
    my ($out_prefix, $ref_splicing_index, $known_tx) = @_;

    my @is_known = (0) x scalar(keys %$ref_splicing_index);

    # For each known splice pattern, check if it is represented by a group of reference transcripts
    # If so, label that group of reference transcripts as "known" 
    foreach my $tx (@$known_tx) {
	my (undef, undef, $chr, $strand, $exons) = @$tx;
	if(@$exons > 2) {
	    my $splicing_key = make_splicing_key($exons);
	    my @target_strands = $SPLICED_STRANDED ? $strand : ("+", "-");
	    foreach my $target_strand (@target_strands) {
		my $i = $ref_splicing_index->{"$chr:$splicing_key:$target_strand"};
		$is_known[$i] = 1 if(defined $i);
	    }	    
	}
    }

    my $fn = "$out_prefix.tx_known.txt";
    open OUT, ">$fn" or die "Failed to open $fn for output";
    print OUT "id\tknown\n";
    for my $i (0..$#is_known) {
	print OUT $i+1, "\t", $is_known[$i], "\n"; 
    }
    close OUT;
}


sub cmp_pred_to_ref {
    my ($pred_tx, $ref_exons, $tx_index, $splice_index, $exon_index) = @_;

    my @pred_tx_counts = (0, 0, 0, 0);
    my ($pred_exon_count, $correct_pred_exon_count, $correct_pred_exon_count_known) = (0, 0, 0);
    my @correct_ref_tx_array = (0) x scalar(keys %$splice_index);
    my @correct_ref_exon_array = (0) x scalar(@$ref_exons);
    my (%seen_tx, %seen_exons);

    # Check each predicted transcript
    foreach my $t1 (@$pred_tx) {

	my (undef, undef, $chr, $strand1, $exons1) = @$t1;

	# We ignore unstranded predictions entirely if options are set to require strand agreement for all features
	next if($SPLICED_STRANDED and $UNSPLICED_STRANDED and $strand1 eq '.');

	# If spliced, compare against reference spliced transcripts
	# We ignore unstranded predictions if options are set to require strand agreement for transcripts
	if(@$exons1 > 2 and !($SPLICED_STRANDED and $strand1 eq '.')) {

	    # Create key that captures intron coords
	    my $splicing_key = make_splicing_key($exons1);
	    
	    # Evaluate transcript unless we have seen the same splice pattern before
	    unless($seen_tx{"$chr:$splicing_key:$strand1"}) {

		$seen_tx{"$chr:$splicing_key:$strand1"} = 1;

		# Find reference transcripts with one or more matching introns
		my $matched = 0;
		my @target_strands = $SPLICED_STRANDED ? $strand1 : ("+", "-");
		foreach my $strand (@target_strands) {

		    #print STDERR "$chr $strand ", $exons1->as_string, "\n";
			
		    my $splice_hits = splice_search($splice_index, $chr, $strand, $exons1);

		    foreach my $hit (@$splice_hits) {

			#print STDERR " hit: $hit\n";

			my $hit_category = 0;
			if($hit eq "$chr:$splicing_key:$strand") {
			    # "Correct": prediction and reference splice patterns match precisely (all introns match)
			    $hit_category = TX_HIT_CORR;
			}
			elsif($hit =~ /$splicing_key/) {
			    # "Partial": prediction is a part of reference splice pattern
			    $hit_category = TX_HIT_PART;
			}
			else {
			    # "Incorrect": one or more junctions match, but patterns are incompatible or
			    # reference is part of prediction
			    $hit_category = TX_HIT_WRONG;
			}

			if($hit_category) {
			    my $t2 = $tx_index->{$hit};
			    $correct_ref_tx_array[$t2] = $hit_category if($correct_ref_tx_array[$t2] < $hit_category);
			    $matched = $hit_category if($matched < $hit_category);
			}
		    }
		}

		# Count predicted transcript
		$pred_tx_counts[$matched]++;

	    }
	}

	# Evaluate each exon
	for(my $i = 0; $i < @$exons1; $i += 2) {

	    # Get exon coords 
	    my $start1 = $exons1->[$i];
	    my $end1 = $exons1->[$i+1];

	    # Check that we have not seen this exon before
	    my $exon_key = "$chr:$start1:$end1:$strand1";
	    next if($seen_exons{$exon_key});
	    $seen_exons{$exon_key} = 1;

	    # Find hits in bin index (unstranded)
	    my @hits = search_bin_index($exon_index, $chr, $start1, $end1);

	    # Evaluate depending on type of hit exon
	    my ($matched_novel, $matched_known);
	    foreach my $hit (@hits) {
		my ($exon_type, undef, $strand2, $start2, $end2, $is_known) = @{$ref_exons->[$hit]};
		if(compare_exons($start1, $end1, $strand1, $start2, $end2, $strand2, $exon_type)) {
		    if($is_known) { $matched_known = 1 } else { $matched_novel = 1 }
		    $correct_ref_exon_array[$hit] = 1;
		}
	    }

	    # Count predicted exon
	    $pred_exon_count++;
	    $correct_pred_exon_count++ if($matched_novel or $matched_known);
	    $correct_pred_exon_count_known++ if($matched_known);
	}
    }

    return PredStats->new(ref_tx => \@correct_ref_tx_array,
			  ref_exons => \@correct_ref_exon_array,
			  pred_tx_count => sum(@pred_tx_counts),
			  pred_tx_wrong => $pred_tx_counts[TX_HIT_WRONG],
			  pred_tx_part => $pred_tx_counts[TX_HIT_PART],
			  pred_tx_corr => $pred_tx_counts[TX_HIT_CORR],
			  pred_exon_count => $pred_exon_count,
			  pred_exon_corr => $correct_pred_exon_count,
			  pred_exon_corr_known => $correct_pred_exon_count_known);
}


sub compare_exons {
    my ($start1, $end1, $strand1, $start2, $end2, $strand2, $exon2_type) = @_;

    my $match = 0;

    if($exon2_type == EXON_SINGLE) {
	if(!$UNSPLICED_STRANDED or $strand1 eq $strand2) {  # Check strand agreement
	    # Check for sufficient overlap (proportion given by $UNSPLICED_COVERAGE)
	    my $left = $start1 > $start2 ? $start1 : $start2;
	    my $right = $end1 < $end2 ? $end1 : $end2;
	    my $cutoff = $UNSPLICED_COVERAGE * ($end2 - $start2);
	    $match = 1 if($right - $left + 1 >= $cutoff);
	}
    }
    elsif(!$SPLICED_STRANDED or $strand1 eq $strand2) {    # Check strand agreement
	# Check for coordinate agreement, according to hit exon type
	if($exon2_type == EXON_INTERNAL) { # Internal exon: perfect agreement required
	    $match = 1 if($start1 == $start2 and $end1 == $end2);
	}
	elsif($exon2_type == EXON_FIRST) { # First exon: internal border agreement required
	    $match = 1 if($end1 == $end2);
	}
	elsif($exon2_type == EXON_LAST) { # Last exon: internal border agreement required
	    $match = 1 if($start1 == $start2);
	}
	else {
	    die "Illegal exon type: $exon2_type";
	}
    }
    
    return $match;
}



## index_tx()
## This function assigns a numerical ID to each splice pattern.
## Splice pattern = unique combination of introns and strand.
## This mapping from splice pattern to ID is stored in %splice_index.
## The mapping from ID to transcripts is stored in @spliced_tx,
## such that $spliced_tx[$id] = [ $tx1, $tx2, ... ],
## where each $txN is an array with transcript info.

sub index_tx {
    my $tx_array = shift;
    
    my (%splicing_index, %splice_index, @spliced_tx);

    my $i = 0;
    foreach my $tx (@$tx_array) {
	my (undef, undef, $chr, $strand, $exons) = @$tx;
	if(@$exons > 2) {

	    # Create "splicing key" (a string that describes the splicing pattern)
	    my $splicing_key = make_splicing_key($exons);
	    $splicing_key = "$chr:$splicing_key:$strand";

	    # Index the splicing pattern if not already done
	    my $j = $splicing_index{$splicing_key};
	    unless(defined $j) {
		# %splicing_index:
		#  key: splicing key
		#  value: index for splicing pattern
		$splicing_index{$splicing_key} = $j = $i;
		$i++;
		# %splice_index:
		#  key: string describing a single splice junction
		#  value: array of splicing keys for splice patterns containing this junction
		foreach my $splice_key (make_splice_keys($chr, $strand, $exons)) {
		    push @{$splice_index{$splice_key}}, $splicing_key;
		}
	    }

	    # Index transcript by splicing pattern
	    push @{$spliced_tx[$j]}, $tx;

	}
    }
    
    return (\%splicing_index, \%splice_index, \@spliced_tx);
}



sub index_exons {
    my $tx_array = shift;
    
    my $exon_index = {};
    my @exon_array;
    my %seen_exons;
    my $nr_exons = 0;  # Counter for unique exons

    # Loop over all transcripts
    foreach my $tx (@$tx_array) {

	# Extract exon coords for transcript
	my (undef, undef, $chr, $strand, $exons) = @$tx;

	# Loop over all exons in transcript
	for(my $i = 0; $i < @$exons; $i += 2) {

	    # Get coords of this exon
	    my $start = $exons->[$i];
	    my $end = $exons->[$i+1];

	    # Determine exon type
	    my $type;
	    if(@$exons == 2) {
		$type = EXON_SINGLE;
	    }
	    elsif($i == 0) {
		$type = EXON_FIRST;
	    }
	    elsif($i == @$exons-2) {
		$type = EXON_LAST;
	    }
	    else {
		$type = EXON_INTERNAL;
	    }

	    # Check that we have not seen an identical exon before
	    my $exon_key = "$type:$chr:$strand:$start:$end";
	    next if($seen_exons{$exon_key});
	    $seen_exons{$exon_key} = 1;

	    # Store this exon in the array and binned index
	    # Element 5 is a flag indicating whether the exon is known (will be set later)
	    $exon_array[$nr_exons] = [$type, $chr, $strand, $start, $end, 0];
	    bin_item($exon_index, $nr_exons, $chr, $start, $end);
	    $nr_exons++;
	}
    }

    return ($exon_index, \@exon_array);
}


sub set_exon_novelty {
    my ($exon_index, $ref_exons, $known_tx) = @_;

    my %seen_exons;

    # Loop over the known transcripts
    foreach my $t1 (@$known_tx) {

	my (undef, undef, $chr, $strand1, $exons1) = @$t1;

	# Evaluate each exon of the known transcript
	for(my $i = 0; $i < @$exons1; $i += 2) {

	    # Get exon coords 
	    my $start1 = $exons1->[$i];
	    my $end1 = $exons1->[$i+1];

	    # Check that we have not seen this exon before
	    my $exon_key = "$chr:$start1:$end1:$strand1";
	    next if($seen_exons{$exon_key});
	    $seen_exons{$exon_key} = 1;

	    # Find hits in bin index (unstranded)
	    my @hits = search_bin_index($exon_index, $chr, $start1, $end1);

	    foreach my $hit (@hits) {

		my ($exon_type, undef, $strand2, $start2, $end2, $is_known) = @{$ref_exons->[$hit]};

		next if($is_known); # If the exon is already classified as known, there is no need to check it again
		next if($strand1 ne $strand2); # Require strand agreement

		# Set each matching hit as known
		if(compare_exons($start1, $end1, $strand1, $start2, $end2, $strand2, $exon_type)) {
		    $ref_exons->[$hit][5] = 1;
		}
		
	    }
	    
	}
    }
}


sub bin_item {
    my ($index, $item, $chr, $start, $end) = @_;
    my $bin = $BINNER->bin_from_coord_range($start, $end);
    push @{$index->{$chr}{$bin}}, [$item, $start, $end];
}


sub search_bin_index {
    my ($index, $chr, $query_start, $query_end) = @_;
    $index = $index->{$chr} or return ();
    my @result;
    foreach my $bin_range ( $BINNER->bin_ranges_from_coord_range($query_start, $query_end) ) {
	foreach my $bin ($bin_range->[0] .. $bin_range->[1]) {
	    if($index->{$bin}) {
		foreach my $target ( @{$index->{$bin}} ) {
		    if($target->[1] <= $query_end and $target->[2] >= $query_start) {
			push @result, $target->[0];
		    }
		}
	    }
	}
    }
    return @result;
}


# Make string with splice junction coords, on the format
# :start1-end1:start2-end2:start3-end3:
# where start1 is the end coordinate of the first exon and end1 is the
# start coordinate of the second exon, and so on
sub make_splicing_key {
    my ($exons) = @_;
    my $key = ":";
    for(my $i = 1; $i < @$exons-2; $i += 2) {
	$key .= $exons->[$i] . '-' . $exons->[$i+1] . ':'; 
    }
    return $key;
}


sub make_splice_keys {
    my ($chr, $strand, $exons) = @_;
    my @keys;
    for(my $i = 1; $i < @$exons-2; $i += 2) {
	push @keys, $chr.':'.$exons->[$i].'-'.$exons->[$i+1].':'.$strand; 
    }
    return @keys;
}

   
# Given a (predicted) transcript structure, find all splicing patterns
# in the reference that contain introns from this structure.
# Returns an array of splice pattern keys.
sub splice_search {
    my ($splice_index, $chr, $strand, $exons) = @_;

    my %hits;

    foreach my $splice_key (make_splice_keys($chr, $strand, $exons)) {
	my $splice_pattern_array = $splice_index->{$splice_key};
	if(defined $splice_pattern_array) {
	    foreach my $pattern (@$splice_pattern_array) {
		$hits{$pattern} = 1;
	    }
	}
    }

    return [ keys %hits ];
}


# Read transcripts from GTF
sub read_gtf {
    my ($fn, $req_gene_type) = @_;

    # Data structure to store result
    my %tx_index;

    # Open file
    if($fn =~ /.gz$/i) {
	$fn = "zcat $fn |";
    }
    open GTF, $fn or die "Failed to open $fn for reading\n";

    # Read GTF file
    while(my $line = <GTF>) {
	# Parse GTF line
	chomp $line;
	my ($chr, $gene_type, $ft_type, $start, $end, undef, $strand, undef, $attrib) = split "\t", $line;
	next unless($ft_type eq 'exon'); # Only consider exons
	next if(defined($req_gene_type) and $req_gene_type ne $gene_type); # Only consider selected gene type if specified
        # Parse attribute string
	my ($tx_id) = $attrib =~ /transcript_id "([\w\.]+)";/;
	unless(defined $tx_id) { die "Failed to parse attribute string: $attrib" }
	# Update transcript info
	if(exists $tx_index{$tx_id}) {
	    $tx_index{$tx_id}[1] eq $gene_type or die "ERROR: conflicting gene types for transcript $tx_id in $fn";
	    $tx_index{$tx_id}[2] eq $chr or die "ERROR: confliciting chromosome for transcript $tx_id in $fn";
	    $tx_index{$tx_id}[3] eq $strand or die "ERROR: confliciting strand for transcript $tx_id in $fn";
	    push @{$tx_index{$tx_id}[4]}, ($start, $end);
	}
	else {
	    $tx_index{$tx_id} = [$tx_id, $gene_type, $chr, $strand, [$start, $end]];
	}
    }
    close GTF or die "Failed to close $fn";

    # Order exons
    foreach my $tx (values %tx_index) {
	$tx->[4] = Genoman::BlockSet->new_from_unsorted($tx->[4]);
    }

    # Return index
    return [values %tx_index];
}


sub count_integers {
    my ($array, $max) = @_;
    my @n;
    if(defined $max) {
	@n = (0) x ($max+1);
    }
    foreach my $x (@$array) {
	die "ERROR: $x > maximum allowed value $max" if($x > $max);
	$n[$x]++;
    }
    return @n;
}


sub count_true {
    my $array = shift;
    my $n = 0;
    foreach my $x (@$array) {
	$n++ if($x);
    }
    return $n;
}

