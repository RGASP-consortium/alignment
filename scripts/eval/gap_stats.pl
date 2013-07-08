#!/bin/env perl

# Script to collect stats about alignment gaps

# Potential TODO: If we would like to add edge distance calculation,
# as now separately implemented in intron_edge_distance.pl, we need to
# add a binary search index to the splices index structure (a fourth
# array element) The edge distances can then be counted and reported
# in make_splice_reports()

# Pragmas and includes
use warnings;
use strict;
use Getopt::Std;


# Constants
use constant {
    # CIGAR operations to count
    CIGAR_D => 0,
    CIGAR_I => 1,
    CIGAR_N => 2,
    CIGAR_D0 => 3,
    CIGAR_D1 => 4,
    CIGAR_I0 => 5,
    CIGAR_I1 => 6,
    CIGAR_S0 => 7,
    CIGAR_S1 => 8,
    # Intron classes
    KNOWN_INTRON => 0,
    SAME_GENE => 1,
    DIFF_GENES => 2,
    ONE_NOVEL => 3,
    BOTH_NOVEL => 4};

my @INTRON_CLASS = ("Known splice", "Same gene", "Different genes", "One novel", "Both novel");
my @CIGAR_OPS = ("D", "I", "N", "D0", "D1", "I0", "I1", "S0", "S1");

my $MAX_ALN_COUNT = 100; # Max alignments to count for splice stats
 # ^ keep this as a variable in case we make it a commandline option at some point


# Do the job and exit
main();
exit(0);


# Function to report usage information
sub usage {
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Collect stats about introns and other gap operation usage

Usage: perl $cmd [options] annot.gtf[.gz] [out-name]

If an out-name is given, the script expects SAM on stdin and will
create four output files:
 os/os.<out-name>.txt    (CIGAR operation size distribution)
 oc/oc.<out-name>.txt    (CIGAR operation counts)
 iac/iac.<out-name>.txt  (intron annotation comparison)
 is/is.<out-name>.txt    (intron size distribution)
The data in the latter two is stratified by read coverage.

If no out-name is given, the script writes the distribution of
annotated intron sizes to stdout.

Options:
-a file.gtf[.gz]   Second annotation file (use for simulated annotation)

EndOfUsage
    exit(1);
}


# Main function
sub main {
    
    my %args;
    getopts('a:', \%args);
    my $artif_annot_fn = $args{'a'};

    @ARGV == 1 or @ARGV == 2 or usage();
    my ($real_annot_fn, $out_name) = @ARGV;

    my $index1 = read_splices($real_annot_fn);
    my $index2 = read_splices($artif_annot_fn) if(defined $artif_annot_fn);

    if(defined $out_name) {
	my ($nr_aln, $op_sizes, $op_counts, $introns) = read_sam();
	make_op_size_report($nr_aln, $op_sizes, $out_name);
	make_op_count_report($nr_aln, $op_counts, $out_name);
	make_splice_reports($nr_aln, $introns, $out_name, $index1, $index2);
    }
    else {
	print_known_splice_sizes($index1);
    }
}


# Read SAM from STDIN and collect stats
sub read_sam {

    # Variables to record stats
    my $nr_aln = 0;
    my (@op_sizes, @op_counts, %introns);

    # Process one alignment per iteration
    while(my $line = <STDIN>) {

	# Parse SAM
	my (undef, $flags, $chr, $pos, undef, $cigar) = split "\t", $line;
	next if($flags & 0x4); # ignore non-alignments

	# Increment counter for number of alignments
	$nr_aln++;

	# Split up CIGAR string and initialize counters for looping over operations
	my @cig_fields = split(/([A-Z])/, $cigar);
	my @op_counts_for_read;
	my $i = 0;

	# Process one CIGAR operation per iteration
	while(@cig_fields) {
	    my $l = shift @cig_fields;
	    my $op = shift @cig_fields;
	    if($op eq 'M') {
		$pos += $l;
	    }
	    else {
		# Determine position of operation (0=first, 1=internal, 2=last)
		my $op_pos; 
		if($i == 0) { # First operation?
		    $op_pos = ($flags & 0x10) ? 2 : 0;
		}
		elsif(@cig_fields == 0) { # Last operation?
		    $op_pos = ($flags & 0x10) ? 0 : 2;
		}
		else {
		    $op_pos = 1;
		}
		# Check which operation it is
		if($op eq 'N') {
		    if($op_pos == 1) {
			$op = CIGAR_N;
			my $loc = "$chr:$pos:$l";
			$introns{$loc}++;
			$pos += $l;
		    }
		    else {
			warn "Warning: terminal N operation in CIGAR string: $cigar"; 
			undef $op;
		    }
		}
		elsif($op eq 'D') {
		    $op = $op_pos == 1 ? CIGAR_D : ($op_pos == 0 ? CIGAR_D0 : CIGAR_D1);
		    $pos += $l;
		}
		elsif($op eq 'I') {
		    $op = $op_pos == 1 ? CIGAR_I : ($op_pos == 0 ? CIGAR_I0 : CIGAR_I1);
		}
		elsif($op eq 'S') {
		    if($op_pos == 0) {
			$op = CIGAR_S0;
		    }
		    elsif($op_pos == 2) {
			$op = CIGAR_S1;
		    }
		    else {
			warn "Warning: internal S operation in CIGAR string: $cigar"; 
			undef $op;
		    }
		}
		else {
		    die "I don't know what to do with CIGAR operation $op";
		}
		# Count operation
		if(defined $op) {
		    $op_counts_for_read[$op]++;
		    $op_sizes[$op]{$l}++;
		}
	    }
	    $i++;
	}

	# Increment counts of operations per read
	for my $i (0..$#CIGAR_OPS) {
	    $op_counts[ $op_counts_for_read[$i] || 0 ][$i]++
	}
    }

    # Return stats
    return ($nr_aln, \@op_sizes, \@op_counts, \%introns);
}


# Make report on CIGAR operation sizes
sub make_op_size_report {
    my ($nr_aln, $counts, $out_name) = @_;

    # Open file and write headers
    mkdir("os") unless(-e "os");
    my $out_fn = "os/os.$out_name.txt";
    open REPORT, ">$out_fn" or die "Failed to open $out_fn for output.\n";
    print REPORT $nr_aln, "\n";
    print REPORT join("\t", "size", @CIGAR_OPS), "\n";

    # Write data
    my @sizes = map { keys %{$counts->[$_]} } (0..$#CIGAR_OPS);
    @sizes = sort_uniq(@sizes);
    foreach my $size (@sizes) {
	print REPORT join("\t", $size, map { $counts->[$_]{$size} || 0 } (0..$#CIGAR_OPS)), "\n";
    }

    # Close file
    close REPORT;
}


# Make report on CIGAR operation counts per read
sub make_op_count_report {
    my ($nr_aln, $counts, $out_name) = @_;

    # Open file and write headers
    mkdir("oc") unless(-e "oc");
    my $out_fn = "oc/oc.$out_name.txt";
    open REPORT, ">$out_fn" or die "Failed to open $out_fn for output.\n";
    print REPORT $nr_aln, "\n";
    print REPORT join("\t", "count", @CIGAR_OPS), "\n";

    # Write data
    for my $i (0..$#$counts) {
	print REPORT $i;
	for my $j (0..$#CIGAR_OPS) {
	    print REPORT "\t", ($counts->[$i][$j] || 0);
	}
	print REPORT "\n";
    }

    # Close file
    close REPORT;
}


# Make report on splices
sub make_splice_reports {
    my ($total_aln, $introns, $out_name, $index1, $index2) = @_;

    my (@class_counts, %size_counts);
    while(my ($intron_id, $nr_aln) = each %$introns) {
	# Get information about intron
	my ($chr, $pos, $len) = split ":", $intron_id;
	my $class1 = classify_intron($chr, $pos-1, $pos+$len, $index1);
	my $class2 = $index2 ? classify_intron($chr, $pos-1, $pos+$len, $index2) : 0;
	# Update per-read counts
	$class_counts[$class1][$class2][0] += $nr_aln;
	$size_counts{$len}[0] += $nr_aln;
	# Update per-splice counts
	$nr_aln = $MAX_ALN_COUNT if($nr_aln > $MAX_ALN_COUNT);
	$class_counts[$class1][$class2][$nr_aln]++;
	$size_counts{$len}[$nr_aln]++;
    }

    my $nr_indices = defined($index2) ? 2 : 1;
    mkdir("iac") unless(-e "ias");
    mkdir("is") unless(-e "is");
    write_splice_class_report("iac/iac.$out_name.txt", \@class_counts, $total_aln, $nr_indices);
    write_splice_size_report("is/is.$out_name.txt", $total_aln, \%size_counts);
}


## Write splice class report
## Output format:
##  cat1
##  cat2
##  number of reads
##  number of introns supported by 1 read
##  number of introns supported by 2 reads
##  etc

sub write_splice_class_report {
    my ($out_fn, $counts, $nr_aln, $nr_indices) = @_;
    
    open REPORT, ">$out_fn" or die "Failed to open $out_fn for output.\n";
    print REPORT $nr_aln, "\n"; # here we could prin the number of indexed splices too...

    print REPORT join("\t", (map { "class$_" } (1..$nr_indices)), "reads", (map { "n$_" } (1..$MAX_ALN_COUNT))), "\n"; 

    if($nr_indices == 1) {
	for my $i (0..$#INTRON_CLASS) {
	    print_splice_report_line(\*REPORT, $INTRON_CLASS[$i], $counts->[$i][0]);
	}
    }
    elsif($nr_indices == 2) {
	for my $i (0..$#INTRON_CLASS) {
	    for my $j (0..$#INTRON_CLASS) {
		print_splice_report_line(\*REPORT, $INTRON_CLASS[$i], $INTRON_CLASS[$j], $counts->[$i][$j]);
	    }
	}
    }
    else {
	die "Unsupported number of indices: $nr_indices";
    }

    close REPORT;
}


##  Write splice sizes report
##  Output format:
##  size
##  number of reads
##  number of introns supported by 1 read
##  number of introns supported by 2 reads
##  etc

sub write_splice_size_report {
    my ($out_fn, $nr_aln, $counts) = @_;
    open REPORT, ">$out_fn" or die "Failed to open $out_fn for output.\n";
    print REPORT $nr_aln, "\n";
    print REPORT join("\t", "size", "reads", map { "n$_" } (1..$MAX_ALN_COUNT)), "\n"; 
    foreach my $size (sort { $a <=> $b } keys %$counts) {
	print_splice_report_line(\*REPORT, $size, $counts->{$size});
    }
    close REPORT;
}


sub print_splice_report_line {
    my ($out, @data) = @_;

    # Counts is the last argument, the rest are meta data
    my $counts = pop @data;

    # Print meta data
    print $out join("\t", @data);

    # Print counts
    for my $i (0..$MAX_ALN_COUNT) {
	print $out "\t", ($counts->[$i] || 0);
    }

    # End line
    print $out "\n";
}


sub print_known_splice_sizes {
    my $index = shift;
    my $splices = $index->[0];
    my %size_counts;
    print "size\tcount\n";
    foreach my $chr (keys %$splices) {
	foreach my $left (keys %{$splices->{$chr}}) {
	    foreach my $right (keys %{$splices->{$chr}{$left}}) {
		$size_counts{$right - $left - 1}++;
	    }
	}
    }
    foreach my $size (sort { $a <=> $b } keys %size_counts) {
	print $size, "\t", $size_counts{$size}, "\n";
    }    
}


sub classify_intron_DBG {
    my ($chr, $left, $right) = @_;   
    my $class = classify_intron(@_);
    print "$chr:$left-$right $class\n";
    return $class;
}


sub classify_intron {
    my ($chr, $left, $right, $index) = @_;

    my ($splices, $left2gene, $right2gene) = @$index;

    if($left2gene->{$chr}{$left}) {  # Check one end first to avoid autovivification
	$splices->{$chr}{$left} or die "index error"; # Paranoia
	if($splices->{$chr}{$left}{$right}) { # This is the most common case by far, so catch it first
	    return KNOWN_INTRON;
	}
	elsif($right2gene->{$chr}{$right}) {
	    if(same_gene($left2gene->{$chr}{$left}, $right2gene->{$chr}{$right})) {
		return SAME_GENE;
	    }
	    else {
		return DIFF_GENES;
	    }
	}
	else {
	    return ONE_NOVEL;
	}
    }
    elsif($right2gene->{$chr}{$right}) {
	return ONE_NOVEL;
    }
    else {
	return BOTH_NOVEL;
    }
}


sub same_gene {
    my ($genes1, $genes2) = @_;
    foreach my $g1 (ref($genes1) ? @$genes1 : $genes1) {
	foreach my $g2 (ref($genes2) ? @$genes2 : $genes2) {
	    return 1 if ($g1 eq $g2);
	}
    }
    return 0;
}


# This function reads a GTF file to create annotation indices:
# %splices{chr}{left}{right}
# %left2gene{chr}{left}
# %right2gene{chr}{right}

sub read_splices {
    my $fn = shift;

    # Read file
    my (%tx_exons, %tx_chr, %tx_gene);
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
	# Get gene and tx id
	my ($gene_id) = $attrib =~ /gene_id "([\w\.]+)"/;
	my ($tx_id) = $attrib =~ /transcript_id "([\w\.]+)"/;
	(defined($gene_id) and defined($tx_id)) or die "Failed to parse attribute string: $attrib";
	# Store gene id in index
	if(exists $tx_gene{$tx_id}) {
	    $tx_gene{$tx_id} eq $gene_id or die "Tx $tx_id belongs to multiple genes";     
	}
	else {
	    $tx_gene{$tx_id} = $gene_id;
	}
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
    my (%splices, %left2gene, %right2gene);
    foreach my $tx_id (keys %tx_exons) {
	my $chr = $tx_chr{$tx_id};
	my $gene_id = $tx_gene{$tx_id};
	my @sorted_coords = sort { $a <=> $b } @{$tx_exons{$tx_id}};
	shift @sorted_coords; # Ignore start of first exon
	pop @sorted_coords;   # Ignore end of last exon
	while(@sorted_coords) {
	    my $left = shift @sorted_coords;
	    my $right = shift @sorted_coords;
	    $splices{$chr}{$left}{$right} = 1;
	    $left2gene{$chr}{$left}{$gene_id} = 1;
	    $right2gene{$chr}{$right}{$gene_id} = 1;
	}
    }

    # Change gene id hashes to scalars or array refs
    foreach my $idx (values(%left2gene), values(%right2gene)) {
	foreach my $coord (keys %$idx) {
	    my @gene_ids = keys %{$idx->{$coord}};
	    $idx->{$coord} = @gene_ids == 1 ? $gene_ids[0] : \@gene_ids;
	}
    }

    # Return index
    return [ \%splices, \%left2gene, \%right2gene ];
}


# Write splice index to file. For debug purposes.

sub write_splices {
    my $index = shift;

    my ($splices, $left2gene, $right2gene) = @$index;
    
    open OUT, ">splices.bed" or die;
    foreach my $chr (keys %$splices) {
	foreach my $left (keys %{$splices->{$chr}}) {
	    foreach my $right (keys %{$splices->{$chr}{$left}}) {
		print OUT join("\t", "chr$chr", $left-1, $right), "\n";
	    }
	}
    }
    close OUT;

    open OUT, ">left.bed" or die;
    foreach my $chr (keys %$left2gene) {
	foreach my $coord (keys %{$left2gene->{$chr}}) {
	    my $g = $left2gene->{$chr}{$coord};
	    $g = join("/",@$g) if ref($g);
	    print OUT join("\t", "chr$chr", $coord-1, $coord, $g), "\n";
	}
    }
    close OUT;

    open OUT, ">right.bed" or die;
    foreach my $chr (keys %$right2gene) {
	foreach my $coord (keys %{$right2gene->{$chr}}) {
	    my $g = $right2gene->{$chr}{$coord};
	    $g = join("/",@$g) if ref($g);
	    print OUT join("\t", "chr$chr", $coord-1, $coord, $g), "\n";
	}
    }
    close OUT;
}


# sort | uniq

sub sort_uniq {
    my %hash = map { $_ =>  undef } @_;
    return(sort {$a <=> $b} keys %hash)
}


