#!/bin/env perl
# Count the number of alignments in introns and intronic repeats
use warnings;
use strict;
use Getopt::Std;
use Genoman::Misc qw(strand_sign_to_numeric strand_numeric_to_sign);
use Genoman::BlockSet;

# Do the job and exit
main();
exit(0);


# Function to report usage information
sub usage {
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Count the number of alignments in introns and intronic repeats

Usage: perl $cmd genes.gtf[.gz] rmsk.txt[.gz]

Arguments:
1. Gene annotation file from Ensembl in GTF format
2. Repeat annotation file from the UCSC Genome Browser database

The script expects SAM on stdin.

EndOfUsage
    exit(1);
}


# Main function
sub main {
    
    @ARGV == 2 or usage();
    my ($ref_fn, $rmsk_fn) = @ARGV;

    # Read reference transcripts
    print STDERR "Reading reference annotation...\n";
    my $ref_tx = read_gtf($ref_fn);

    # Read repeats
    print STDERR "Reading repeats...\n";
    my $repeats = read_rmsk($rmsk_fn);

    # Index introns
    print STDERR "Indexing introns...\n";
    my $introns = index_introns($ref_tx);

    # Process alignments
    print STDERR "Processing SAM data from STDIN...\n";
    process_sam($introns, $repeats);
}


sub process_sam {
    my ($introns, $repeats) = @_;

    # This is where we store the counts.
    my $intronic_non_repeat = 0;
    my $intronic_repeat = 0;
    my $total = 0;

    while(my $line = <STDIN>) {

	# Split line
	chomp $line;
	my @sam = split "\t", $line;

	# Ignore non-alignments
	next if($sam[1] & 0x4);

	# Count alignment
	$total++;

	# Get coords
	my $chr = $sam[2];
	my $start = $sam[3];
	my $span = calc_aln_span($sam[5]);
	my $end = $start + $span - 1;

	# Check if fully intronic
	next if(calc_overlap($introns, $chr, $start, $end) < $span);

	# Check if overlaps repeats
	if(calc_overlap($repeats, $chr, $start, $end)) {
	    $intronic_repeat++;
	} else {
	    $intronic_non_repeat++;
	}

	# Progress indicator
	$total++;
	if($total % 100000 == 0) { print STDERR " processed $total alignments\n"; }
    }

    print STDERR " processed $total alignments (all)\n";

    print $intronic_non_repeat, "\t", $intronic_repeat, "\t", $total, "\n";
}


sub index_introns {   
    my $tx_list = shift;

    # Get all exons and gene bounds for each chr
    my %exons;
    my %gene_bounds;
    foreach my $tx (@$tx_list) {
	my ($chr, $exons) = @$tx[2,4];
	push @{$exons{$chr}}, @$exons;
	push @{$gene_bounds{$chr}}, ( $exons->[0], $exons->[-1] );
    }

    # Determine introns for each chr
    my %introns;
    foreach my $chr (keys %gene_bounds) {
	my $genic = Genoman::BlockSet->new_from_unsorted($gene_bounds{$chr});
	my $non_exonic = Genoman::BlockSet->new_from_unsorted($exons{$chr})->inverse();
	$introns{$chr} = $genic->intersection($non_exonic);
    }

    return \%introns;
}


sub read_gtf {
    my ($fn, $valid_chr) = @_;

    # Hashes to store indices
    my %tx_info;

    # Open file
    if($fn =~ /.gz$/i) {
	$fn = "zcat $fn |";
    }
    open(GTF, $fn) or die "Failed to open $fn for reading\n";

    # Read GTF file
    while(my $line = <GTF>) {
	# Parse GTF line
	chomp $line;
	my ($chr, $tx_type, $ft_type, $start, $end, undef, $strand, undef, $attrib) = split "\t", $line;
	next unless($ft_type eq 'exon'); # Only consider exons
	$strand = strand_sign_to_numeric($strand);
        # Parse attribute string
	my ($tx_id) = $attrib =~ /transcript_id "([^"]+)"/;
	unless(defined($tx_id)) { die "Failed to parse attribute string: $attrib" }
	# Update transcript info
	if(exists $tx_info{$tx_id}) {
	    $tx_info{$tx_id}[1] eq $tx_type or die "ERROR: conflicting transcript types ($tx_info{$tx_id}[0], $tx_type) for transcript $tx_id in $fn";
	    $tx_info{$tx_id}[2] eq $chr or die "ERROR: conflicting chromosome for transcript $tx_id in $fn";
	    $tx_info{$tx_id}[3] == $strand or die "ERROR: conflicting strand for transcript $tx_id in $fn\n";
	    push @{$tx_info{$tx_id}[4]}, ($start, $end);
	}
	else {
	    $tx_info{$tx_id} = [$tx_id, $tx_type, $chr, $strand, [$start, $end]];
	}
    }
    close(GTF) or die "Failed to close file $fn";

    # Finalize index by collapsing exons
    foreach my $tx (values %tx_info) {
	$tx->[4] = Genoman::BlockSet->new_from_unsorted($tx->[4]);
    }

    # Return indices
    return [ @tx_info{ sort keys %tx_info } ];
}


sub read_rmsk {
    my $fn = shift;

    # Open file
    if($fn =~ /.gz$/i) {
	$fn = "zcat $fn |";
    }
    open(RMSK, $fn) or die "Failed to open $fn for reading\n";

    # Read and index repeats
    my %repeats;
    while(my $line = <RMSK>) {
	my (undef, undef, undef, undef, undef, $chr, $start, $end) = split "\t", $line;
	$start++;
	push @{$repeats{$chr}}, ($start, $end);
    }

    close(RMSK) or die "Failed to close file $fn";

    # Finalize index by collapsing exons
    foreach my $chr (keys %repeats) {
	$repeats{$chr} = Genoman::BlockSet->new_from_unsorted($repeats{$chr});
    }   

    return \%repeats;
}


sub calc_aln_span {
    my ($cigar) = @_;
    my $span = 0;
    my @cig_fields = split(/([A-Z])/, $cigar);
    while(@cig_fields) {
	my $l = shift @cig_fields;
	my $op = shift @cig_fields;
	if($op eq 'M' or $op eq 'D' or $op eq 'N') {
	    $span += $l;
	}
	elsif($op ne 'I' and $op ne 'S') {
	    die "I don't know what to do with CIGAR operation $op";
	}
    }
    return $span;
}


sub calc_overlap {
    my ($index, $chr, $start, $end) = @_;
    $index = $index->{$chr} or return 0;
    return $index->slice($start, $end)->sum();
}


