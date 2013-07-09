#!/bin/env perl

# Calculate the distribution of number of mismatches per alignment stratified by mean read quality

use warnings;
use strict;

sub usage {
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Calculate the distribution of number of mismatches per alignment
stratified by mean read quality

Usage: $cmd <read_folder> <input.bam> <phred_coding>

Arguments:
1. Folder containing the read files (and preferably no other files)
   compressed with gzip. This script expects paired-end data, and read
   number should be indicated by _1 or _2 in file names.
2. BAM file to analyze
3. A number e.g. 33 or 64 indicating the quality score coding.

IMPORTANT: the read identifiers in the bam file and the fastq files
must be the same! Sometimes if the read identifier in the fastq file
contains spaces only the first part of the identifier will be in
the bam file.

EndOfUsage
    exit;
}

usage() unless(@ARGV == 3);

# First we compute the average quality for each read
my $fa_folder = $ARGV[0];
my $bam_file = $ARGV[1];
my $phred = $ARGV[2];

my %read_seqs_1 = ();
my %read_seqs_2 = ();

my %scores = ();

opendir(DIR,$fa_folder);
my @files = readdir(DIR);
closedir(DIR);

# Read fastq files, calculate average quality and store in hash
foreach my $file (@files) {
  next if $file =~ /^\./;
  next unless $file =~ /.(txt|fq|fastq)+(.gz)?$/;
  my ($mate) = $file =~ /.*_([12]).*.gz/;
  my $command = "cat";
  print STDERR "$file\t$mate\n";
  if ($file =~ /.gz/) { $command = "zcat"; }
  open FA, "$command $fa_folder$file |" or die "Could not open $fa_folder$file with $command\n";
  while (my $fa_line = <FA>) {
    chomp $fa_line;
    my $id = substr $fa_line, 1;
    $id = substr $id, 0, -2;
    $fa_line = <FA>;
    $fa_line = <FA>;
    my $qual_seq = <FA>;
    chomp $qual_seq;
    my $score = 0;
    map { $score += ord($_)-$phred } split("",$qual_seq);
    my $avg_score = $score / length($qual_seq);
    if ($mate == 1) { $read_seqs_1{$id} = int($avg_score); }
    else { $read_seqs_2{$id} = int($avg_score); }
	if (!exists($scores{int($avg_score)})) {
      $scores{int($avg_score)} = 1;
    }
    else {
      $scores{int($avg_score)}++;
    }
  }
  close FA;
}

# Initialize array of hashes for each possible read quality
my @mm_count = ();
for (my $i = 0; $i <= 40; $i++) {
  %{$mm_count[$i]} = ();
}

# Read bam file and count number of mismatches for each alignment.
# Alignment against N in reference or alignment does not count as mismatch.
# Make sure not to count double when N is aligned against N.
open BAM, "samtools view $bam_file | " or die "Could not open $bam_file with samtools view\n";
while (my $line = <BAM>) {
  my @values = split("\t",$line);
  if ($values[1] & 0x4) { next; }
  if ($values[1] & 0x100) { next; }
  if ($values[2] eq 'chrM') { next; }
  my $mate = 1;
  if ($values[1] & 0x80) { $mate = 2; }
  my @cig_fields = split(/([A-Z])/, $values[5]);
  my $avg_score = 0;
  if (($mate == 1 && !exists($read_seqs_1{$values[0]})) || ($mate == 2 && !exists($read_seqs_2{$values[0]}))) {
    print STDERR "$line\n$values[0]\n";
    last;
  }
  if ($mate == 1) { $avg_score = $read_seqs_1{$values[0]}; }
  else { $avg_score = $read_seqs_2{$values[0]}; }
  my $indels = 0;
  my $pos = 0;
  my $numN = 0;
  my $numN_ref = 0;
  my $s = 0;
  my %numN_pos = ();
  my %numN_pos_ref = ();
  while(@cig_fields) {
    my $l = shift @cig_fields;
    my $op = shift @cig_fields;
    if ($op eq 'I' || $op eq 'D') {
      $indels += $l;
    }
    if ($op eq 'S') {
      $s = $l;
    }
    if ($op eq 'M') {
      my $sub_seq = substr $values[9], $pos, $l;
      my $n = $sub_seq =~ tr/N//;
      $numN += $n;
      if ($n > 0) {
 	my $p = index($sub_seq,'N');
 	while ($p != -1) {
 	  $numN_pos{$pos + $p - $s} = 1;
 	  $p = index($sub_seq, 'N', $p + 1);
 	}
      }
    }
    if ($op eq 'I' || $op eq 'M' || $op eq 'S') {
      $pos += $l;
    }
  }
  my ($md) = $line =~ m/MD:Z:(\S+)/g;
  $md =~ s/\^\D+/:/g;
  $numN_ref = $md =~ tr/N//;
  $pos = 0;
  if ($numN_ref > 0) {
    my @md_fields = split(/([a-zA-Z:])/,$md);
    while (@md_fields) {
      my $l = shift @md_fields;
      $pos += $l;
      if ($#md_fields == -1) { next; }
      my $mm = shift @md_fields;
      if ($mm eq 'N') { $numN_pos_ref{$pos} = 1; }
      if ($mm ne ':') { $pos++; }
    }
  }
  $numN = 0;
  # N-N is a mismatch for samtools make sure it is only counted exactly once and not twice
  for my $key (keys %numN_pos) {
    if (!exists($numN_pos_ref{$key})) {
      $numN++;
    }
  }
  my ($nm) = $line =~ m/NM:i:(\d+)/g;
  $nm = $nm - $indels - $numN - $numN_ref;
  if ($nm < 0) {
    print STDERR $line;
    print STDERR "nm=$nm, indels=$indels, numN=$numN, numN_ref=$numN_ref\n";
    $nm = 0;
  }
  if (exists $mm_count[$avg_score]{$nm}) {
    $mm_count[$avg_score]{$nm}++;
  }
  else {
    $mm_count[$avg_score]{$nm} = 1;
  }
}

# Report number of mismatches for each mean read quality
for (my $i = 0; $i <= 40; $i++) {
  foreach my $key (sort keys %{$mm_count[$i]}) {
    print "$i\t$key\t$mm_count[$i]{$key}\n";
  }
}
