#!/bin/env perl

# Filter simulated transcript GTF file to exclude silent transcripts

use warnings;
use strict;
use Getopt::Std;

# Parse command line args

my %args;
getopts("m", \%args);
my $mature_only = $args{m};

# Report usage information?
unless(@ARGV == 1) {

    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Filter simulated transcript GTF file to exclude silent transcripts

Usage: perl $0 [-m] expr.txt <in.gtf >out.gtf

The expression file given on the command line should be a
tab-delimited (compact) quantification file.

By default, all transcripts with exonic or intronic read counts are
included. Use option -m to only include transcripts expressed in
mature form (i.e. with exonic counts).

EndOfUsage
    exit(1);
}

my $fpkm_fn = $ARGV[0];

# Read FPKM for all transcripts
open FPKM, $fpkm_fn or die "Failed to open $fpkm_fn";
my @read_counts;
my $line = <FPKM>;
die "header error in $fpkm_fn" unless($line eq "id\tsize\texons\treads.exonic\treads.intronic\tfpkm\n");
while(my $line = <FPKM>) {
    chomp $line;
    my @fields = split "\t", $line;
    $read_counts[ $fields[0] ] = $mature_only ? $fields[3] : $fields[3] + $fields[4];
}

# Filter GTF
while(my $line = <STDIN>) {
    my @fields = split "\t", $line;
    my ($tx_id) = $fields[8] =~ /transcript_id \"(\d+)\";/;
    print $line if($read_counts[$tx_id] > 0);
}
