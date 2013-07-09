#!/bin/env perl

# From a GTF file, extract the type for each gene
# Expects GTF on STDIN.
# Output: a table associating gene IDs with gene types

use warnings;
use strict;

my %gene2type;

while(my $line = <STDIN>) {
    my @fields = split "\t", $line;
    my $type = $fields[1];
    my ($gene_id) = $fields[8] =~ /gene_id \"(\w+)\"/;
    if(exists $gene2type{$gene_id}) {
	die unless($gene2type{$gene_id} eq $type);
    }
    else {
	$gene2type{$gene_id} = $type;
    }
}

foreach my $gene_id (sort keys %gene2type) {
    print $gene_id, "\t", $gene2type{$gene_id}, "\n";
}
