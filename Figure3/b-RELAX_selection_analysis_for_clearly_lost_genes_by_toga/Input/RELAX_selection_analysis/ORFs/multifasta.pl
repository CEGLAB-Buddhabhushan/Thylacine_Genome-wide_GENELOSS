#!/usr/bin/perl
use strict;
use warnings;

# Get the multi-FASTA file from the command-line argument
my $multifastafile = $ARGV[0];
open(my $MULTIFASTA, "<", $multifastafile) || die "File not found: '$multifastafile'.\n";

my %sequence;
my ($header, $seq) = ("", "");

# Read the multi-FASTA file and store each sequence with its header
while (<$MULTIFASTA>) {
    chomp;
    if (/^>(.*)/) {  # If the line is a header
        if ($seq) { 
            $sequence{$header} = $seq;
        }
        $header = $1;
        $seq = '';
    } else {
        $seq .= $_;
    }
}

# Store the last sequence if the file doesn't end with a header line
$sequence{$header} = $seq if $seq;

# Process each sequence
foreach my $name (keys %sequence) {
    my $newseq = $sequence{$name};
    
    # Check if the sequence starts with "ATG" and follows the triplet pattern
    if ($newseq =~ /^ATG(?:[ATGCN]{3})*$/) {
        print ">$name\n$newseq\n";
    }
}

close $MULTIFASTA;
