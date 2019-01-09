#!/usr/bin/perl 
use warnings;
use strict;


unless (exists $ARGV[1]) {print "\n[file] [seq prefix]\n";exit;}

open FASTA, $ARGV[0];
my $count = 1;

{
local $/ = '>'; 
<FASTA>;                                             # throw away the first line 'cos will only contain ">"

while (<FASTA>) 
	{	
    	chomp $_;
    	my ($seq_id, @sequence) = split "\n";            # split the fasta input into Id and sequence
    	my $fasta_sequence = join '',@sequence;          # reassembles the sequence
	print "\>$ARGV[1]\_$count\n$fasta_sequence\n";
	$count++;
	}
}
close FASTA;