#!/usr/bin/perl 
use warnings;
use strict;
use Getopt::Long;

my $usage = "============================================================================================
USAGE:
-f		fasta file
-s		Salmon TPM file from Transrate e.g xxx.quant.sf
-t		Minimum TPM e.g. 1
-l		Minimum transcript length in bp
============================================================================================
" ;

my ($fasta, $salmon, $min_tpm, $min_length);

GetOptions(
    'f|fasta:s'     => \$fasta,
    's|salmon:s'     => \$salmon,
	't|tpm:s'     => \$min_tpm,
	'l|length:s'     => \$min_length,
	);

if( ! defined $fasta) {
print "$usage\nWARNING: Cannot proceed without input file fasta file\n\n"; exit;
}
if( ! defined $salmon) {
print "$usage\nWARNING: Cannot proceed without input salmon quantification file\n\n"; exit;
}
if( ! defined $min_tpm) {
print "$usage\nWARNING: Cannot proceed without minimum tpm value\n\n"; exit;
}
if( ! defined $min_length) {
print "$usage\nWARNING: Cannot proceed without minimum transcript length\n\n"; exit;
}


open TPM, $salmon;
my %tpm_lookup;
while (<TPM>)	
	{	
	chomp $_;
	unless ($_ =~ /^\#/)
		{
		my @data = split '\t', $_;
		if ($data[2] >= $min_tpm && $data[1] >= $min_length)
			{
			$tpm_lookup{$data[0]} = $data[2];
			}
		}
	}



open FASTA, $fasta;
open OUT, ">$fasta\.tpm$min_tpm\.minlength$min_length\.fa";
{
local $/ = '>'; 
<FASTA>;                                             # throw away the first line 'cos will only contain ">"
while (<FASTA>) 
	{	
    	chomp $_;
    	my ($seq_id, @sequence) = split "\n";            # split the fasta input into Id and sequence
	my $fasta_sequence = join '',@sequence;          # reassembles the sequence
    	if (exists $tpm_lookup{$seq_id})
		{
		print OUT "\>$seq_id\n$fasta_sequence\n";
		}
	}
}
close FASTA;
close OUT;