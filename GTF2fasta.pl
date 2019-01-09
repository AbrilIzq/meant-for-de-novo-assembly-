#!/usr/bin/perl 
use warnings;
use strict;


use Getopt::Long;
# Richard Emes University of Nottingham 2016

my $usage = "
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
R D Emes University of Nottingham 2016

pull genes from a reference genome GFF file and a fasta file of resequenced DNA.
copes with soft clipping meaning reads map over ends of transcript/genome [in this case will return to end of reference fasta]

USAGE:
-f	fasta file of genomic DNA or transcriptome
-g	gtf file where columns 4/5 are start end points
-o	output
-n	name of entity to parse e.g \"transcript\" \"exon\" etc must match column 3 of GTF
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
" ;


my ($file, $gtf, $name, $out);

GetOptions(
    'f|fasta:s'     => \$file,
    'g|gtf:s'   => \$gtf,
	'n|name:s'   => \$name,
	'o|output:s'   => \$out,	
                 );


if( ! defined $file) {
print "$usage\nWARNING: Cannot proceed without fasta file to process\n\n"; exit;
}
if( ! defined $gtf) {
print "$usage\nWARNING: Cannot proceed without genomic DNA file\n\n"; exit;
}
if( ! defined $name) {
print "$usage\nWARNING: Cannot proceed without name of entity to parse\n\n"; exit;
}
if( ! defined $out) {
print "$usage\nWARNING: Cannot proceed without output file\n\n"; exit;
}
####################################################

my %lookup; # hash of arrays key is fasta sequence name array contains details in single string
## read GTF 
open FILE, "<$gtf";
while (<FILE>)
{
chomp $_;
unless ($_ =~ /^\#/)
	{
	my $line = $_;
	my @data = split '\t', $line;
	
	if ($data[2] eq $name)
		{
		my $name = $data[0];
		my $start = $data[3];
		my $end = $data[4];
		my $details = $data[8];

		$details =~ s/\"//g;
		$details =~ s/\;//g;
		$details =~ s/\'//g;
		if ($details =~ /(.*?)\s$/) {$details = $1};
		my $hash_details = $start."@".$end."@".$details;
		push(@{$lookup{$name}}, $hash_details);
		}
	}
}
close FILE;



open OUT, ">$out";
# read in fasta file 
my $fasta_sequence;

{
	open FASTA, "<$file";
	{
	local $/ = '>'; 
	<FASTA>;									# throw away the first line 'cos will only contain ">"

	while (<FASTA>) 
		{	
		chomp $_;
		my ($seq_id, @sequence) = split "\n";	# split the fasta input into Id and sequence
		$fasta_sequence = join '',@sequence;	# reassembles the sequence
		my $seq_length = length $fasta_sequence;
		
		if (exists $lookup{$seq_id})
			{
			foreach (@{$lookup{$seq_id}})
				{
				chomp $_;
				my ($start,$end,$details) = split '\@', $_;
				
				if ($end > $seq_length){$end = $seq_length;} # in case of soft clipping of reads map over the end of contig
				if ($start < 1){$start = 1} # in case of soft clipping of reads map over the end of contig
				$details =~ s/gene_id //g;
				$details =~ s/ transcript_id /_/g;
				my $length = ($end-$start)+1; 	# to account for substr start and end
				my $sub_start = $start-1 ;		# to account for substr start and end
				my $seq = substr($fasta_sequence, $sub_start, $length); #$seq, start, length of substring)
				print OUT "\>$details\n$seq\n";
				}
			}
		}
	close FASTA;
}
}
close OUT;
