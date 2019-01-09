#!/usr/bin/perl 
use warnings;
use strict;
use Getopt::Long;

my $usage = "============================================================================================
USAGE:
-f	fasta file
-c	cd-est-hit cluster xxx.est.clstr
-m	Minimum number of sequences in cluster to retain
-o	Out file name
-g	should minimum cluster size be from different input files (Y/N) Assumes that sequences are named xxx_yyy where yyy is a number
============================================================================================
" ;

my ($fasta, $cluster, $min_seq, $output, $groups);

GetOptions(
    'f|fasta:s'     => \$fasta,
    'c|cluster:s'     => \$cluster,
	'm|min:s'     => \$min_seq,
	'o|output:s'     => \$output,
	'g|group:s'     => \$groups,
	);

if( ! defined $fasta) {
print "$usage\nWARNING: Cannot proceed without input file fasta file\n\n"; exit;
}
if( ! defined $cluster) {
print "$usage\nWARNING: Cannot proceed without cd-est-hit cluster file\n\n"; exit;
}
if( ! defined $min_seq) {
print "$usage\nWARNING: Cannot proceed without minimum number in cluster\n\n"; exit;
}
if( ! defined $output) {
print "$usage\nWARNING: Cannot proceed without output file name\n\n"; exit;
}
if( ! defined $groups) {
print "$usage\nWARNING: Cannot proceed without group (-g) should minimum cluster size be from different input files (Y/N)\n\n"; exit;
}
$groups = uc($groups);
unless ($groups =~ /[YN]/){print "$usage\nWARNING: Cannot proceed without group (-g) should minimum cluster size be from different input files (Y/N)\n\n"; exit;}



my $array_size_needed = $min_seq -1 ; # array size needed will be 1 less than minimum due to array index starting at 0

my %cluster_lookup;
	
if ($groups eq "N")
	{
	open CLUSTER, $cluster;
	{
	local $/ = '>Cluster'; 
	<CLUSTER>;                                             # throw away the first line 'cos will only contain ">"
	while (<CLUSTER>)	
		{	
		chomp $_;
		my ($cluster_id, @sequences) = split "\n";            # split the fasta input into Id and sequence

		if (exists $sequences[$array_size_needed])
			{
			foreach(@sequences)
				{
				chomp $_;
				if ($_ =~ /.*?\>(.*?)\.\.\.\s\*/)
					{
					$cluster_lookup{$1} = 1;
					}
				}
			}
		}
	}
	close CLUSTER;
}



if ($groups eq "Y")
	{
	open CLUSTER, $cluster;
	{
	local $/ = '>Cluster'; 
	<CLUSTER>;                                             # throw away the first line 'cos will only contain ">"
	while (<CLUSTER>)	
		{	
		chomp $_;
		my ($cluster_id, @sequences) = split "\n";            # split the fasta input into Id and sequence

		if (exists $sequences[$array_size_needed])
			{
			my $longest = "NA";
			my @groups = ();
			foreach(@sequences)
				{
				chomp $_;
				if ($_ =~ /.*?\>(.*?)\_\d+\.\.\./)
					{
					push @groups, $1;
					}
				if ($_ =~ /.*?\>(.*?)\.\.\.\s\*/)
					{
					$longest = $1;
					}
				}
				my @uniq_groups = uniq_array(@groups);
				
				if (exists $uniq_groups[$array_size_needed])
					{
					$cluster_lookup{$longest} = 1;
					}
			}
		}
	}
	close CLUSTER;
}




open FASTA, $fasta;
open OUT, ">$output";
{
local $/ = '>'; 
<FASTA>;                                             # throw away the first line 'cos will only contain ">"
while (<FASTA>) 
	{	
    	chomp $_;
    	my ($seq_id, @sequence) = split "\n";            # split the fasta input into Id and sequence
		my $fasta_sequence = join '',@sequence;          # reassembles the sequence
    	if (exists $cluster_lookup{$seq_id})
		{
		print OUT "\>$seq_id\n$fasta_sequence\n";
		}
	}
}
close FASTA;
close OUT;



sub uniq_array{
##### make a unique list from the @genes array
my @in = @_;
my %seen = ();
my @uniq = ();
foreach (@in)
{chomp $_;
unless ($seen{$_}){
$seen{$_} = 1;
push (@uniq, $_);
}
}

return @uniq;
}