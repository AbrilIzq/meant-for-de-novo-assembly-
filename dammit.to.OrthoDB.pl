#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $usage = "
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Parse OrthoDB
USAGE : 

-s species to match with Homo sapiens etc CHECK SPELLINGS quote if space e.g. \"Homo sapiens\"
-d dammit output final gff3 file
-o orthoDB file e.g ODB8_EukOGs_genes_ALL_levels.txt
-t orthoDB taxonomy level Vertebrata, Metazoa etc 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
" ;


my ($species, $dammit, $ortho, $taxa);

GetOptions(
        'd|dammit:s'     => \$dammit,
        's|species:s'     => \$species,
        'o|ortho:s'     => \$ortho,
        't|taxa:s'     => \$taxa
);


if( ! defined $dammit) {
print "\n\nWARNING: Cannot proceed without dammit file to process\n$usage\n\n"; exit;
}
if( ! defined $species) {
print "\n\nWARNING: Cannot proceed without species to output\n$usage\n\n"; exit;
}
if( ! defined $ortho) {
print "\n\nWARNING: Cannot proceed without orthoDB file\n$usage\n\n"; exit;
}
if( ! defined $taxa) {
print "\n\nWARNING: Cannot proceed without orthoDB taxa\n$usage\n\n"; exit;
}
############################################################################
############################################################################



############################################################################
## Process the dammit file for LAST hits
############################################################################
print "processing Dammit file\n";
my %dammit_lookup; 
my %last_lookup; #list of LAST IDS

open DAMMIT, $dammit;
while (<DAMMIT>)
	{
	chomp $_;
	my @data = split '\t', $_;
	my $transcript = $data[0];
	if ($data[1] eq "LAST")
		{
		my $details = $data[8];
		
			if ($details =~ /ID\=.*?Target\=(.*?)[\s\t].*/)
				{
				$dammit_lookup{$transcript} = $1;
				$last_lookup{$1} = 1;
				}
		}
	}
close DAMMIT;


print "processing OrthoDB file\n";

my %ortholookup; ### need to add array for hash of all ORthoDB hits for target species
my %description_lookup;
my %matched_eukogs;



open ORTHO, $ortho;
while (<ORTHO>)
	{
	chomp $_;
	my @data = split'\t', $_;
	
	my $ortho_description = "NA";
	if (exists $data[5]) # data[5] = description
		{
		$ortho_description = $data[5];
		}
	
	if ($data[0] =~ /.*?\:$taxa/)
		{
		if ($data[4] eq $species) 
			{
			push @{$ortholookup{$data[1]}}, $data[2]; # data[2] = protein id hash of arrays of all matching orthologs in species of interest
			$description_lookup{$data[1]} = $ortho_description; # data[1] = KOG id
			}
		}
	if (exists $last_lookup{$data[2]})
		{
		$matched_eukogs{$data[2]} = $data[1];
		} 
	}
close ORTHO;


my $species_out = $species;
$species_out =~ s/ /_/g;


## match transcripts per matching orthoKOG
my %transcript_count_per_KOG; # hash of arrays of transcripts 
my %unique_KOGS;

foreach my $key (sort(keys %dammit_lookup))
	{
		my $KOG = "NA";
		my $id = $dammit_lookup{$key}; # all new transcripts
		
		if (exists $matched_eukogs{$id})
			{
			$KOG = $matched_eukogs{$id};
			$unique_KOGS{$KOG} = 1;
			
			if (exists $transcript_count_per_KOG{$KOG})
				{
				my $count = $transcript_count_per_KOG{$KOG};
				$count++;
				$transcript_count_per_KOG{$KOG} = $count;
				}
			else 
				{
				$transcript_count_per_KOG{$KOG} = 1;
				}
			}
	}

	

open OUT, ">$dammit\_$taxa\_$species_out\.orthodb.out";
print OUT "Transcript\tOrthologs\tEukOG\t$species_out\tDescription\tEukOG $species Count\tEukOG Transcriptome Count\n";

foreach my $key (sort(keys %dammit_lookup))
	{
		my $id = $dammit_lookup{$key};

		my $KOG = "NA";
		my $matched_spp = "NA";
		my $description = "NA";
		my $species_KOG_count = "NA";
		my $transcripts_per_KOG_count = "NA";	
		
		if (exists $matched_eukogs{$id})
			{
			$KOG = $matched_eukogs{$id};
			}
		if (exists $ortholookup{$KOG})
			{
			my @orthos = @{${ortholookup{$KOG}}};
			$matched_spp = join ',' , @orthos;
			$species_KOG_count = scalar @orthos;
			}
		if (exists $description_lookup{$KOG})
			{
			$description = $description_lookup{$KOG};
			}
		if (exists $transcript_count_per_KOG{$KOG})
			{
			$transcripts_per_KOG_count = $transcript_count_per_KOG{$KOG};
			}
	print OUT "$key\t$id\t$KOG\t$matched_spp\t$description\t$species_KOG_count\t$transcripts_per_KOG_count\n";
	}
close OUT;


print "preparing cluster file\n";

open RERUN, "$dammit\_$taxa\_$species_out\.orthodb.out";
open OUT2, ">$dammit\_$taxa\_$species_out\.orthodb.KOGs";

my %orthologs_lookup_kogs;
my %transcripts_lookup_kogs; # hash of arrays of all transcripts in KOG
my %description_lookup_kogs;
my %kogs_unique;
<RERUN>;
while (<RERUN>)
	{
	chomp $_;
	my @data = split '\t', $_;
	my $transcript = $data[0];
	my $kog = $data[2];
	my $orthologs = $data[3];
	my $description = $data[4];
	
	$kogs_unique{$kog} = 1;
	
	$orthologs_lookup_kogs{$kog} = $orthologs;
	$description_lookup_kogs{$kog} = $description;
	push @{$transcripts_lookup_kogs{$kog}}, $transcript;
	}
close RERUN;

print OUT2 "OrthoDB_EuKOG\tDescription\t$species_out\_IDs\ttranscriptome_IDS\tOrtholog_count\tTrancript_count\tOrtho_type\($species_out\:Transcriptome\)\n";
foreach my $key (sort(keys %kogs_unique))
	{
	my @transcripts = @{${transcripts_lookup_kogs{$key}}};
	my $description = $description_lookup_kogs{$key};
	my $transcripts_out = join ',', @transcripts;
	my $transcript_count = scalar @transcripts;
	
	my $orthologs_out = $orthologs_lookup_kogs{$key};
	my $orthologs_count = 0;
	if ($orthologs_out eq "NA")
		{
		$orthologs_count = 0;
		}
	else {
		  $orthologs_count = ($orthologs_out =~ tr/,//); # count commas
		  $orthologs_count++; # add one to account for last ortholog e.g one,two = two values but one comma
			}
	
	my $type = "NA";
      
	if ($orthologs_count == 1 && $transcript_count == 1) {$type = "1:1";}
	elsif ($orthologs_count > 1 && $transcript_count == 1) {$type = "m:1";}
	elsif ($orthologs_count == 1 && $transcript_count > 1) {$type = "1:m";}
	elsif ($orthologs_count > 1 && $transcript_count > 1) {$type = "m:m";}
    elsif ($orthologs_count == 0 ) {$type = "unique";}  
	
	print OUT2 "$key\t$description\t$orthologs_out\t$transcripts_out\t$orthologs_count\t$transcript_count\t$type\n";
	
	}
close OUT2;




