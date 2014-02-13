#!/usr/bin/perl

use lib "/usr/analysis/software/ensembl-api/ensembl/modules";
use lib "/usr/analysis/software/ensembl-api/ensembl-compara/modules";
use lib "/usr/analysis/software/ensembl-api/ensembl-variation/modules";
use lib "/usr/analysis/software/ensembl-api/ensembl-functgenomics/modules";

use Bio::EnsEMBL::Registry;

my $transcript_id = $ARGV[0];

#Bio::EnsEMBL::Registry->load_registry_from_db(-host=>'ensembldb.ensembl.org',-user=>'anonymous');
Bio::EnsEMBL::Registry->load_registry_from_db(-host=>'gwar-dev.mc.vanderbilt.edu',-user=>'script_access',-pass=>'bushlabrocks');

$transcript_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Human','Core','Transcript');

$transcript = $transcript_adaptor->fetch_by_stable_id($transcript_id);
my $q_strand = 0;
my $t_stable_id = $transcript->stable_id();
my $gene = $transcript->get_Gene();
my $g_stable_id = $gene->stable_id();
my $t_start = $transcript->start();
my $t_end = $transcript->end();
my @exons = @{ $transcript->get_all_Exons() };
my $cDNA = $transcript->spliced_seq();
my $slice = $transcript->slice();
my $s_name = $slice->seq_region_name();
my $s_start = $slice->start();
my $s_end = $slice->end();
my $peptide = $transcript->translate()->seq();
my $peplength = length($peptide);
my @locs = $transcript->pep2genomic(0,$peplength-1);

for ($i=0; $i<$peplength; $i++) {
	my @query = $transcript->pep2genomic($i,$i);
	my $q_start = @query[0]->start();
	# Adjust the end to match BED range definitions
	# (inclusive start, exclusive end)	
	my $q_end = @query[0]->end()+1;

	# Prevent erroneous Ensembl results
	if ($q_start < $s_start or $q_end > $s_end) {next;}

	# If pep2genomic returns a gap, infer the strand
	# from the previous codon. i.e. don't overwrite
	# the strand variable from the previous iteration.
	if (ref @query[0] != "Bio::EnsEMBL::Mapper::Gap") {
		my $q_strand = @query[0]->strand();
	}
	my $amino_acid = substr($peptide,$i,1);
	my $pep_index = $i+1;
	print "$t_stable_id\t$g_stable_id\t$pep_index\t$amino_acid\t$q_start\t$q_end\tchr$s_name\t$q_strand\t$peptide\n";
}
