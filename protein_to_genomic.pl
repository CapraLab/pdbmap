#!/usr/bin/perl

use lib "/usr/analysis/software/ensembl-api/ensembl/modules";
use lib "/usr/analysis/software/ensembl-api/ensembl-compara/modules";
use lib "/usr/analysis/software/ensembl-api/ensembl-variation/modules";
use lib "/usr/analysis/software/ensembl-api/ensembl-functgenomics/modules";

use Bio::EnsEMBL::Registry;

my $pdb = $ARGV[0];
my $chain = $ARGV[1];
my $unp = $ARGV[2];
my $species = $ARGV[3];

Bio::EnsEMBL::Registry->load_registry_from_db(-host=>'ensembldb.ensembl.org',-user=>'anonymous');

$transcript_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species,'Core','Transcript');

open(GenomicCoords,">>GenomicCoords.tab");
open(PDBTranscript,">>PDBTranscript.tab");

@transcripts = @{ $transcript_adaptor->fetch_all_by_external_name($unp,'uniprot%')};
foreach my $transcript (@transcripts) {
	my $t_stable_id = $transcript->stable_id();
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
		my $q_end = @query[0]->end();

		# If pep2genomic returns a gap, infer the strand
		# from the previous codon. i.e. don't overwrite
		# the strand variable from the previous iteration.
		if (ref @query[0] != "Bio::EnsEMBL::Mapper::Gap") {
			my $q_strand = @query[0]->strand();
		}
		my $amino_acid = substr($peptide,$i,1);
		my $pep_index = $i+1;
		print GenomicCoords "$t_stable_id\t$pep_index\t$amino_acid\t$q_start\t$q_end\t$s_name\t$q_strand\n";
	}

	print PDBTranscript "$pdb\t$chain\t$t_stable_id\n";
}
close(GenomicCoords);
close(PDBTranscript);
