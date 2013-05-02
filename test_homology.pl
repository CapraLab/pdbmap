#!/usr/bin/perl

use lib "/usr/analysis/software/ensembl-api/ensembl/modules";
use lib "/usr/analysis/software/ensembl-api/ensembl-compara/modules";
use lib "/usr/analysis/software/ensembl-api/ensembl-variation/modules";
use lib "/usr/analysis/software/ensembl-api/ensembl-functgenomics/modules";

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_registry_from_db(-host=>'ensembldb.ensembl.org',-user=>'anonymous');

$gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor('sus_scrofa','Core','Gene');
@genes = @{ $gene_adaptor->fetch_all_by_external_name('Ssc.55804','unigene')};
foreach my $gene (@{genes}) {
	my $stable_id = $gene->stable_id();
	print "\n",$stable_id,":\n";
	$member_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','GeneMember');
	my $member = $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE',$stable_id);
	$transcript_adaptor = Bio::EnsEMBL::Registry->get_adaptor($member->taxon()->common_name,'Core','Transcript');
	@transcripts = @{ $transcript_adaptor->fetch_all_by_Gene($member->get_Gene())};
	my $orig_transcript = @transcripts[0];
	print "Original Transcript: ",$orig_transcript->stable_id(),"\n",$orig_transcript->translate()->seq(),"\n";
	my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','Homology');
	my $homologies = $homology_adaptor->fetch_all_by_Member($member);

	foreach $homology (@{$homologies}) {
		foreach $member (@{$homology->get_all_Members()}) {
			$taxon = $member->taxon();
			if (lc($taxon->common_name) eq 'human') {
				print "\nCommon Name: ", $taxon->common_name,"\n";
				print "Binomial: ", $taxon->binomial,"\n";
				print "Protein: ", $member->stable_id,"\n";
				print "Classification: ", $taxon->classification,"\n";
				my $transcript = $member->get_Transcript();
				print "\nHuman Transcript: ",$transcript->stable_id(),"\n";
				print $transcript->translate()->seq(),"\n";
			}
		}
	}
}