#!/usr/bin/perl

use lib "/usr/analysis/software/ensembl-api/ensembl/modules";
use lib "/usr/analysis/software/ensembl-api/ensembl-compara/modules";
use lib "/usr/analysis/software/ensembl-api/ensembl-variation/modules";
use lib "/usr/analysis/software/ensembl-api/ensembl-functgenomics/modules";

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_registry_from_db(-host=>'gwar-dev.mc.vanderbilt.edu',-user=>'script_access',-pass=>'bushlabrocks');
$slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "slice" );
$transcript_adaptor = Bio::EnsEMBL::Registry->get_adaptor("human" ,'core','transcript');

my $bed =  $ARGV[0];
open(FIN,$bed);
while (<FIN>) {
	chomp;
	($chr,$start,$end,$name) = split("\t");

	$slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $start, $end );
	if (ref($slice) eq Bio::EnsEMBL::Slice) {
		@transcripts = @{ $transcript_adaptor->fetch_all_by_Slice($slice) };
		foreach my $transcript (@transcripts) {
			my $trans_id = $transcript->stable_id();
			print STDOUT "$chr\t$start\t$end\t$name\t$trans_id\n";
		}
	}
}

