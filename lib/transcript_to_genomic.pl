#!/usr/bin/env perl

# use lib "/dors/capra_lab/opt/ensembl_87/ensembl/modules";

use Bio::EnsEMBL::Registry;
my $registry = 'Bio::EnsEMBL::Registry';

my $instances_found = $registry->load_all(); # pens the config file in the environment variable ENSEMBL_CONF

if ($instances_found < 1) {
  print STDERR "Fundamental problems with EnsEMBL database connection setup.\n";
  print STDERR "Environment Variable ENSEMBL_REGISTRY=" . $ENV{"ENSEMBL_REGISTRY"} . "\n";
  exit 1;
}

# Command line argument: transcript ID
my $transcript_id = $ARGV[0];
# Open connection with local Ensembl database (read-only)
# Bio::EnsEMBL::Registry->load_registry_from_db(-host=>'vgi01.accre.vanderbilt.edu',-user=>'script_access',-pass=>'capralab');
# Bio::EnsEMBL::Registry->load_registry_from_db(-host =>'useastdb.ensembl.org',-user =>'anonymous');
# Create a transcript adaptor from the transcript ID
$transcript_adaptor = Bio::EnsEMBL::Registry->get_adaptor('human','core','Transcript');
$transcript = $transcript_adaptor->fetch_by_stable_id($transcript_id);

# Check that the transcript query was successful
if (!defined $transcript) {
  print STDERR "Not a valid human transcript ID: $transcript_id\n";
  exit 1; # clearly indicate that no transcript was found
}

# Extract general transcript data
my $q_strand = $transcript->strand();
my $t_stable_id = $transcript->stable_id();
my $gene = $transcript->get_Gene();
my $g_stable_id = $gene->stable_id();
my $slice = $transcript->slice();
my $s_name = $slice->seq_region_name();
my $s_start = $slice->start();
my $s_end = $slice->end();
my $translation = $transcript->translate();
my $p_stable_id = $transcript->translation()->stable_id();
my $biotype = $transcript->biotype();
if ($biotype ne "protein_coding") {
	print STDERR "$transcript_id is not protein coding.\n";
	exit 1;
}
if (!defined $translation) {
	print STDERR "No translation for $transcript_id at $s_name:$s_start-$s_end\n";
	exit 1; # clearly indicate that no transcript was found
}
my $peptide = $translation->seq();
my $peplength = length($peptide);

# For each amino acid in the peptide sequence, map to its genomic location
for ($i=1; $i<=$peplength; $i++) {
	my $pep_index = $i; # Name appropriately
	my @query = $transcript->pep2genomic($pep_index,$pep_index);
	my $q_start = @query[0]->start();
	my $q_end = @query[0]->end()+1; # Adjust for exclusive end
	# Catch erroneous Ensembl results
	if ($q_start < $s_start or $q_end > $s_end) {next;}
	# Pull the rescode from the peptide sequence
	my $amino_acid = substr($peptide,$pep_index-1,1);
	print "$t_stable_id\t$p_stable_id\t$g_stable_id\t$pep_index\t$amino_acid\t$q_start\t$q_end\tchr$s_name\t$q_strand\n";
}
