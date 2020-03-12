#!/usr/bin/env perl 

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

if (1) {
my $instances_found = $registry->load_all(); # pens the config file in the environment variable ENSEMBL_CONF

if ($instances_found < 1) {
 print STDERR "Fundamental problems with EnsEMBL database connection setup.\n";
 print STDERR "Environment Variable ENSEMBL_REGISTRY=" . $ENV{"ENSEMBL_REGISTRY"} . "\n";
 exit 1;
}
}



# Command line argument: transcript ID
my $transcript_id = $ARGV[0];
# Open connection with local Ensembl database (read-only)
# Bio::EnsEMBL::Registry->load_registry_from_db(-host=>'vgi01.accre.vanderbilt.edu',-user=>'script_access',-pass=>'capralab'); #,verbose=>1);
# Bio::EnsEMBL::Registry->load_registry_from_db(-host =>'ensembldb.ensembl.org',-user =>'anonymous', -port => 3337); #, verbose => 1);
# Create a transcript adaptor from the transcript ID
$transcript_adaptor = Bio::EnsEMBL::Registry->get_adaptor('human','core','Transcript');
$transcript = $transcript_adaptor->fetch_by_stable_id($transcript_id);
print $transcript->translate()->seq()."\n";

if (!defined $transcript) {
  print STDERR "Not a valid human transcript ID: $transcript_id\n";
  exit 1; # clearly indicate that no transcript was found
}

my $trv_adaptor = $registry->get_adaptor('human','variation','TranscriptVariation');
if (!defined $trv_adaptor) {
  print STDERR "trv_adaptor fail";
  exit 1; # clearly indicate that no transcript was found
}
my $transcript_list = [$transcript];
my @trvs = @{ $trv_adaptor->fetch_all_by_Transcripts([$transcript]) };

while (my $tv = shift @trvs) {
  my $tvas = $tv->get_all_alternate_TranscriptVariationAlleles();
  foreach my $tva (@{$tvas}) {
  print "Variant ", $tv->variation_feature->variation_name, " allele ", $tva->variation_feature_seq;

  }

}



print $transcript->translate()->seq()."\n";
