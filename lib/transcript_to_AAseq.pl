#!/usr/bin/env perl 

use Bio::EnsEMBL::Registry;

# Command line argument: transcript ID
my $transcript_id = $ARGV[0];
# Open connection with local Ensembl database (read-only)
Bio::EnsEMBL::Registry->load_registry_from_db(-host=>'vgi01.accre.vanderbilt.edu',-user=>'script_access',-pass=>'capralab'); #,verbose=>1);
# Bio::EnsEMBL::Registry->load_registry_from_db(-host =>'useastdb.ensembl.org',-user =>'anonymous', verbose => 1);
# Create a transcript adaptor from the transcript ID
$transcript_adaptor = Bio::EnsEMBL::Registry->get_adaptor('human','core','Transcript');
$transcript = $transcript_adaptor->fetch_by_stable_id($transcript_id);

if (!defined $transcript) {
  print STDERR "Not a valid human transcript ID: $transcript_id\n";
  exit 1; # clearly indicate that no transcript was found
}


print $transcript->translate()->seq()."\n";
