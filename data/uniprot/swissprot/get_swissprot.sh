#!/bin/sh

# Download all human swiss-prot entries
wget -q -N --reject -nH -nd --timeout=100000 ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz

# Reduce to human-only proteins
./filter_uniprot.py uniprot_sprot.dat.gz > uniprot_sprot_human.dat

# Remove the gzip archive
rm -f uniprot_sprot.dat.gz
