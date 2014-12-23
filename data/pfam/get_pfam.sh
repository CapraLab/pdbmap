#!/bin/sh

# Downloads the PFAM-to-PDB mapping,
# Online repo updated daily

# Download the mapping file
wget -N --reject -nH -nd --timeout=100000 ftp://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/pdb_pfam_mapping.txt
