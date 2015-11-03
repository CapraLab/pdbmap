#!/bin/sh
#
# Project        : PDBMap
# Filename       : get_pfam.sh
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-17
# Description    : Downloads the PFAM-to-PDB mapping
#
# Once uploaded to mysql, use this query to annotate the
# pfam entries with UniProt accession numbers:
#  UPDATE pfam a
#  left join Chain b
#  on a.pdbid=b.structid and a.chain=b.chain and b.label='uniprot-pdb' and biounit=0 and model=0
#  SET a.unp=b.unp;

# Download the mapping file
#wget -q -N --reject -nH -nd --timeout=100000 ftp://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/pdb_pfam_mapping.txt
wget -q -N --reject -nH -nd --timeout=100000 -O pdb_pfam_mapping.txt http://www.rcsb.org/pdb/rest/hmmer?file=hmmer_pdb_all.txt
