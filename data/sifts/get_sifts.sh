#!/bin/sh
#
# Project        : PDBMap
# Filename       : get_sifts.sh
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-17
# Description    : Downloads a local copy of SIFTS pdb->uniprot residue mapping

wget -q -N --reject -nH -nd --timeout=100000 ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz
gunzip -c pdb_chain_uniprot.tsv.gz > pdb_chain_uniprot.tsv
