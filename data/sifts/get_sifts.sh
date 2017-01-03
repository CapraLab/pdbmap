#!/bin/sh
#
# Project        : PDBMap
# Filename       : get_sifts.sh
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2017-01-03
# Description    : Downloads a local copy of SIFTS pdb->uniprot residue mapping

#wget -q -N --reject -nH -nd --timeout=100000 ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz
#gunzip -c pdb_chain_uniprot.tsv.gz > pdb_chain_uniprot.tsv

wget -r -N --no-parent --reject -nH -nd --timeout=100000 --tries=100 -P xml ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml

