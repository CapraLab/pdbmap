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

# Download the mapping file
wget -q -N --reject -nH -nd --timeout=100000 ftp://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/pdb_pfam_mapping.txt
