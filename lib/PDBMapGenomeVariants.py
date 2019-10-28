#!/usr/bin/env python
#
# Project        : PDBMap
# Filename       : PDBMapGenomeVariants.py
# Author         : Chris Moth, based on 2014 work of R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : chris.moth mike.sivley@vanderbilt.edu
# Date           : 2019-0618-
# Description    : Class PDBMapGenomeVariants accesses and capsulates
#                : everything we can know about variant data sets such as
#                  clinvar, exac, COSMIC, etc
#=============================================================================#

# See main check for cmd line parsing
import sys
import os
import re
import subprocess
from Bio.Data import IUPACData
from lib import PDBMapSQLdb

import logging
LOGGER = logging.getLogger(__name__)

class PDBMapGenomeVariants():
    """Set of functions to query GenomicConsequence and other SQL tables
       of PDBMap tha store clinvar/COSMIC/etc variants"""

    @staticmethod
    def query_missense_variants(ensembl_transcripts,label: str) -> (bool,list):
        """Return list of 3-tuples (transcript_pos, ref_amino_acid, alt_amino_acid)
           Parameters:
                ensembl_transript: ENST_....  transcript identifier(s)
                label: cosmic/CLINVAR/exac etc label for SQL retrieval"""
        #                               0           1              2
        query  = """SELECT distinct protein_pos,ref_amino_acid,alt_amino_acid FROM pdbmap_v13.GenomicConsequence\nWHERE """

        if isinstance(ensembl_transcripts,list) and len(ensembl_transcripts) == 1:
            ensembl_transcripts = ensembl_transcripts[0]

        if isinstance(ensembl_transcripts,str):
            query += "transcript='%s' "%ensembl_transcripts
        else:
            assert isinstance(ensembl_transcripts,list), "ensembl_transcripts parameter must either be a single ENST... string or list of such strings"
            query += "transcript in (" + ','.join("'{0}'".format(enst_transcript) for enst_transcript in ensembl_transcripts) + ') '

        query += """ and label=%(label)s and consequence LIKE '%%missense_variant%%' and length(ref_amino_acid)=1 and length(alt_amino_acid)=1 """

        variant_list = []
        with PDBMapSQLdb() as sql_connection:
            sql_connection.execute(query,{'label': label})
            rows = sql_connection.fetchall()

        for row in rows:
           variant_list.append(  (int(row[0]),row[1],row[2]))
            
        return (True,variant_list)
    
# Main check
if __name__== "__main__":
  LOGGER.critical("Class definition. Should not be called from command line.")
  sys.exit(1)
