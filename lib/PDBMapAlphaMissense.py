#!/usr/bin/env python
#
# Project        : PDBMap
# Filename       : PDBMapAlphaMissense.py
# Author         : Chris Moth
# Organization   : Meiler Lab
#                : Vanderbilt University Medical Center
# Email          : chris.moth@vanderbilt.edu
# Date           : 2019-0618-
# Description    : Class PDBMapTranscriptUniprot.py accesses and encapsulates
#                : everything we can know about one ENST....  transcript id
#=============================================================================#

# See main check for cmd line parsing
import sys
import os
import re
import subprocess
from lib.PDBMapSQLdb import PDBMapSQLdb
from typing import Optional

import logging
LOGGER = logging.getLogger(__name__)

class PDBMapAlphaMissense():
    """Encapsulate SQL calls to fetch alphamissense scores"""

    @staticmethod
    def score_from_uniprot_id(uniprot_id: str, protein_variant: str) -> Optional[float]:
        """Return alphamissense score or None if not found in database"""
        # Do NOT freeze SQL with a % wildcard with no string
        if not uniprot_id:
            return None
        query = """\
SELECT am_pathogenicity FROM AlphaMissense_aa_substitutions
WHERE uniprot_id = %(uniprot_id)s and protein_variant=%(protein_variant)s"""
        with PDBMapSQLdb() as sql_connection:
            sql_connection.execute(query,{'uniprot_id': uniprot_id.split('-')[0], 'protein_variant': protein_variant})
            row = sql_connection.fetchone()
            if not row or len(row) != 1:
                return None
            return float(row[0])
        

    @staticmethod
    def score_from_ENSEMBL_isoform_id(ensembl_enst_id: str, protein_variant: str) -> Optional[float]:
        """Return alphamissense score or None if not found in database"""
        # Do NOT freeze SQL with a % wildcard with no string
        if not ensembl_enst_id:
            return None
        query = """\
SELECT am_pathogenicity FROM AlphaMissense_isoforms_aa_substitutions
WHERE transcript_id like %(ensembl_enst_id)s and protein_variant=%(protein_variant)s"""
        with PDBMapSQLdb() as sql_connection:
            sql_connection.execute(query,{'ensembl_enst_id': ensembl_enst_id.split('.')[0] + '%', 'protein_variant': protein_variant})
            row = sql_connection.fetchone()
            if not row or len(row) != 1:
                return None
            return float(row[0])
        
# Main check
if __name__== "__main__":
  LOGGER.critical("Class definition. Should not be called from command line.")
  sys.exit(1)
