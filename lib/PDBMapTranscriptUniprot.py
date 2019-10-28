#!/usr/bin/env python
#
# Project        : PDBMap
# Filename       : PDBMapTranscriptUniprot.py
# Author         : Chris Moth, based on 2014 work of R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : chris.moth mike.sivley@vanderbilt.edu
# Date           : 2019-0618-
# Description    : Class PDBMapTranscriptUniprot.py accesses and encapsulates
#                : everything we can know about one ENST....  transcript id
#=============================================================================#

# See main check for cmd line parsing
import sys
import os
import re
import subprocess
from Bio.Data import IUPACData
from lib.PDBMapSQLdb import PDBMapSQLdb
from lib.PDBMapTranscriptBase import PDBMapTranscriptBase

import logging
LOGGER = logging.getLogger(__name__)

class PDBMapTranscriptUniprot(PDBMapTranscriptBase):
    """Encapsulate Uniprot transcript IDs and their mappings to amino acid sequences
       via uniparc_id from the idmapping file, or sql queries"""

    def __init__(self,uniprot_id):
        # Define transcript, gene, and sequence
        self._uniprot_id = uniprot_id
        super().__init__(None,None) # Set self._aa_seq and self._uniparc_id to None for starters

    def __repr__(self):
        return self._uniprot_id

    @property
    def uniprot_id(self):
        return self._uniprot_id

    @property
    def id(self):
        return self._uniprot_id

    def load_aa_seq_from_sql(self):
        if not self._aa_seq:
            query = """
SELECT unp,uniparc,fasta FROM Idmapping 
LEFT JOIN Uniparc on Uniparc.uniparc = Idmapping.ID  
where Id_Type = 'UniParc' and unp = %(unp)s"""
            with PDBMapSQLdb() as sql_connection:
                sql_connection.execute(query,{'unp': self.uniprot_id})
                row = sql_connection.fetchone()
                if not row or len(row) != 3:
                    return (False,"For bad unp %s, SQL query failed"%self.uniprot_id)
                self._uniparc_id = str(row[1]).strip()
                self._aa_seq = str(row[2]).strip()
        return (True,self.aa_seq)

    @property
    def aa_seq(self):
        if not self._aa_seq:
            (result,aa_seq) = self.load_aa_seq_from_sql()
            if not result: # Then aa_seq contains an error message
                LOGGER.critical(aa_seq)
                raise LookupError(aa_seq)
        return self._aa_seq
    
# Main check
if __name__== "__main__":
  LOGGER.critical("Class definition. Should not be called from command line.")
  sys.exit(1)
