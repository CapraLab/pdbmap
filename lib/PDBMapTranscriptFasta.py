#!/usr/bin/env python
#
# Project        : PDBMap
# Filename       : PDBMapTranscriptBase.py
# Author         : Chris Moth, based on 2014 work of R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : chris.moth mike.sivley@vanderbilt.edu
# Date           : 2019-0618-
# Description    : Class PDBMapTranscriptBase.py accesses and encapsulates
#                : everything we can know about one ENST....  transcript id
#=============================================================================#

# See main check for cmd line parsing
import sys
import os
import re
import subprocess
from Bio import SeqIO
from Bio.Data import IUPACData
from lib import PDBMapSQLdb
from lib import PDBMapTranscriptBase

import logging
LOGGER = logging.getLogger(__name__)

class PDBMapTranscriptFasta(PDBMapTranscriptBase):
    """Like PDBMapTranscriptBase, but construts with a fasta string"""

    @property
    def id(self):
        return self._fasta_id

    def __init__(self,fasta_string):
        if not fasta_string:
            raise Exception("You must supply a fasta_string to create this transcript")

        if fasta_string[0] != '>':
            LOGGER.info("Treating fasta_string %s as filename"%fasta_string)
            with open(fasta_string,'r') as f:
                fasta_string = f.read()

        import io
        fasta_handle = io.StringIO(fasta_string)
        for record in SeqIO.parse(fasta_handle, 'fasta'):
            super().__init__(record.seq)
            self._fasta_id = record.id
            break

        assert self._fasta_id and self.aa_seq,"Fasta string lacks id and/or sequence:\n%s"%fasta_string
