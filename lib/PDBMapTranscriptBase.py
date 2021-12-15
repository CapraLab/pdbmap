#!/usr/bin/env python
#
# Project        : PDBMap
# Filename       : PDBMapTranscriptBase.py
# Author         : Chris Moth, based on 2014 work of R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : chris.moth@vanderbilt.edu, formerly mike.sivley@vanderbilt.edu
# Date           : 2019-0618-
# Description    : Class PDBMapTranscriptBase.py accesses and encapsulates
#                : everything common to transcript child classes (Ensembl, Uniprot, Fasta)
# =============================================================================s

import sys
from typing import Tuple
from typing import Optional
from Bio.Data import IUPACData
from lib import PDBMapSQLdb

import logging

LOGGER = logging.getLogger(__name__)


class PDBMapTranscriptBase:
    """Base class for all PDBMapTranscript* Amino Acid (aa_seq) transcript representations """

    def __init__(self, aa_seq=None, uniparc_id=None):
        # Define transcript, gene, and sequence
        self._uniparc_id = uniparc_id  # Later populate with uniparc_id
        self._aa_seq = aa_seq  # Single letter amino acid letters of transcript

    # Quick sanity check to ensure returned aa_seq makes some sense
    @staticmethod
    def _aa_seq_check(aa_seq: str) -> Tuple[bool, Optional[str]]:
        all_letters_ok = all(aa in IUPACData.protein_letters for aa in aa_seq)
        if all_letters_ok:
            return True, None
        return False, "Bad AA letters are %s" % str([aa not in IUPACData.protein_letters for aa in aa_seq])

    @property
    def md5sum(self):
        import hashlib
        return hashlib.md5(self.aa_seq.encode('UTF-8')).hexdigest()

    @property
    def aa_seq(self) -> str:
        if self._aa_seq:
            return self._aa_seq

        if (not self._aa_seq) and self._uniparc_id:
            (success, error_msg) = self.load_aa_seq_from_uniparc()
            if not success:
                raise Exception("Transcript Base class cannot return aa_seq from uniparc ID %s\nReason: %s" % (
                    self._uniparc_id, error_msg))
            if success:
                return self._aa_seq

        raise Exception("Transcript Base class cannot return aa_seq when neither aa_seq nor uniparc_id is provided")

    @property
    def len(self):
        return len(self.aa_seq)

    def load_uniparc_from_md5sum(self) -> Tuple[bool, str]:
        query = """
SELECT uniparc,fasta Uniparc
where md5sum = %(md5sum)s"""
        with PDBMapSQLdb() as db:
            db.execute(query, {'md5sum': self.md5sum})
            row = db.fetchone()
            if not row or len(row) != 2:
                return False, "Query failed to return row with 2-tuple"
            self._uniparc_id = str(row[0]).strip()
            self._aa_seq = str(row[1]).strip()
            return True, self._uniparc_id

    @property
    def uniparc_id(self):
        if self._uniparc_id:
            return self._uniparc_id
        return self.load_uniparc_from_md5sum()

    def load_aa_seq_from_uniparc(self):
        query = """
SELECT md5sum,fasta FROM Uniparc
where uniparc = %(uniparc_id)s"""
        with PDBMapSQLdb() as db:
            db.execute(query, {'uniparc_id': self._uniparc_id})
            row = db.fetchone()
            if not row or len(row) != 2:
                return False, "Query failed to return row with 3-tuple"
            # self._md5sum = str(row[0]).strip()
            self._aa_seq = str(row[1]).strip()
        return True, self._aa_seq

    @staticmethod
    def describe_transcript_differences(t1, t2):
        """In English, succinctly describe AA seq differences in 2 PDBMapTranscripts"""
        if len(t1.aa_seq) != len(t2.aa_seq):
            return "Transcript 1 (%s %s) has len=%d.  Transcript 2 (%s %s) has len=%d" % (
                type(t1), t1.id, len(t1.aa_seq),
                type(t2), t2.id, len(t2.aa_seq))

        differences = ""
        for tpos in range(len(t1.aa_seq)):
            if t1.aa_seq[tpos] != t2.aa_seq[tpos]:
                differences += "%4d %s vs. %s\n" % (tpos + 1, t1.aa_seq[tpos], t2.aa_seq[tpos])
        return differences

    @staticmethod
    def analyze_transcript_differences(t1, t2):
        """Return # of single point variants, and percent of mismatches, between 2 PDBMapTranscripts"""
        if len(t1.aa_seq) != len(t2.aa_seq):
            # If transcripts are not the same length, return the longer length, and 100% mismatch
            return max(len(t1.aa_seq), len(t2.aa_seq)), 1.0

        differences = 0
        for tpos in range(len(t1.aa_seq)):
            if t1.aa_seq[tpos] != t2.aa_seq[tpos]:
                differences += 1
        return differences, float(differences) / float(len(t1.aa_seq))

    # Main check


if __name__ == "__main__":
    LOGGER.critical("Class definition. Should not be called from command line.")
    sys.exit(1)
