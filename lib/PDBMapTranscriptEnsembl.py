#!/usr/bin/env python
#
# Project        : PDBMap
# Filename       : PDBMapTranscriptEnsembl.py
# Author         : Chris Moth, based on 2014 work of R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : chris.moth mike.sivley@vanderbilt.edu
# Date           : 2019-0618-
# Description    : Class PDBMapTranscriptEnsembl.py accesses and encapsulates
#                : everything we can know about one ENST....  transcript id
# =============================================================================#

# See main check for cmd line parsing
import sys
import os
import re
import subprocess
from Bio.Data import IUPACData
from lib import PDBMapTranscriptBase
from typing import Tuple

import logging

LOGGER = logging.getLogger(__name__)


class PDBMapTranscriptEnsembl(PDBMapTranscriptBase):
    """Encapsulate ENST... transcript IDs and their mappings to amino acid sequences
       and genomic locations"""
    # One global variable we'd like to tease out is the 
    # ENSEMBL_REGISTRY and specific dbname (containing Genome and Ensembl API version)

    _ensembl_registry = None
    _ensembl_api_version = None
    _ensembl_genome_version = None

    @classmethod
    def __init__ensembl_registry__(cls):
        if cls._ensembl_registry:  # Don't re-initialize
            return
        # One global variable we'd like to tease out is the
        # ENSEMBL_REGISTRY and specific dbname (containing Genome and Ensembl API version)
        cls._ensembl_registry = os.getenv("ENSEMBL_REGISTRY")
        if not cls._ensembl_registry:
            fail_message = "ENSEMBL_REGISTRY environment variable requried for PERL api configuration"
            LOGGER.critical(fail_message)
            sys.exit(fail_message)
        cls._ensembl_api_version = None
        cls._ensembl_genome_version = None
        if cls._ensembl_registry:
            with open(cls._ensembl_registry, 'r') as file:
                for line in file:
                    if '-dbname' in line and 'homo_sapiens_core_' in line and '=>' in line:
                        # Tease out ensembl API and Genome from the
                        # dbname => homo_sapiens_core_API_GENOME   string in the registry file
                        match_object = re.match(".*dbname.*homo_sapiens_core_(\\d*)_(\\d*).*", line)
                        if match_object:
                            cls._ensembl_api_version = int(match_object.group(1))
                            cls._ensembl_genome_version = int(match_object.group(2))

        if not (cls._ensembl_api_version and cls._ensembl_genome_version):
            LOGGER.critical(
                "registry file %s lacks a -dbname => 'homo_sapiens_core_API_GENOME' entry" % cls._ensembl_registry)
            sys.exit(1)

    def __init__(self, ensembl_ENST_transcript_id):
        """An ENSTnnnnnn ENSEMBL transcript ID is required to instantiate"""
        # Define transcript, gene, and sequence
        if not PDBMapTranscriptEnsembl._ensembl_registry:
            PDBMapTranscriptEnsembl.__init__ensembl_registry__()

        self._ensembl_ENST = ensembl_ENST_transcript_id
        self.ensembl_ENSG = None  # ENSG* gene identifer mapped to this ENST transcript
        self.ensembl_ENSP = None  # ENSP* protein id mapped to this ENST transcript
        self.ensembl_chromosome = None  # Chromosome 'chr*' returned from the ensembl API
        self.ensembl_strand = None  # 1 or -1, depending on forward or reverse strand

        # genome_locations is a dictionary leyed of transcript location
        # Sequence | key:   transcript seqid
        # Sequence | value: (aa_letter,start,end,strand)
        self.ensembl_sequence = {}  # ensembl_sequence[store genomic sequence keyed on genome location
        super().__init__(None, None)  # Set self._aa_seq and self._uniparc_id to None for starters

    @property
    def id(self):
        return self._ensembl_ENST

    @property
    def aa_seq(self):
        if not self._aa_seq:
            (result, aa_seq) = self.load_aa_seq_from_ENSEMBL()
            if not result:  # Then aa_seq contains an error message
                msg = "ENSMBL API failed to return aa_seq for %s: %s" % (self._ensembl_ENST, aa_seq)
                LOGGER.critical(msg)
                return None
        return self._aa_seq

    def is_standard_chromosome(self):
        if not self.ensembl_chromosome:
            LOGGER.error("Chromosome chr* string must be first loaded from ENSMBL Api call")
            return False

        chrom_suffix = ''
        if self.ensembl_chromosome[0:3] == 'chr':
            chrom_suffix = self.ensembl_chromosome[3]
        # Make sure the chromosome is chr1 to chr22 OR chrX, chrY, chrMT
        # If transcript on a haplotype chromosome, ignore
        return ((chrom_suffix.isdigit() and int(chrom_suffix) <= 22 and chrom_suffix[0] >= '1') or
                (chrom_suffix in ['X', 'Y', 'MT']))

        #    LOGGER.warn("Ignoring non-standard chromosome %s returned from \'%s\' command output"%(chrom,cmd))

    def _run_perl_script(self, script_filename: str) -> subprocess.CompletedProcess:
        """
        Run one of the ...pl perl scripts that interfaces to the ENSEMBL Perl API.
        The one argument is simply the ENST transcript id
        """
        cmd_list = ["%s" % script_filename, self._ensembl_ENST]
        cmd_string = ' '.join(cmd_list)
        LOGGER.info("Launching PERL script to access ENSEMBL: %s" % cmd_string)
        completed_process = subprocess.run(cmd_list, encoding='utf-8', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if completed_process.returncode != 0:
            LOGGER.warning("Exit Code %d returned from: %s\nstderr: %s\nstdout:%s" % (
                completed_process.returncode,
                cmd_string,
                completed_process.stderr.rstrip(), completed_process.stdout.rstrip()))
        return completed_process

    # Quick sanity check to ensure returned aa_seq makes some sense
    def _aa_seq_check(self) -> Tuple[bool, str]:
        if not self._aa_seq:  # MAjor fail of calling software
            LOGGER.critical("AAseq not yet loaded")
            sys.exit(1)

        valid_aa_letters = set(IUPACData.protein_letters)  # Most standard 20
        valid_aa_letters.add('U')  # Deal with Selenomet
        all_letters_ok = all(aa in valid_aa_letters for aa in self.aa_seq)
        if all_letters_ok:
            return True, self.aa_seq
        return False, "Bad AA letters are %s" % str([aa for aa in self.aa_seq if aa not in valid_aa_letters])

    def load_aa_seq_from_ENSEMBL(self) -> Tuple[bool, str]:
        """Run lib/transcript_to_AAseq.pl to retrieve aa_seq from ENST... id"""
        if self._aa_seq:
            LOGGER.info("aa_seq already loaded.  Returning from load_aa_seq_from_ENSEMBL immediately")
            return True, self.aa_seq
        completed_process = self._run_perl_script("transcript_to_AAseq.pl")
        if completed_process.returncode != 0:
            return False, completed_process.stderr.rstrip() + "(GRCh%s v %s)" % (
                PDBMapTranscriptEnsembl._ensembl_genome_version, PDBMapTranscriptEnsembl._ensembl_api_version)
        self._aa_seq = completed_process.stdout.rstrip()
        return self._aa_seq_check()

    def load_chromosome_location(self):
        """ Load genomic locations for the _ensembl_ENST via the ENSEMBL PERL API """
        if self.ensembl_sequence and self.ensembl_chromosome:
            LOGGER.info("chromosome location already loaded.  Returning from load_chromosome_location immediately")
            return True, self.ensembl_chromosome, self.ensembl_sequence

        completed_process = self._run_perl_script("transcript_to_genomic.pl")
        if completed_process.returncode != 0:
            return False, completed_process.stderr.rstrip(), None
        lineno = 1
        for line in completed_process.stdout.split('\n'):
            if len(line) < 1:  # It's nothing
                continue
            LOGGER.debug("Line %3d: %s", lineno, line)
            lineno += 1
            if line.startswith('#'):  # It's a comment
                continue
            fields = line.split('\t')
            if len(fields) < 2:  # Let's not be upset with a blank line
                continue
            if len(fields) < 9:  # Let's report on lines that don't make sense
                LOGGER.warning("Incomplete line returned from transcript_to_genomic.pl:\n[%s]" % line)
                continue

            # The left column of returned values should match our requested ENST transcript id
            transcript = fields[0]
            assert transcript == self._ensembl_ENST

            # The ENSP Protein identifier shold be the same for all liness returned from the PERL API
            if self.ensembl_ENSP:
                assert (self.ensembl_ENSP == fields[1].strip())
            else:
                self.ensembl_ENSP = fields[1].strip()

            # The ENSP Gene identifier shold be the same for all lines returned by the PERL API
            if self.ensembl_ENSG:
                assert (self.ensembl_ENSG == fields[2].strip())
            else:
                self.ensembl_ENSG = fields[2].strip()

            # The ENSEMBL chrnn should be the same for all returned lines from the PERL API
            # And, it must be one of the "standard" human chromosomes
            if self.ensembl_chromosome:
                assert (self.ensembl_chromosome == fields[7].strip())
            else:
                # If the returned chromosome is not a "standard" human , then we need to return failure to the caller
                self.ensembl_chromosome = fields[7].strip()
                if not self.is_standard_chromosome():
                    bad_chr_format_message = "Unrecognized chr:%s  in line %d: %s" % (
                        self.ensembl_chromosome,
                        lineno,
                        line)

                    LOGGER.warning(bad_chr_format_message)
                    return False, bad_chr_format_message, None

            if self.ensembl_strand:
                assert (self.ensembl_strand == int(fields[8]))
            else:
                self.ensembl_strand = int(fields[8])

            transcript_index = int(fields[3])  # 1 to N - the index of the aa transcript
            aa_letter = fields[4].upper()
            if aa_letter not in IUPACData.protein_letters:
                aa_letter = 'X'  # replace non-standard amino acids with X
            start = int(fields[5])  # Genomic position on the chromosome where this AA starts
            end = int(fields[6])  # Genomic position on the chromosome where this AA ends

            # It is convenient to build up the AA sequence without doing the other PERL API call
            self.ensembl_sequence[transcript_index] = (aa_letter, start, end)

        if not self.ensembl_sequence:
            return False, "No interpretable data returned from ENSEMBL PERL API transcript_to_genomic.pl: %s" % \
                        completed_process.stdout.rstrip(), None

        residue_numbers = sorted(self.ensembl_sequence.keys())
        if residue_numbers[0] != 1:
            return False, "Amino Acid for residue #1 was not returned by transcript_to_genomic.pl", None
        if residue_numbers[-1] != len(residue_numbers):
            return False, "Last residue number returned by transcript_to_genomic.pl does not match residue count", None

        # Update the internal aa_seq string with these aa_letters
        self._aa_seq = ''.join([self.ensembl_sequence[transcript_index][0] for transcript_index in residue_numbers])

        return True, self.ensembl_chromosome, self.ensembl_sequence


# Main check
if __name__ == "__main__":
    LOGGER.critical("Class definition. Should not be called from command line.")
    sys.exit(1)
