#!/usr/bin/env python
#
# Project        : PDBMap
# Filename       : PDBMapGnomad.py
# Author         : Chris Moth
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : chris.moth@vanderbilt.edu
# Date           : 2020-10-27
# =============================================================================#

""" Class PDBMapVep() returns the current set of Gnomad variants, for a given
    ENST ENSEMBL transcript identifier"""

# See main check for cmd line parsing
# import sys
import os
# import gzip
import logging
# import string
# import vcf
import pandas as pd
from typing import Dict, Union
import subprocess as sp
from lib.PDBMapGlobals import PDBMapGlobals
# from lib.PDBMapSQLdb import PDBMapSQLdb
from lib.PDBMapVEP import PDBMapVEP
from lib.PDBMapTranscriptEnsembl import PDBMapTranscriptEnsembl
# from lib import bed  # PyVCF emulator for BED files

LOGGER = logging.getLogger(__name__)
SQL_EXTENDED_LOCK_WAIT = "SET innodb_lock_wait_timeout=200;\n"


class PDBMapGnomad:
    """ Class PDBMapVep() returns the current set of Gnomad variants, for a given
        ENST ENSEMBL transcript identifier"""

    def __init__(self, config_dict: Dict[str, str] = None):
        """Construct a PDBMapGnomad() interface class from a dictionary containing:
           definition strings needed to run the variant effect predictor
           Also gnomad_dir and gnomad_filename_template
         """

        self._config_dict = config_dict
        if not self._config_dict:
            LOGGER.info("Initializing PDBMapVEP with PDBMapGlobals.config dictionary")
            self._config_dict = PDBMapGlobals.config
        if not self._config_dict or 'gnomad_dir' not in self._config_dict:
            raise Exception("gnomad_dir= is missing from the config file(s), or invalid")
        if not self._config_dict or 'gnomad_filename_template' not in self._config_dict:
            raise Exception("gnomad_filename_template= is missing from the config file(s), or invalid")

    def retrieve_gnomad_missense_variants(self,
                                          transcript: PDBMapTranscriptEnsembl,
                                          vep_echo_filename=None,
                                          output_directory=None) -> Union[pd.DataFrame, None]:
        """
        Given an Ensembl Transcript, determine the range of referenced genomic coordinates
        and query the Gnomad grch38 liftover database source files directly via tabix.
        Then, run those variants through the variant effect predictor to filter for missense variants.


        @transcript: Ensembl Transcript
        @vep_echo_filename: Optional

        @return: Dataframe of missense Gnomad variants variants
        """

        # Call the ENSEMBL PERL API to load genomic coordinates for the ensembl transcript
        (success, chrom, genomic_coordinates) = transcript.load_chromosome_location()
        if not success or len(genomic_coordinates) < 1:
            LOGGER.warning("Unable to load chromosome locations for %s", transcript.id)
            return None

        # The Gnomad supplied vcf filename omits chr in chrnn/chrX/chrY
        gnomad_vcf_file = os.path.join(
            self._config_dict['gnomad_dir'],
            self._config_dict['gnomad_filename_template'] % chrom[3:]
        )
        if not os.path.exists(gnomad_vcf_file):
            LOGGER.warning("Unable to open %s GNOMAD genomic coordinates for %s", gnomad_vcf_file, transcript.id)
            return 0;


        # The min and max coordinates fall out of the returned array from ENSEMBL PERL API
        # Because of splice variants, I take care to NOT assume the coordinate range is
        # always increases
        min_coordinate = min(genomic_coordinates[1][1],genomic_coordinates[1][2])
        max_coordinate = max(genomic_coordinates[1][1],genomic_coordinates[1][2])
        for index in range(2,len(genomic_coordinates)+1):
            min_coordinate = min(min_coordinate,genomic_coordinates[index][1],genomic_coordinates[index][2])
            max_coordinate = max(max_coordinate,genomic_coordinates[index][1],genomic_coordinates[index][2])

        LOGGER.info("Chrom: %s Coordinate range [%d,%d].  Excerpting vcf segment next with tabix.",
                    chrom, min_coordinate, max_coordinate)



        # Excerpt the tiny (by comparison) chrom region referenced by the Ensembl transcript
        gnomad_vcf_excerpt_filename = "%s_%d_%d_tabix.vcf" % (chrom[3:], min_coordinate, max_coordinate)
        if output_directory:
            gnomad_vcf_excerpt_filename = os.path.join(
                output_directory, gnomad_vcf_excerpt_filename)

        tabix_command = ['tabix', '-h',
                         gnomad_vcf_file,
                         "%s:%d-%d" % (chrom, min_coordinate, max_coordinate)]

        LOGGER.info("Running tabix with:\n%s" % ' '.join(tabix_command))
        with open(gnomad_vcf_excerpt_filename, 'wb') as gnomad_vcf_excerpt_f:
            completed_process = sp.run(tabix_command, stdout=gnomad_vcf_excerpt_f, timeout=240)
            completed_process.check_returncode()

        # Sanity check things a bit but counting the lines excerpted from the Chrom coordinate range.
        linecount = 0
        with open(gnomad_vcf_excerpt_filename, 'r') as gnomad_vcf_excerpt_f:
            for _ in gnomad_vcf_excerpt_f:
                linecount += 1

        LOGGER.info("%d lines in tabix excerpted file %s" % (linecount, gnomad_vcf_excerpt_filename))

        # This chunk of code is on hold, ebecause we get the MAFs from after the VEP call now.
        # MAFs are critical to dealing with Gnomad consequences.  We extract those from the tabix excerpt
        # gnomad_df = pd.DataFrame(columns=['gene', 'chrom', 'pos', 'maf'])
        # pdbmap_vep = PDBMapVEP(self._config_dict)
        # vcf_reader = pdbmap_vep.vcf_reader_from_file_without_vep(gnomad_vcf_excerpt_filename,prepend_chr=False)
        # for vcf_record in pdbmap_vep.yield_completed_vcf_records(vcf_reader):
        #        gnomad_df = gnomad_df.append(
        #            {'gene': vcf_record.CSQ[0]['SYMBOL'],
        #             'chrom': vcf_record.CHROM,
        #             'pos': vcf_record.POS,
        #             # 'mutation': "%s%s%s" % (CSQ['Ref_AminoAcid'], CSQ['Protein_position'], CSQ['Alt_AminoAcid']),
        #             'maf': float(vcf_record.INFO['AF']) # Minor allele frequency
        #             }, ignore_index=True)

        LOGGER.info("Running VEP to update transcript annotations. Will filter for %s" % transcript.id)
        pdbmap_vep = PDBMapVEP(self._config_dict)

        vcf_reader = pdbmap_vep.vcf_reader_from_file_supplemented_with_vep_outputs(
            gnomad_vcf_excerpt_filename,
            vep_echo_filename=vep_echo_filename,
            prepend_chr=False)

        # We now reduce the gnomad columns from over 4000 to these here:
        df_gnomad_missense = pd.DataFrame(columns=['gene', 'chrom', 'pos', 'transcript', 'Ref_AminoAcid',
                                   'Protein_position', 'Alt_AminoAcid', 'maf'])

        # Filter the VEP output records for missense variants, and build up a dataframe of all the collected
        # missense variants, which will be returned to the caller.
        for vcf_record in pdbmap_vep.yield_completed_vcf_records(vcf_reader):
            for CSQ in vcf_record.CSQ:
                vep_transcript = CSQ['Feature']
                if str(vep_transcript).strip() == transcript.id and str(CSQ['Consequence']).lower().find(
                        'missense') != -1:
                    df_gnomad_missense = pd.concat([df_gnomad_missense,pd.DataFrame(
                        [{'gene': CSQ['SYMBOL'],
                         'chrom': vcf_record.CHROM,
                         'pos': vcf_record.POS,
                         'transcript': CSQ['Feature'],
                         'Ref_AminoAcid': CSQ['Ref_AminoAcid'],
                         'Alt_AminoAcid': CSQ['Alt_AminoAcid'],
                         'Protein_position': CSQ['Protein_position'],
                         'maf': float(vcf_record.INFO['AF'])  # Minor allele frequency
                         }])], ignore_index=True)
        LOGGER.info("%d raw Gnomad missense variants excerpted from VCF fragment for %s" ,
                    len(df_gnomad_missense), transcript.id)
        return df_gnomad_missense
