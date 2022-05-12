#!/usr/bin/env python
#
# Project        : PDBMap
# Filename       : PDBMapSIFTSdb.py
# Author         : 2019 Chrs Moth
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : chris.moth@vanderbilt.edu
# Date           : 2019-12-19

# See main check for cmd line parsing
import sys,os,csv
from collections import OrderedDict
import logging
from Bio import pairwise2
# unused and deprecated: from Bio.SubsMat import MatrixInfo as matlist
from Bio.PDB import Structure
from Bio.SeqUtils import seq1
from lib import PDBMapSQLdb
from typing import Dict,List,Tuple

LOGGER = logging.getLogger(__name__)

class PDBMapSIFTSdb():
    """ Provide SQL queries and other processing related to SIFTS-sourced tables in pdbmap_v*
        Before this class can be used, SQL tables must have first been created with
        the create_schema_sifts*.sql files.  Population (INSERTs) into the tables must
        have been completed with ../scripts/sifts_parser.py

        As documented in those prerequisites, SIFTS source data is on the web under
        https://www.ebi.ac.uk/pdbe/docs/sifts 
    """

    @staticmethod
    def sifts_best_unps(structure: Structure) -> Dict[str, str]:
        """Return a dictionary mapping each chain ID in a the Structure
           to the uniprot identifier that SIFTS
           has determined to be 'best' for each chain"""
    
        chain_to_best_unp = {}
    
        query = "SELECT pdbid,uniprot_acc,mapping_pdb_chain FROM sifts_mappings_pdb_uniprot_best_isoforms where pdbid=%(pdbid)s ORDER BY mapping_pdb_chain"
        with PDBMapSQLdb() as db:
            db.activate_dict_cursor()
            db.execute(query,{'pdbid': structure.id})
            for row in db.fetchall():
                chain_to_best_unp[row['mapping_pdb_chain']] = row['uniprot_acc']
    
        return chain_to_best_unp

    @staticmethod
    def pdbs_chains_for_all_unp_isoforms(unp: str) -> List[Tuple[str,str]]:
        query = "SELECT DISTINCT pdbid,mapping_pdb_chain FROM sifts_mappings_pdb_uniprot_all_isoforms where uniprot_acc LIKE %s ORDER BY pdbid, mapping_pdb_chain"
        pdbs_chains = []
        with PDBMapSQLdb() as db:
            db.activate_row_cursor()
            db.execute(query,(unp.split('-')[0]+'%',))
            for row in db.fetchall():
                pdbs_chains.append((row[0],row[1]))
    
        return pdbs_chains

    @staticmethod
    def pdbs_chains_coverage_for_unp(
        unp: str, 
        unp_is_canonical: bool,
        max_unp_start:int = 2000000000, 
        min_unp_end:int = 1) -> List[Dict[str,str]]:
        """Query sifts for pdb IDs and chains that are aligned to a uniprot identifer
           unp:              Uniprot identifier.  Ideally supplied with a dash
           unp_is_canonical: If True, both dashed and un-dashed unps will be searched
           max_unp_start,min_unp_end: Change from defaults to restrict query to certain coverage"""

        query = """
SELECT * from 
(SELECT 
pdbid,mapping_pdb_chain,
MIN(mapping_unp_start) AS min_unp_start,
MAX(mapping_unp_end)  AS max_unp_end
FROM pdbmap_v14.sifts_mappings_pdb_uniprot_all_isoforms WHERE """

        if unp_is_canonical and '-' in unp: # Ask SIFTS to search for dashless canoncical unp form as well
            query += "(uniprot_acc = %s OR uniprot_acc = %s)"
            values_tuple = (unp,unp.split('-')[0],max_unp_start,min_unp_end)
        else: # We have an isoform-specific unp with a dash, or canonical unp without a dash
            query += "(uniprot_acc = %s)"
            values_tuple = (unp,max_unp_start,min_unp_end)
        query += """
GROUP BY pdbid, mapping_pdb_chain
) minMaxTranscriptResiduesCovered
WHERE minMaxTranscriptResiduesCovered.min_unp_start <= %s and minMaxTranscriptResiduesCovered.max_unp_end >= %s
ORDER BY pdbid, mapping_pdb_chain"""
        pdbs_chains = []
        with PDBMapSQLdb() as db:
            db.activate_dict_cursor()
            db.execute(query,values_tuple)
            for row in db.fetchall():
                pdbs_chains.append(row)

        return pdbs_chains

