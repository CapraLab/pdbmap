#!/usr/bin/env python
# Project        : PDBMap
# Filename       : PDBMapModbase
# Author         : Chris Moth
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2017-09-07
# Description    : PDBMapStructure equivalent for SwissModel models
#                  based on Mike Sivley's PDBMapModel.py

import argparse
# import json

import sys,os,re
from typing import List
from typing import Dict
from typing import Tuple


from lib.PDBMapSQLdb import PDBMapSQLdb
from lib.PDBMapTranscriptBase import PDBMapTranscriptBase

class PDBMapModbase2020():
    def __init__(self,config_dict):
        assert 'modbase2020_dir' in config_dict,"The key modbase2020_dir is required in the passed dictionar - but not found"
        self.modbase2020_dir = config_dict['modbase2020_dir']

    def ensp2modelids(self,ensp_id: str) -> List[str]:
        """Return a list of (model_id,database_id) tuples for an ENSP id"""
        ensp_model_matches = []
        with PDBMapSQLdb() as db:
            rows_selected =  db.execute("SELECT model_id,database_id FROM Modbase2020 WHERE database_id LIKE %(ensp_id)s",{'ensp_id': ensp_id+'%%'})
            if rows_selected:
                row_tuples = db.fetchall()
                ensp_model_matches = [(row[0],row[1]) for row in row_tuples]
            assert rows_selected == len(ensp_model_matches)
        return ensp_model_matches if ensp_model_matches else []

    def transcript2modelids(self,transcript: PDBMapTranscriptBase) -> List[Tuple[str,str]]:
        """Return a list of (model_id,database_id) tuples for any transcript, via its amino acid sequence"""
        transcript_model_matches = []
        with PDBMapSQLdb() as db:
            rows_selected =  db.execute("SELECT model_id,database_id FROM Modbase2020 WHERE seq_id = %s",(transcript.md5sum + transcript.aa_seq[0:4] + transcript.aa_seq[-4:],))
            if rows_selected:
                row_tuples = db.fetchall()
                transcript_model_matches = [(row[0],row[1]) for row in row_tuples]
            assert rows_selected == len(transcript_model_matches)
        return transcript_model_matches if transcript_model_matches else []

    def transcript2summary_rows(self,transcript: PDBMapTranscriptBase) -> List[Dict]:
        """Return a list of Dictionaries that contain the rows of the Modbase summary file that match a sequence of interest"""
        modbase2020_summary_rows = []
        with PDBMapSQLdb() as db:
            db.activate_dict_cursor()
            rows_selected_count =  db.execute(
                "SELECT * FROM Modbase2020 WHERE seq_id = %s",
                (transcript.md5sum + transcript.aa_seq[0:4] + transcript.aa_seq[-4:],))
            if rows_selected_count:
                summary_rows = db.fetchall()
                # Convert the summary_rows Tuple of Dicts returned by SQL to a list of dicts
                modbase2020_summary_rows = [summary_row for summary_row in summary_rows]
            assert rows_selected_count == len(modbase2020_summary_rows)
        return modbase2020_summary_rows if modbase2020_summary_rows else []

    def get_info(self,model_id_database_id_tuple: Tuple[str,str]) -> Dict[str, str]:
        """Return the entire row of information from the summary file for a specific 
        (model_id,database_id) Tuple"""

        with PDBMapSQLdb() as db:
            db.activate_dict_cursor()
            # OK - so we only use the model_id part of the tuple :)
            rows_selected =  db.execute("SELECT * FROM Modbase2020 WHERE model_id = %s",(model_id_database_id_tuple[0],))
            assert rows_selected == 1
            fetched_data = db.fetchone()
            return fetched_data

    def get_coord_file(self,modelid: str) -> str:
        return os.path.join(self.modbase2020_dir,modelid+".pdb.xz")



