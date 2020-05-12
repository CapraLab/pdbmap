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


from lib.PDBMapSQLdb import PDBMapSQLdb

class PDBMapModbase2016():
    def __init__(self,config_dict):
        assert 'modbase2016_dir' in config_dict,"The key modbase2016_dir is required in the passed dictionar - but not found"
        self.modbase2016_dir = config_dict['modbase2016_dir']

    def ensp2modelids(self,ensp_id: str) -> List[str]:
        """Return a list of modelids for an ENSP id"""
        ensp_model_matches = []
        with PDBMapSQLdb() as db:
            rows_selected =  db.execute("SELECT database_id FROM Modbase2016 WHERE database_id LIKE %(ensp_id)s",{'ensp_id': ensp_id+'%%'})
            if rows_selected:
                row_tuples = db.fetchall()
                ensp_model_matches = [row[0] for row in row_tuples]
            assert rows_selected == len(ensp_model_matches)
        return ensp_model_matches if ensp_model_matches else []

    def get_info(self,modelid: str) -> Dict[str, str]:
        with PDBMapSQLdb() as db:
            db.activate_dict_cursor()
            rows_selected =  db.execute("SELECT * FROM Modbase2016 WHERE database_id = %s",(modelid,))
            assert rows_selected == 1
            fetched_data = db.fetchone()
            return fetched_data

    def get_coord_file(self,modelid: str) -> str:
        return os.path.join(self.modbase2016_dir,modelid+".pdb.gz")



