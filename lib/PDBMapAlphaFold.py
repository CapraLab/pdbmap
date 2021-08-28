#!/usr/bin/env python
# Project        : PDBMap
# Filename       : PDBMapAlphaFold
# Author         : Chris Moth
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : chris.moth@vanderbilt.edu
# Date           : 2021-08-01
# Description    : Manage browsing of downloaded AlphaFold models

import argparse
# import json

import sys,os,re
from typing import List
from typing import Dict
from typing import Tuple

from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

from lib.PDBMapTranscriptBase import PDBMapTranscriptBase

class PDBMapAlphaFold():
    def __init__(self,config_dict):
        assert 'alphafold_dir' in config_dict,"The key alphafold_dir is required in the passed dictionar - but not found"
        self._alphafold_dir = config_dict['alphafold_dir']

    def first_window_modelid(self,canonical_unp: str) -> str:
        return "AF-%s-F1"%canonical_unp.split('-')[0]

    def best_covering_model(self, canonical_unp: str, sequence_len: int, variant_of_interest_pos: int) -> Tuple[int,int,str]:
        # The simpleent case is when the transcript has length 2700 or less, and we thus have only one model to fret about
        if sequence_len <= 2700:
            return (1,sequence_len,self.first_window_modelid(canonical_unp))

        # For longer transcripts, alpha-fold is going to create models in a 1400-residue sliding window of every 200
        # the -F1- model will contain 1-1400, -F2- will contain 201-1600, etc.
        # Figure out _which_ of the 200-long windows the variant of interest lies in
        window_200 = (variant_of_interest_pos-1) // 200
        window_200_last = (sequence_len-1) // 200
        model_window_first = window_200 - 3
        model_window_last = window_200 + 3 
        # Slide model selector to right if we fell 
        # off left edge of the transcript
        if model_window_first < 0:
            model_window_last -= model_window_first
            model_window_first = 0
        elif model_window_last > window_200_last:
            model_window_first -= (model_window_last - window_200_last)
            model_window_last = window_200_last

        assert model_window_last - model_window_first == 6,"Problems with finding best covering alpha model on large transcript for position %d"%variant_of_interest_pos
        model_id = "AF-%s-F%d"%(canonical_unp.split('-')[0],model_window_first + 1)
        model_seq_start = 1 + model_window_first * 200
        model_seq_end = (model_window_last + 1) * 200
        if model_seq_end > sequence_len:
            model_seq_end = sequence_len

  
        return (model_seq_start, model_seq_end, model_id)
            
    def get_coord_filename(self,modelid: str) -> str:
        # The model ID is AF-UNPUNP-model_v1_.cif
        # Hence, this works
        return os.path.join(self._alphafold_dir,'human',modelid[3:5],modelid[5:7],modelid + "-model_v1.cif.gz")

    def renumber_windowed_model(self,structure: Structure, alphafold_mmCIF_dict:Dict ) -> Structure:
        # Grab the Alphafold dictionary entry that descrives the residue range in the structure
        seq_db_align_begin = int(alphafold_mmCIF_dict['_ma_target_ref_db_details.seq_db_align_begin'][0])
        seq_db_align_end = int(alphafold_mmCIF_dict['_ma_target_ref_db_details.seq_db_align_end'][0])

        # start empty
        renumbered_structure = Structure(structure.id)
        for model in structure:
            renumbered_model = Model(model.id)
            for chain in model:
                transcript_residue_number = seq_db_align_begin
                renumbered_chain = Chain(chain.id)
                for residue in chain:
                    renumbered_residue = residue.copy()
                    renumbered_residue.id = (' ',transcript_residue_number,' ')
                    # The above copy routines fail to copy disorder properly - so just wipe out all notion of disorder
                    for atom in renumbered_residue:
                        atom.disordered_flag = 0
                    renumbered_residue.disordered = 0
                    renumbered_chain.add(renumbered_residue)
                    transcript_residue_number += 1

                assert transcript_residue_number == seq_db_align_end + 1
                renumbered_model.add(renumbered_chain)

            renumbered_structure.add(renumbered_model)
        return renumbered_structure

