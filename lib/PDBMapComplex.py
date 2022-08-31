#!/usr/bin/env python
# Project        : PDBMap
# Filename       : PDBMapComplex
# Author         : Chris Moth
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : chris.moth@vanderbilt.edu
# Date           : 2022-08-08
# Description    : Manage loading of multichain complexes, alignments to transcripts, variant hooks
#                  Excerpted/factored out of - pathprox3.py

import os
import re
import sys
import gzip
import lzma
import string
import json
import pandas as pd
import numpy as np
import subprocess as sp
import warnings
import platform

from typing import Dict
from typing import List
from typing import Tuple

from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

import warnings

from lib import PDBMapTranscriptBase
from lib import PDBMapTranscriptUniprot
from lib import PDBMapTranscriptFasta
from lib import PDBMapTranscriptEnsembl
from lib import PDBMapAlignment
from lib import PDBMapAlphaFold
from lib import PDBMapSwiss
from lib import PDBMapModbase2020
from lib import PDBMapProtein
from lib import PDBMapGlobals
from lib import PDB36Parser

import logging

from pdbmap.lib.PDBMapAlignment import sifts_best_unps

LOGGER = logging.getLogger(__name__)


class PDBMapComplex:
    """
    Load a structure from the filesystem and adorn it with:
        uniprot transcripts
        ensembl transcripts
        alignments of chain to the AA seqences
        a renumbered structure, based on the alignments
        placed spheres variants of interest
    """

    def __init__(self,
                 structure_type: str,
                 structure_id: str,
                 single_chain: str = None,
                 user_chain_to_transcript: Dict[str, PDBMapTranscriptBase] = None,
                 info_func=(lambda info: info),  # By default, asserts just echo back the passed string
                 usermodel_filename: str = None
                 ):

        """
        Load a structure and note previously chain2transcript
        @param structure_type:  One of 'biounit', 'pdb', 'alphafold', 'swiss', 'modbase', 'usermodel'
            'biounit': A 4 character pdb id must be given for structure id
            'pdb':     A 4 character pdb id must be given for structure id
            'swiss':   A swissmodel ID must be given for structure id

        @param structure_id:    4 char pdb id, swiss model id, alpha fold id, modbase ID

        @param single_chain: If provided, all other chains will be stripped away from the Complex

        @param user_chain_to_transcript: A mapping of chain IDs to Transcript as provided by a user
                                         Example: The -chainBunp=ABC123 syntax of the pathprox command line

        @param info_func:   Function that returns it's single string argument, after having a chance to log it.

        @param usermodel_filename only neededd for structure_type === 'usermodel'
        """

        self.user_chain_to_transcript = user_chain_to_transcript
        self.chain_to_transcript = {}  # mapping of PDB chain letters to transcripts (usually uniprot)
        self.transcript_to_chains = {}  # A reverse dictionary where a uniprot transcript ID maps to a list of chains
        self.chain_to_alignment = {}
        self._unp_to_ENST_transcripts = {}  # For the specific uniprot IDs in the complex, match to list of ENST
        self.ensembl_transcript_to_cosmis: Dict[str, pd.DataFrame] = {}
        self.is_biounit = False
        self.structure_type = structure_type
        self.structure_id = structure_id
        self.single_chain = single_chain
        self.structure: Structure = None
        self.renumbered_structure: Structure = None
        self.source_coordinates_filename = ""
        self.info_func = info_func
        self.usermodel_filename = usermodel_filename

        self._label = self.structure_id  # Typicaly a complex user will override this at a later point

        if self.structure_type == 'biounit' or self.structure_id == 'pdb':
            self._load_pdb()
        elif self.structure_type == 'swiss':
            self._load_swiss()
        elif self.structure_type == 'alphafold':
            self._load_alphafold()
        elif self.structure_type == 'modbase':
            self._load_modbase()
        elif self.structure_type == 'usermodel':
            self._load_usermodel()
        else:
            PDBMapGlobals.exit("structure_type must be one of biounit, pdb, alphafold, swiss, modbase, usermodel")

        # DROP Models after the first one.  Typically this only is in play for NMR ensembles
        # Preprocess structural properties
        # Pathprox cannot use models beyond the 0th one - so drop those right away
        _model_list = list(self.structure.get_models())
        # Drop all the models after the first one.
        if len(_model_list) > 1:
            for _model in _model_list[1:]:
                self.structure.detach_child(_model.id)
            LOGGER.info("Removed %d extraneous models from %s", len(_model_list) - 1,
                        self.structure.id)

        # DROP chains from 0th model If user only
        # interested in a single chain, trim the structure now.
        if single_chain:
            if single_chain not in self.structure[0]:
                msg = "You requested restriction to chain %s.  But, this is not in structure %s" % (
                    single_chain, self.structure.id)
                LOGGER.critical(msg)
                PDBMapGlobals.exit(msg)
            for chain in list(self.structure[0].get_chains()):
                if chain.id != single_chain:
                    LOGGER.info('Dropping chain %s from structure %s (args.chain=%s)',
                                chain.id,
                                self.structure.id,
                                self.single_chain)
                    self.structure[0].detach_child(chain.id)

        if self.structure_type == 'biounit' or self.structure_id == 'pdb':
            self._align_pdb_chains_to_best_uniprot_transcripts_with_sifts()
        else:
            self._align_model_chains_per_user_chain_to_transcript()

        # In a multimer, multiple chains often associate with the same chain.
        # Create that reverse structure lookup structure from transcript ID to list of chain IDs
        self.transcript_to_chains = {}
        for chain_letter in self.chain_to_transcript:
            chains = self.transcript_to_chains.get(self.chain_to_transcript[chain_letter].id, [])
            chains.append(chain_letter)
            self.transcript_to_chains[self.chain_to_transcript[chain_letter].id] = chains

        assert self.structure is not None, "%s %s was not loaded.  Likely due to software bug.  See logs" % (
            self.structure_type, self.structure_id)

    @staticmethod
    def _open_for_extension(file_name: str):
        _, _coordinates_filename_extension = os.path.splitext(file_name)
        if _coordinates_filename_extension == '.gz':
            _open_function, _open_mode = (gzip.open, 'rt')
        elif _coordinates_filename_extension == '.xz':
            _open_function, _open_mode = (lzma.open, 'rt')
        else:
            _open_function, _open_mode = (open, 'r')

        return _open_function(file_name, _open_mode)

    def _load_structure_and_mmcif_dict(self, is_cif_format: bool) -> None:
        # With self.source_coordinates_filename initialized beforehand, call this routine to load
        # a structure file in pdb or cif format.  On return, self.structure will be
        # populated with a Biopython structure.  If the source file is cif format
        # we will also init self.mmcif_dict

        try_parser_36 = False
        self.mmcif_dict = {}
        self.structure = None

        with self._open_for_extension(self.source_coordinates_filename) as fin:
            warnings.filterwarnings('ignore', category=PDBConstructionWarning)
            # source_coordinates_filename_pieces = self.source_coordinates_filename.split('.')
            if is_cif_format:
                _mmcif_parser = MMCIFParser(QUIET=True)
                self.structure = _mmcif_parser.get_structure(self.structure_id, fin)
                self.mmcif_dict = _mmcif_parser._mmcif_dict
            else:
                try:
                    self.structure = PDBParser().get_structure(self.structure_id, fin)
                except ValueError:
                    try_parser_36 = True  # We will make a last ditch effort to read this because of alpha in int columns

            warnings.resetwarnings()

        if try_parser_36:
            LOGGER.critical("Attempting Hybrid36 PDB Parser for %s" % self.source_coordinates_filename)
            with self._open_for_extension(self.source_coordinates_filename) as fin:
                warnings.filterwarnings('ignore', category=PDBConstructionWarning)
                self.structure = PDB36Parser().get_structure(self.structure_id, fin)
                warnings.resetwarnings()

    def _load_pdb(self) -> None:
        """From the pdb ID (previously set as self.structure_id,
           load self.structure"""
        assert self.structure_type == 'pdb' or self.structure_type == 'biounit'
        self.is_biounit = False
        self.structure = None
        self.mmcif_dict = {}
        # Ultimately, we still need to build the biounit file from the raw deposition
        # outselves.  But I've not written that yet.
        if self.structure_type == 'biounit':
            # Let's hope to load the biounit from mmCIF file
            is_cif_format = True
            self.source_coordinates_filename = os.path.join(
                PDBMapGlobals.config['pdb_dir'],
                "biounit",
                "mmCIF",
                "divided",
                self.structure_id.lower()[1:3],  # middle 2 letters of id are a directory layer
                "%s-assembly1.cif.gz" % self.structure_id.lower())
            # Often there is no cif biounit, so we will see if there is a .pdb1.gz biounit next...
            if not os.path.exists(self.source_coordinates_filename):
                is_cif_format = False
                self.source_coordinates_filename = os.path.join(
                    PDBMapGlobals.config['pdb_dir'],
                    "biounit",
                    "PDB",
                    "divided",
                    self.structure_id.lower()[1:3],
                    "%s.pdb1.gz" % self.structure_id.lower())
            if os.path.exists(self.source_coordinates_filename):
                self._load_structure_and_mmcif_dict(is_cif_format)  # init self.structure= and perhaps self.mmcif_dct
                self.is_biounit = True  # Only case where this is true
            else:
                LOGGER.info(
                    "The biounit file %s was not found.  Attempting normal .pdb read" % os.path.basename(
                        self.source_coordinates_filename))

        # Either structure_type is NOT 'biounit' OR biounit file was not found
        # So we next attempt to load a deposited PDB
        if not self.structure:
            # In July 2022, there are NO structures that are .pdb but NOT .cif
            # We shall load from the .cif file
            self.source_coordinates_filename = os.path.join(
                PDBMapGlobals.config['pdb_dir'],
                "structures",
                "divided",
                "mmCIF",
                self.structure_id.lower()[1:3],
                "%s.cif.gz" % self.structure_id.lower())

            self._load_structure_and_mmcif_dict(is_cif_format=True)

    def _align_model_chains_per_user_chain_to_transcript(self):
        """
        Given the presents of self.user_chain_to_transcript
        attempt to _align_ each PDB chain to this transcript using PDBMapAlignment

        @return None:
        """
        self.chain_to_transcript = {}

        assert len(self.user_chain_to_transcript) >= 1, \
            "At least one chain to transcript mapping required on command line to use a model structure"

        # Below we'll assign the user-transcript to any other non-HETATM
        # chains in the swiss model.
        assigned_chain_id = next(iter(self.user_chain_to_transcript))

        # Swiss _could_ be a homo-oliomer - and we want to capture that!
        # by pointing the additional chain IDs at the same alignment
        for chain in self.structure[0]:
            if chain.id not in self.chain_to_transcript:
                first_residue = next(iter(self.structure[0][chain.id]))
                # Swiss model seems to put some HETATMs in some of the chains.  DO NOT align to those!
                # Recall the Biopython residue ID is a triple, that starts with a ' ' in non-Hetero cases
                if first_residue.id[0] == ' ':
                    self.chain_to_transcript[chain.id] = self.user_chain_to_transcript[assigned_chain_id]

        # Now that chain_to_transcript is initialized, create alignment dictionaries to connect model residues
        # to transcript positions, and vice-versa
        for chain_letter in self.chain_to_transcript:
            # For swiss models, the '_' chain is always for HETATM or RNA/DNA ligand type entries - SKIP!
            if self.structure_type == 'swiss' and chain_letter == '_':
                continue
            alignment = PDBMapAlignment()
            (success, message) = alignment.align_trivial(
                self.chain_to_transcript[chain_letter],
                self.structure,
                chain_id=chain_letter)

            if success:
                LOGGER.info("%s to chain %s Trivial align successful for residues in [%d,%d]\n%s",
                            self.chain_to_transcript[chain_letter].id,
                            chain_letter,
                            next(iter(alignment.seq_to_resid)),  # First alignment
                            next(reversed(alignment.seq_to_resid)),  # Last alignment
                            alignment.aln_str)
            elif self.structure_type == 'usermodel':  # Hmm...  Apparently we are OK with a user-model doing this call
                (success, message) = alignment.align_biopython(self.chain_to_transcript[chain_letter], self.structure,
                                                               chain_id=chain_letter)
                if success:
                    LOGGER.warning("TRIVIAL ALIGNMENT OF USER MODEL FAILED\n" +
                                   "%s to chain %s Biopython successful for residues in [%d,%d]\n%s",
                                   self.chain_to_transcript[chain_letter].id,
                                   chain_letter,
                                   next(iter(alignment.seq_to_resid)),
                                   next(reversed(alignment.seq_to_resid)),
                                   alignment.aln_str)

            assert success, self.info_func("Unable to align %s %s chain %s to %s" % (
                self.structure_type,
                self.structure_id,
                chain_letter,
                self.chain_to_transcript[chain_letter].id))
            self.chain_to_alignment[chain_letter] = alignment

    def _align_pdb_chains_to_best_uniprot_transcripts_with_sifts(self):
        assert self.structure_type == 'pdb' or self.structure_type == 'biounit'
        self.chain_to_alignment = {}
        sifts_chain_to_best_unp = sifts_best_unps(self.structure)
        for chain_letter in sifts_chain_to_best_unp:
            if chain_letter not in self.structure[0]:
                continue  # It is A-OK to have a biounit which lacks chains from sifts
            if chain_letter in self.chain_to_transcript:
                best_unp_transcript = self.chain_to_transcript[
                    chain_letter]  # Use the transcript inited from command line
            else:
                sifts_best_uniprot_id = sifts_chain_to_best_unp[chain_letter]
                ensemble_transcript_ids = []
                # Before we get too excited about an isoform specific ID
                # If there is no cross-ref'd ENST transcript from idmapping, then revert to canonical
                if PDBMapProtein.unp2uniparc(sifts_best_uniprot_id) is None:
                    LOGGER.warning(
                        "The best uniprot ID for chain %s is %s, but it has not been reviewed by uniprot.  Considering canonical",
                        chain_letter,
                        sifts_best_uniprot_id)
                else:  # The 'best' uniprot ID is Swiss-curated.  Does it have ensembl_transcript_ids in idmapping file?
                    ensemble_transcript_ids = PDBMapProtein.unp2enst(sifts_best_uniprot_id)
                    if not ensemble_transcript_ids:
                        LOGGER.critical(
                            "Unp %s is best for chain %s. HOWEVER it lacks ENST*.. so reverting to canonical",
                            sifts_chain_to_best_unp[chain_letter], chain_letter)
                if not ensemble_transcript_ids:
                    # To get the canonical uniprot ID, we strip the '-' off any uniprot ID we have in hand
                    # And we ask for the "best" unp.  That will give us n A123456-nn format identifier if indeed
                    # the dashed form exists for this uniprot ID
                    canonical_unp = PDBMapProtein.best_unp(sifts_chain_to_best_unp[chain_letter].split('-')[0])
                    # If this unp is NOT in our curated/reviewed set from uniprot, we need to skip this chain
                    if PDBMapProtein.unp2uniparc(canonical_unp) is None:
                        LOGGER.warning(
                            "Chain %s has canonical unp %s, but it has not been reviewed by uniprot.  Skipping" % (
                                chain_letter, canonical_unp))
                        continue
                    sifts_chain_to_best_unp[chain_letter] = canonical_unp
                    LOGGER.critical("For chain %s, unp now=%s", chain_letter, sifts_chain_to_best_unp[chain_letter])

                best_unp_transcript = PDBMapTranscriptUniprot(sifts_chain_to_best_unp[chain_letter])

            alignment = PDBMapAlignment()
            (success, message) = alignment.align_trivial(best_unp_transcript, self.structure, chain_id=chain_letter)
            if success:
                LOGGER.info("%s to chain %s Trivial align successful for residues in [%d,%d]\n%s" % (
                    best_unp_transcript.id, chain_letter, next(iter(alignment.seq_to_resid)),
                    next(reversed(alignment.seq_to_resid)), alignment.aln_str))
            else:
                LOGGER.info("Unable to trivially align transcript %s to %s.%s" % (
                    best_unp_transcript.id, self.structure_id, chain_letter))
                is_canonical = PDBMapProtein.isCanonicalByUniparc(best_unp_transcript.id)
                if is_canonical:
                    (success, message) = alignment.align_sifts_canonical(best_unp_transcript,
                                                                         self.structure,
                                                                         chain_id=chain_letter)
                if not success:
                    # Then we have to attempt a PDB alignment using the sequence info embedded in the .cif file

                    # If we started from a pdb or an incomplete biounit dict, then we fetch the full mmCIF file here
                    if self.mmcif_dict and PDBMapAlignment.seq_resid_required_keys_present(self.mmcif_dict):
                        _mmcif_dict = self.mmcif_dict
                    else:
                        mmcif_structure_filename = os.path.join(PDBMapGlobals.config['pdb_dir'], 'structures',
                                                                'divided',
                                                                'mmCIF',
                                                                self.structure_id.lower()[1:3],
                                                                "%s.cif.gz" % self.structure_id)
                        if os.path.exists(mmcif_structure_filename):
                            if self.structure_type == 'biounit':
                                LOGGER.info("Biounit mmcif files lack keys necessary for alignment")
                            LOGGER.info("(re)Loading to get alignment keys %s", mmcif_structure_filename)
                            mmCIF_parser = MMCIFParser(QUIET=True)
                            with gzip.open(mmcif_structure_filename, 'rt') as structure_fin:
                                _temp_mmcif_structure = mmCIF_parser.get_structure(self.structure_id.lower(),
                                                                                   structure_fin)
                                _mmcif_dict = mmCIF_parser._mmcif_dict
                            pdb_seq_residue_id_xref = PDBMapAlignment.create_pdb_seq_resid_xref(_mmcif_dict)

                            (success, message) = alignment.align_sifts_isoform_specific(best_unp_transcript,
                                                                                        _temp_mmcif_structure,
                                                                                        pdb_seq_residue_id_xref,
                                                                                        chain_id=chain_letter,
                                                                                        is_canonical=is_canonical)
                    if not success:
                        LOGGER.warning("Unable to align with sifts: %s", message)
                        LOGGER.warning("Attempting to align chain %s with biopython call", chain_letter)
                        success, message = alignment.align_biopython(best_unp_transcript, self.structure, chain_letter)
                        if not success:
                            LOGGER.critical("Also Unable to align with biopython: %s", message)

            assert success, message
            self.chain_to_alignment[chain_letter] = alignment
            self.chain_to_transcript[chain_letter] = best_unp_transcript

    def structure_renumber_per_alignments_and_lose_disorder(self) -> None:
        """
        Create self.renumbered_structure and copy in renumbered resides in chains (and models)
        of the same name as self.structure.
        """

        # Create the empty renumbered_structure to receive the
        # renumbered residues
        self.renumbered_structure = Structure(self.structure.id)

        # Typically, we will have removed all models beyond self.structure[0] at this point.
        # But we might as wlel be thorough.
        for model in self.structure:
            renumbered_model = Model(model.id)
            for chain in model:
                if chain.id not in self.chain_to_alignment:
                    # Only howl about unaligned chains first time through
                    if model.id == 0:
                        LOGGER.warning("*** %s.%s has not been aligned to a human uniprot transcript.",
                                       self.structure.id, chain.id)

                    # This is ugly - but we are just pointing at the old chain
                    # It solves the problem of having disordered residues that are not brought over by .copy()
                    renumbered_model.add(chain)
                else:
                    alignment = self.chain_to_alignment[chain.id]
                    renumbered_chain = Chain(chain.id)
                    for residue in chain:
                        if residue.id in alignment.resid_to_seq:
                            renumbered_residue = residue.copy()
                            renumbered_residue.id = (' ', alignment.resid_to_seq[residue.id], ' ')
                            # The above copy routines fail to copy disorder properly
                            # - so just wipe out all notion of disorder
                            for atom in renumbered_residue:
                                atom.disordered_flag = 0
                            renumbered_residue.disordered = 0
                            renumbered_chain.add(renumbered_residue)
                    renumbered_model.add(renumbered_chain)

            self.renumbered_structure.add(renumbered_model)

    @staticmethod
    def uniprot_and_ensembl_close_enough(
            uniprot_transcript: PDBMapTranscriptUniprot,
            ensembl_transcript: PDBMapTranscriptEnsembl) -> bool:
        """
        Uniprot and Ensembl transcripts can drift.  That's OK if only in a few positions.
        Here we isolate the "are two transcripts close enough?" logic so that we can share the code easily
        between initial pipeline planning and later pathprox transcript acquisition
        """
        if not ensembl_transcript.aa_seq:
            LOGGER.warning("Ensembl transcript %s mapped to %s has no associated aa_seq.  Skipping",
                           ensembl_transcript.id,
                           uniprot_transcript.id)
            return False

        differences, pct_different = PDBMapTranscriptBase.analyze_transcript_differences(
            uniprot_transcript, ensembl_transcript)

        if differences == 0:  # Awesome - the usual case where uniprot and ENST match
            LOGGER.info("Ensembl transcript %s has same aa_seq as transcript %s", ensembl_transcript.id,
                        uniprot_transcript.id)
        elif differences == 1:  # Let's not kill pathprox if only one variant between uniprot and ENST
            LOGGER.warning(
                "Transcripts vary in one position: %s" % PDBMapTranscriptBase.describe_transcript_differences(
                    uniprot_transcript, ensembl_transcript))
        elif pct_different < 0.01:  # Similarly let Pathprox continue if we have a few variants off - but less than 1% of the sequence
            LOGGER.warning(
                "Transcripts vary in multiple positions: %s" % PDBMapTranscriptBase.describe_transcript_differences(
                    uniprot_transcript, ensembl_transcript))
        else:
            LOGGER.warning("%s and %s AA sequences differ markedly:\n%s" % (
                uniprot_transcript.id, ensembl_transcript.id,
                PDBMapTranscriptBase.describe_transcript_differences(uniprot_transcript, ensembl_transcript)))
            return False
        return True

    def uniprot_to_ENST_transcripts(self, uniprot_transcript: PDBMapTranscriptUniprot) -> List[PDBMapTranscriptEnsembl]:
        """

        @param uniprot_transcript: The Uniprot Transcript for which matching ENST Ensembl Transcripts are desired
        @return: A filtered list of ensembl transcripts, which is taken from the Uniprot IDMapping first, and then
                filtered to ensure that the Ensembl transcript AA lengths match the uniprot length, and variantions
                are minimal
        """
        # If we already fetched this before, then return that result
        if not self._unp_to_ENST_transcripts:
            self._unp_to_ENST_transcripts = {}
        if uniprot_transcript.id in self._unp_to_ENST_transcripts:
            return self._unp_to_ENST_transcripts[uniprot_transcript.id]

        # Create an empty list
        self._unp_to_ENST_transcripts[uniprot_transcript.id] = []

        # We need to load the ENST transcripts for this uniprot ID
        # And we need to filter to make sure that we match AA sequence
        unfiltered_ensemble_transcript_ids = PDBMapProtein.unp2enst(uniprot_transcript.id)

        # Let's now iterate through all the possible ENST associated with the unp
        # and bin them for exact matches
        # Let's STOP killing pathprox in the event that ONE of the transcripts x-ref'd in IDMapping
        # has a big problem
        _ENST_identical_transcripts = []

        for ensembl_transcript_id in unfiltered_ensemble_transcript_ids:
            ensembl_transcript = PDBMapTranscriptEnsembl(ensembl_transcript_id)
            if not PDBMapComplex.uniprot_and_ensembl_close_enough(uniprot_transcript, ensembl_transcript):
                continue

            # To arrive here we have reasonably matching AA sequences between the uniprot and ENST ensembl transcript
            # So, add the ENST transcript to the list representing the chain
            self._unp_to_ENST_transcripts[uniprot_transcript.id].append(ensembl_transcript)
            # if ensembl_transcript not in ENST_transcripts:
            #    ENST_transcripts.append(ensembl_transcript)

        if not self._unp_to_ENST_transcripts[uniprot_transcript.id]:
            LOGGER.critical("No matching ENSMBL transcripts were found for transcript %s",
                            str(uniprot_transcript.id))

        return self._unp_to_ENST_transcripts[uniprot_transcript.id]

    def _load_swiss(self) -> None:
        """
        Load a swiss model given prior setting of self.structure_id
        @return:
        """
        assert self.structure_type == 'swiss'
        PDBMapSwiss.load_swiss_INDEX_JSON(PDBMapGlobals.config['swiss_dir'], PDBMapGlobals.config['swiss_summary'])
        self.source_coordinates_filename = PDBMapSwiss.get_coord_file(self.structure_id)
        self._load_structure_and_mmcif_dict(is_cif_format=False)

    def _load_modbase(self):
        assert self.structure_type == 'modbase'
        modbase20 = PDBMapModbase2020(PDBMapGlobals.config)
        self.source_coordinates_filename = modbase20.get_coord_file(self.structure_id)
        if os.path.exists(self.source_coordinates_filename):
            self._load_structure_and_mmcif_dict(is_cif_format=False)
        else:
            msg = "Modbase model %s not found in %s" % (self.structure_id, self.source_coordinates_filename)
            LOGGER.critical(msg)
            PDBMapGlobals.exit(msg)

    def _load_alphafold(self):
        assert self.structure_type == 'alphafold'
        alpha_fold = PDBMapAlphaFold(PDBMapGlobals.config)
        self.source_coordinates_filename = alpha_fold.get_coord_filename(self.structure_id)
        if os.path.exists(self.source_coordinates_filename):
            self._load_structure_and_mmcif_dict(is_cif_format=True)
            self.alpha_fold_local_metrics = alpha_fold.pLDDTs(self.mmcif_dict)
            # If we are dealing with not a -F1-... but a -F2- or higher, then we need to reumber the 1..N by the window*200
            self.structure = alpha_fold.renumber_windowed_model(self.structure, self.mmcif_dict)
        else:
            msg = "Alpha fold model %s not found in %s" % (self.structure_id, self.source_coordinates_filename)
            LOGGER.critical(msg)
            PDBMapGlobals.exit(msg)

    def _load_usermodel(self):
        assert self.structure_type == 'usermodel'

        self.source_coordinates_filename = self.usermodel_filename
        is_cif_format = (self.source_coordinates_filename.find("cif") != -1)

        self._load_structure_and_mmcif_dict(is_cif_format)

    @staticmethod
    def parse_user_chain_to_transcript(args_remaining: str) -> Dict[str, str]:
        """
        Sometimes a user will know precisely which transcript(s) they would like to associate with
        specific protein chains.  This routine gathers those user specifications from the command line,
        after other aguments are gathered, with formats like --chainBunp=P12435 or ==chainC1enst=ENST00000001243567
        See pathprox3.py for examples of how this fits inwith argument parsing generaly.

        @param args_remaining:   # Left over from the argparser of main entry point, after known standard command line
        @return: use_chain_to_transcript dictionary
        """
        user_chain_to_transcript: Dict[str, str] = {}

        chain_x_unp_re = re.compile('^--chain(.{1,3})unp=(.*)$')
        chain_x_enst_re = re.compile('^--chain(.{1,3})enst=(.*)$')
        chain_x_fasta_re = re.compile('^--chain(.{1,3})fasta=(.*)$')
        transcript = None
        for arg in args_remaining:
            re_match = chain_x_unp_re.match(arg)
            chain_letter = None

            if re_match:  # A uniprot ID transcript ID was explicitly assigned to the chain
                chain_letter = re_match.group(1)
                uniprot_id = re_match.group(2)
                transcript = PDBMapTranscriptUniprot(uniprot_id)
                # if chain_letter in chain_to_transcript and transcript.id == chain_to_transcript[chain_letter].id:
                #    LOGGER.warning("argument %s is redundant.  Pathprox automatically assigned/aligned")
                LOGGER.info("Successful load of Uniprot transcript %s", uniprot_id)
            else:
                re_match = chain_x_enst_re.match(arg)
                if re_match:  # An ENSEMBL transcript ID was explicitly assigned to the chain
                    chain_letter = re_match.group(1)
                    ensembl_transcript_id = re_match.group(2)
                    transcript = PDBMapTranscriptEnsembl(ensembl_transcript_id)
                    LOGGER.info("Successful load of Ensembl transcript %s", ensembl_transcript_id)
                else:
                    re_match = chain_x_fasta_re.match(arg)
                    if re_match:  # A fasta amino acid string was explicitly assigned to the chain
                        chain_letter = re_match.group(1)
                        fasta = re_match.group(2)
                        transcript = PDBMapTranscriptFasta(fasta)
                        LOGGER.info("Successful load of fasta transcript %s" % fasta)
                    else:
                        exit_message: str = "Command line argument %s is invalid" % arg
                        LOGGER.critical(exit_message)
                        PDBMapGlobals.exit(exit_message)

            if transcript:
                user_chain_to_transcript[chain_letter] = transcript
        return user_chain_to_transcript

    def assign_chain_letters_where_blank(self) -> None:
        """
        Iterate over chains in the complex, and assign letter codes when
        blank chain IDs are found.  A cleanup for old modbase, mainly.
        @return None:
        """

        letters = list(string.ascii_uppercase + string.ascii_lowercase)
        next_chain_id_subscript = 0  # We start with 'A'
        for chain in list(self.structure.get_chains()):
            if chain.id in ('', ' '):
                # No way that we have 52 unlabelled chains...
                # But don't use a chian ID already in the structure.
                while letters[next_chain_id_subscript] in self.structure[0]:
                    next_chain_id_subscript += 1

                next_chain_id = letters[next_chain_id_subscript]
                LOGGER.info('Renaming blank/missing chain ID to %s', next_chain_id)
                chain.id = next_chain_id
                next_chain_id_subscript += 1

    @staticmethod
    def residue_center_of_mass(residue: Residue) -> np.ndarray:
        """
        Given a BioPython Residue, iterate over all the heavy (non-hydrogen) atom locations to
        calculate a center of mass (without regard for varying mass of Carbon/Nitrogen/Oxygen/Sulfur
        @param residue:
        @return: The x/y/z mean of the heavy atom coordinates inside the residue
        """
        if residue is None:
            return np.array([np.NaN, np.Nan, np.Nan])

        # 2022 August 5.  I have reworked this to avoid counting hydrogens  (nmr/high res PDBs)
        # and I also am averaging coordinates of heavy atoms which are disorded instead of simply recounting them
        residue_atoms_xyzs = []
        for atom in residue.get_list():
            if atom.element == 'H':  # Skip Hydrogen
                continue

            if atom.is_disordered():
                disordered_atom_xyzs = [np.array(disordered_atom.get_coord())
                                        for disordered_atom in atom.disordered_get_list()]
                atom_mean_xyz = np.array(disordered_atom_xyzs).mean(axis=0)
            else:
                atom_mean_xyz = np.array(atom.get_coord())

            residue_atoms_xyzs.append(atom_mean_xyz)
        return np.array(residue_atoms_xyzs).mean(axis=0)

    # Construct a dataframe from the coordinate file
    def unp_pos_if_aligned(self, chain_id: str, residue_id: Tuple) -> int:
        """
        Given a Chain ID, and residue ID in the chain, return the 1..n amino acide number of the transcript
        that was aligned to the chain and the residue

        It is often the case (with PDB files) that resolved amino acids are not aligned to transcript postions
        In those cases, return 0

        @param chain_id:
        @param residue_id:
        @return:
        """

        # If chain 'X' is not aligned to a transcript at all, then return 0
        if chain_id not in self.chain_to_alignment:
            return 0

        _alignment = self.chain_to_alignment[chain_id]

        # If the Biopython residue ID triple was never aligned to a transcript
        # position, return 0
        if residue_id not in _alignment.resid_to_seq:
            return 0

        # Success, we can say that a particular PDB residue is aligned to a transcript position
        return int(_alignment.resid_to_seq[residue_id])

    def create_residue_centroids_dataframe_from_aligned_chains(self) -> pd.DataFrame:
        """
        In one iteration through all the residues the entire loaded comple, we create a new dataframe consisting of
        6 columns
        1: 'chain'. The chain ID of the residue, found via residue.get_parent().id
        2: 'resid'.  The residue ID Tuple of the residue
        3: 'unp_pos' The amino acid position in the uniprot transcript, determined by prior alignment of chain
        4-6: the x/y/z center of mass of the residue
        @return:
        """
        # The "for r..." list comprehension includes only non-hetero residues with (blank,#,insert code) i.e. (not r.id(0).strip())
        # and only then of they are aligned to transcript positions i.e. unp_pos_if_aligned... != None
        _dataframe_rows = []
        for chain in self.structure[0]:
            if chain.id not in self.chain_to_alignment:
                continue
            for residue in chain:
                if len(residue.id[0].strip()) == 0:  # Ignore Hetero residues (len != 0)
                    _uniprot_transcript_pos = self.unp_pos_if_aligned(chain.id, residue.id)
                    if _uniprot_transcript_pos:
                        _dataframe_rows.append(
                            [chain.id,  # 1
                             residue.id,  # 2
                             _uniprot_transcript_pos  # 3
                             ] + list(PDBMapComplex.residue_center_of_mass(residue))  # 456
                        )

        return pd.DataFrame(_dataframe_rows, columns=["chain", "resid", "unp_pos", "x", "y", "z"])

    @staticmethod
    def _load_one_cosmis_file(cosmis_filename: str) -> Tuple[pd.DataFrame, float, float, float]:
        """

        @param cosmis_filename:
        @return: Dataframe of the COSMIS scores from the file,  alpha_parameter, average, stddev
        """
        # you have to see thees files - lots of COMMENTS at top, blank lines
        headers_re = re.compile('^#POS *SEQ *SCORE *QQ-INTERVAL *STD *MSA *DATA *$')
        float_regex = '([-+]?[0-9]*\.?[0-9]*)'
        spaces = ' *'
        int_regex = '([0-9]*)'
        cosmis_data_re = re.compile(
            '^' + spaces + int_regex +  # POS
            spaces + '([A-Z])' +  # SEQ
            spaces + float_regex +  # SCORE
            spaces + '[' + spaces + float_regex +  # QQ-INTERVAL_low
            ',' + spaces + float_regex + ']' +  # QQ-INTERVAL_high
            spaces + float_regex +  # STD
            spaces + int_regex +  # MSA
            '/100[ \n]*$'
        )

        alpha_parameter = None
        average = None
        std = None

        alpha_parameter_re = re.compile('^#The alpha parameter ' + float_regex + '$')
        average_re = re.compile('^#Average = ' + float_regex + '$')
        standard_deviation_re = re.compile('^#Standard Deviation = ' + float_regex + '$')

        cosmis_line_list = []
        LOGGER.info('Attempting load of cosmis scores from %s', cosmis_filename)
        with open(cosmis_filename) as cosmis_f:
            headers_seen = False
            for line in cosmis_f.readlines():
                # The data format is a bit ill-defined.  Let's not start parsing lines until we know
                # we have mostly gone through the headers.
                if not headers_seen:  # this is a comment in Bian's output format
                    if line.startswith('#') and headers_re.match(
                            line):  # Then wcosmis_line_liste've match the #POS SEQ etc next-to=last header
                        headers_seen = True
                else:
                    if alpha_parameter is None:
                        m = alpha_parameter_re.match(line)
                        if m:
                            alpha_parameter = float(m.group(1))

                    else:  # Parse out the components of what s likely s data line
                        m = cosmis_data_re.match(line)
                        if m:
                            cosmis_dict = {'pos': int(m.group(1)),
                                           'seq': str(m.group(2)),
                                           'score': float(m.group(3)),
                                           'qq-interval_low': float(m.group(4)),
                                           'qq-interval_high': float(m.group(5)),
                                           'std': float(m.group(6)),
                                           'msa': int(m.group(7))}
                            cosmis_line_list.append(cosmis_dict)
                        else:
                            m = average_re.match(line)
                            if m:
                                average = float(m.group(1))
                            else:
                                m = standard_deviation_re.match(line)
                                if m:
                                    std = float(m.group(1))

        assert alpha_parameter, "Alpha Parameter was not found in in %s" % cosmis_filename
        assert average is not None, "Average was not found in %s" % cosmis_filename
        assert std, "Standard Deviation was not found in %s" % cosmis_filename
        assert len(cosmis_line_list) == cosmis_line_list[-1]['pos']
        cosmis_df = pd.DataFrame(cosmis_line_list)
        return cosmis_df, alpha_parameter, average, std

    def load_cosmis_scores(self):
        for chain in self.structure[0]:
            if chain.id not in self.chain_to_transcript:
                LOGGER.warning("Chain %s does not represent a Human uniprot ID.  Skipping COSMIS load", chain.id)
                continue
            ensembl_transcripts = self.uniprot_to_ENST_transcripts(self.chain_to_transcript[chain.id])
            for ensembl_transcript in ensembl_transcripts:
                # We kiad Bian Li's normalized cosmis scores (average-0, stddev=1)
                if ensembl_transcript.id not in self.ensembl_transcript_to_cosmis:
                    cosmis_norm_rates_filename = os.path.join(PDBMapGlobals.config['cosmis_dir'],
                                                              "%s_norm_rates.txt" % ensembl_transcript.id)
                    df_cosmis_scores, cosmis_alpha_parameter, cosmis_average, cosmis_stddev = \
                        self._load_one_cosmis_file(cosmis_norm_rates_filename)
                    self.ensembl_transcript_to_cosmis[ensembl_transcript.id] = df_cosmis_scores
                    # df_cosmis_scores, cosmis_alpha_parameter, cosmis_average, cosmis_stddev

                # cosmis_orig_rates_filename = os.path.join(PDBMapGlobals.config['cosmis_dir'],
                # "%s_orig_rates.txt" % ensembl_transcript.id)

    def write_cosmis_scores(self, cosmis_scores_json_filename: str):
        cosmis_scores_output_dict = {}
        for chain in self.structure[0]:
            if chain.id not in self.chain_to_transcript:
                LOGGER.warning("Chain %s does not represent a Human uniprot ID.  Skipping COSMIS load", chain.id)
                continue
            ensembl_transcripts = self.uniprot_to_ENST_transcripts(self.chain_to_transcript[chain.id])
            for ensembl_transcript in ensembl_transcripts:
                df_cosmis_scores = self.ensembl_transcript_to_cosmis[ensembl_transcript.id]
                # Each chain will now point to a dictionary which maps the residue numbers to a dictionary of
                # cosmis values (seq/score/qq-interval_low and _high, std, msa
                cosmis_scores_output_dict[chain.id] = df_cosmis_scores.set_index('pos').T.to_dict('dict')
                break
        # cosmis_scores_json_filename ="cosmis_scores.json"
        with open(cosmis_scores_json_filename, 'w') as f:
            json.dump(cosmis_scores_output_dict, f)
        LOGGER.info("Cosmis scores written to %s", cosmis_scores_json_filename)

    # Supplement with any requested variant datasets from SQL
    # 2019 October - loop over chain_to_transcript
    # convert to lists of ENST
    # load variants in a per-transcript manner

    def chain_to_ENST_transcripts(self, chain_id) -> List[PDBMapTranscriptEnsembl]:
        """Given a chain ID from the structure, return a list of ensembl transcripts which
           can be used to query genomic databases for known variants.

           Requires dict chain_to_transcript to map chain_id to a transcript """

        if chain_id not in self.chain_to_transcript:
            LOGGER.critical(
                "Chain %s is not associated with any transcript ID.  Variants cannot be loaded from SQL for chain %s",
                chain_id, chain_id)
            return []

        transcript = self.chain_to_transcript[chain_id]

        # Uniprot identifiers are needed to query the PDBMapProtein idmapping file
        uniprot_ids = []  # Empty list of uniprot identifiers

        enst_transcripts = []
        if type(transcript) == PDBMapTranscriptUniprot:
            uniprot_ids = [transcript.id]
        elif type(transcript) == PDBMapTranscriptEnsembl:
            enst_transcripts = [transcript]
            uniprot_ids = PDBMapProtein.enst2unp(transcript.id)

        for unp in uniprot_ids:
            uniprot_transcript = PDBMapTranscriptUniprot(unp)
            ensemble_transcript_ids = PDBMapProtein.unp2enst(unp)
            # Let's now iterate through all the possible ENST associated with the unp
            # and bin them for exact matches
            # Let's STOP killing pathprox in the event that ONE of the transcripts x-ref'd in IDMapping
            # has a big problem
            _ENST_identical_transcripts = []

            for ensembl_transcript_id in ensemble_transcript_ids:
                ensembl_transcript = PDBMapTranscriptEnsembl(ensembl_transcript_id)
                if ensembl_transcript not in enst_transcripts:
                    if PDBMapComplex.uniprot_and_ensembl_close_enough(uniprot_transcript, ensembl_transcript):
                        enst_transcripts.append(ensembl_transcript)

        if not enst_transcripts:
            LOGGER.critical("No ENSMBL transcripts were referenced to chain %s transcripts %s",
                            chain_id, str(uniprot_ids))
        return enst_transcripts

    @property
    def label(self) -> str:
        return self._label

    @label.setter
    def label(self, label) -> None:
        self._label = label

    @property
    def _base_pdb_filename(self) -> str:
        return "%s_%s_%s" % (
            self.label, self.structure.id, 1 if self.is_biounit else 0)

    @property
    def renumbered_pdb_filename(self) -> str:
        return "%s_renum.pdb" % self._base_pdb_filename

    @property
    def renumbered_cif_filename(self) -> str:
        return "%s_renum.cif" % self._base_pdb_filename

    @property
    def renumbered_pdb_filenameSS(self) -> str:
        return "%s_renumSS.pdb" % self._base_pdb_filename  # Renumbered + chimera HELIX/SHEET notations

    @property
    def renumbered_cif_filenameSS(self) -> str:
        return "%s_renumSS.cif" % self._base_pdb_filename  # Renumbered + chimera HELIX/SHEET notations

    @property
    def original_pdb_filename(self) -> str:
        return "%s.pdb" % self._base_pdb_filename

    @property
    def original_cif_filename(self) -> str:
        return "%s.cif" % self._base_pdb_filename

    def write_renumbered(self, file_format: str):
        if file_format == 'mmCIF':
            LOGGER.info("Writing renumbered CIF to file %s", self.renumbered_pdb_filename)
            cifio = MMCIFIO()
            cifio.set_structure(self.renumbered_structure)
            cifio.save(self.renumbered_cif_filename)
        elif file_format == 'pdb':
            LOGGER.info("Writing renumbered PDB to file %s", self.renumbered_pdb_filename)
            pdbio = PDBIO()
            pdbio.set_structure(self.renumbered_structure)
            pdbio.save(self.renumbered_pdb_filename)
        elif file_format == 'chimera2':
            # We now have a "clean" renumbered PDB - but for improved ngl visualization, helps to
            # let chimera add some secondary structure annotation
            script_dir = "%s/bin" % os.path.dirname(os.path.abspath(sys.argv[0]))

            script = '"%s/writePdbWithSecondary.py %s %s"' % (
                script_dir, self.renumbered_pdb_filename, self.renumbered_pdb_filenameSS)
            cmd = "TEMP=$PYTHONPATH; unset PYTHONPATH; %s --nogui --silent --script %s; export PYTHONPATH=$TEMP" % (
                PDBMapGlobals.config['chimera_headless'], script)
            # Allow Mac OSX to use the GUI window
            if platform.system() == 'Darwin':
                cmd = "TEMP=$PYTHONPATH; unset PYTHONPATH; chimera --silent --script %s; export PYTHONPATH=$TEMP" % \
                      script
            try:
                LOGGER.info("Running Chimera script: %s", cmd)
                # status  = os.system(cmd)
                completed_process = sp.run(cmd, capture_output=True, shell=True, text=True)

                if completed_process.stdout and len(completed_process.stdout) > 0:
                    LOGGER.info("chimera stdout: %s", completed_process.stdout)
                if completed_process.stderr and len(completed_process.stderr) > 0:
                    LOGGER.info("chimera stderr: %s", str(completed_process.stderr))
                if completed_process.returncode != 0:
                    raise Exception("Chimera process returned non-zero exit status.")
            except Exception as e:
                LOGGER.exception("Chimera failed")
                raise
        else:
            PDBMapGlobals.exit("write_renumbered format=%s unrecognized" % file_format)

    def write_original(self, file_format='mmCIF') -> None:
        assert file_format == 'mmCIF', "Only mmCIF supported for write_original at the moment"
        #
        # Write the original structure to the output directory
        #
        LOGGER.info("Writing original structure to file %s", self.original_cif_filename)
        cifio = MMCIFIO()

        cifio.set_structure(self.structure)
        cifio.save(self.original_cif_filename)
