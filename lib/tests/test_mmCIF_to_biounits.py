#!/usr/bin/env python
"""Very specific tests of transcript alignments"""

import pytest
import warnings
import logging
import pprint

from Bio.PDB import *
from lib.mmCIF_to_biounits import mmCIF_to_biounits

LOGGER = logging.getLogger()
# warnings.simplefilter('ignore', PDBConstructionWarning)

def test_mmCIF_to_biounits():
    """
    6DWU should return two biounit file

    First, align using the new sifts isoform specific mechanism.
    Second, aling with biopython Needleman-Wunsch
    Third, attempt a terrible alignment
    """

    # Alignment 1: Sifts isoform specific alignment
    filename_6DWU = PDBList().retrieve_pdb_file('6DWU',file_format='mmCif',pdir='/tmp',overwrite=True)
    mmcif_parser = MMCIFParser(QUIET=True)
    structure_6DWU = mmcif_parser.get_structure('6DWU',filename_6DWU)
    # import pdb; pdb.set_trace()
    biounits = mmCIF_to_biounits(mmcif_parser._mmcif_dict)
    assert len(biounits) == 2,"Two biounits should have been created - but only got %d"%len(biounits)


test_mmCIF_to_biounits() 

