#!/usr/bin/env python
"""Very specific tests of transcript alignments"""

import pytest
import sys
import warnings
import logging
import pprint

from Bio.PDB import *
from lib.mmCIF_to_biounits import mmCIF_to_biounits

# warnings.simplefilter('ignore', PDBConstructionWarning)

def exercise_mmCIF_to_biounits(pdb_id):
    LOGGER.info("Attempting biounit creation for %s"%pdb_id)
    filename = PDBList().retrieve_pdb_file(pdb_id,file_format='mmCif',pdir='/tmp',overwrite=True)
    LOGGER.info("Successful download of %s to %s",pdb_id,filename)
    sys.exit(1)
    mmcif_parser = MMCIFParser(QUIET=True)
    structure_6DWU = mmcif_parser.get_structure('pdb',filename)
    biounits_list= mmCIF_to_biounits(mmcif_parser._mmcif_dict)
    assert type(biounits_list[0]) == dict
    LOGGER.info("%d biounits returned for %s",len(biounits_list),pdb_id)

    # Can we save the biounit???
    mmcif_io = MMCIFIO()
    mmcif_io.set_dict(biounits_list[0])
    mmcif_io.save('/tmp/%s1.cif'%pdb_id)
    return biounits_list


def test_mmCIF_to_biounits():
    # biounits_list = exercise_mmCIF_to_biounits('6dw1')
    # assert len(biounits) == 2,"Two biounits should have been created - but only got %d"%len(biounits)

    biounits_list = exercise_mmCIF_to_biounits('4rkk')
    assert len(biounits_list) == 1
    biounits_list = exercise_mmCIF_to_biounits('3c70')
    assert len(biounits_list) == 1

FORMAT='%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
logging.basicConfig(format=FORMAT,level=logging.DEBUG)
LOGGER = logging.getLogger()
LOGGER.info("Starting...")
test_mmCIF_to_biounits() 


