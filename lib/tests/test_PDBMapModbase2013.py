#!/usr/bin/env python
import pytest
import warnings
import logging
import pprint

from Bio.PDB import *
# from Bio import PDBConstructionWarning
from lib.PDBMapGlobals import PDBMapGlobals
from lib.PDBMapTranscriptEnsembl import PDBMapTranscriptEnsembl
from lib.PDBMapModbase2013 import PDBMapModbase2013

LOGGER = logging.getLogger()
# warnings.simplefilter('ignore', PDBConstructionWarning)

config = PDBMapGlobals.config
config['dbname'] = 'pdbmap_v14'

def test_modbase2013():
    modbase = PDBMapModbase2013(config)
    modelids = modbase.ensp2modelids('ENSP00000343764')
    assert len(modelids) > 140, "Unable to find common modbase ID in database"
    assert type(modelids[0]) == str,"Wrong datatype in returned model Ids"

