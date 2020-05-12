#!/usr/bin/env python
import pytest
import warnings
import logging
import pprint

from Bio.PDB import *
# from Bio import PDBConstructionWarning
from lib.PDBMapGlobals import PDBMapGlobals
from lib.PDBMapTranscriptEnsembl import PDBMapTranscriptEnsembl
from lib.PDBMapModbase2016 import PDBMapModbase2016

LOGGER = logging.getLogger()
# warnings.simplefilter('ignore', PDBConstructionWarning)

config = PDBMapGlobals.config
config['dbname'] = 'pdbmap_v14'

def test_modbase2016():
    modbase = PDBMapModbase2016(config)
    modelids = modbase.ensp2modelids('ENSP00000343764')
    assert len(modelids) > 200, "Unable to find common modbase ID in database"
    assert type(modelids[0]) == str,"Wrong datatype in returned model Ids"

