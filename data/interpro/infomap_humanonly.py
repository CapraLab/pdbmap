#!/usr/bin/env python3
"""
Chris Moth 2017-Nov 30 and 2022-Sept 09
Extract ONLY human entries from the ~600GB match_complete.xml
The resulting match_human.xml file in the same directory (specify in -c configfile)
is used to help generate domain graphics for isoform-specific uniprot ACs in the pipeline
output reports

The strategy is to iterate through each node of the match_complete.xml,
skip those that are non-human, and then write everything else into match_human.xml
"""

import argparse
import configparser
import os
from lib import PDBMapProtein

# The lxml library fits our need to iterate through the huge match_complete XML file
# https://lxml.de/parsing.html
from lxml import etree

import logging
from logging.handlers import RotatingFileHandler
import gzip
from typing import Generator

logging.basicConfig(format='%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
                    datefmt='%d-%m-%Y:%H:%M:%S', )

LOGGER = logging.getLogger()
LOGGER.setLevel(logging.INFO)

fh = RotatingFileHandler("infomap_humanonly.log", maxBytes=(1048576 * 5), backupCount=7)
formatter = logging.Formatter('%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s')
fh.setFormatter(formatter)
fh.setLevel(logging.INFO)
LOGGER.addHandler(fh)

# Setup the Config File Parser
conf_parser = argparse.ArgumentParser(add_help=False)
conf_parser.add_argument("-c", "--conf_file",
                         help="Specify config file", metavar="FILE", required=True)
args, remaining_argv = conf_parser.parse_known_args()
conf_file = args.conf_file
LOGGER.info('Reading config file: %s', conf_file)
config = configparser.SafeConfigParser()
config.read([conf_file])
defaults = {}
defaults.update(dict(config.items("Genome_PDB_Mapper")))
LOGGER.info('Interpro directory: %s', defaults['interpro_dir'])

xml_input_file = os.path.join(defaults['interpro_dir'], 'match_complete.xml.gz')
xml_output_file = os.path.join(defaults['interpro_dir'], 'match_humanonly.xml')

LOGGER.info("Writing hard-coded header to output file: %s", xml_output_file)
xml_output = open(xml_output_file, "w+")
xml_output.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
xml_output.write('<!DOCTYPE interpromatch SYSTEM "match_complete.dtd">\n')
xml_output.write('<interpromatch>\n')

# Load the idmapping entires for human proteins
PDBMapProtein.load_idmapping(defaults['idmapping'])
# PDBMapSwiss.load_swiss_INDEX_JSON(defaults['swiss_dir'],defaults['swiss_summary'])

xml_added = 0
xml_skipped = 0
last_added = 0
last_skipped = 0


def iterate_xml() -> Generator[etree._Element, None, None]:
    """
    With each iteration of this generator
    we return a depth=1 subtree of the XML for exanimation by the caller.
    """

    LOGGER.info("Opening %s for iterative parsing", xml_input_file)
    xml_input_file_gz = gzip.open(xml_input_file, "rb")

    # Open up a parser that does not attempt to load the entire giant file into RAM
    doc = etree.iterparse(xml_input_file_gz, events=('start', 'end'))
    _, root = next(doc)
    start_tag = None
    for event, _element in doc:
        # This logic is not clear - but it ensures
        # that only "protein" and other depth=1
        # xml sub-trees are "yieleded" below
        # Importantly, this approach prevents memory leaks
        if event == 'start' and start_tag is None:
            start_tag = _element.tag
        if event == 'end' and _element.tag == start_tag:
            yield _element
            start_tag = None
            root.clear()


for element in iterate_xml():
    if element.tag == 'protein':
        # import pdb; pdb.set_trace()
        unp = element.get('id')
        # If the base part of the uniprot identifier is in the KB based on our load of idmapping, then this is human
        isHumanUnp = PDBMapProtein.unp2uniprotKB(unp.split('-')[0])
        if isHumanUnp:
            xml_output.write(etree.tostring(element).decode('utf-8'))
            LOGGER.info("Added human unp %s to xml", unp)
            xml_added += 1
        else:
            xml_skipped += 1
    else:
        xml_output.write(etree.tostring(element).decode('utf-8'))
        LOGGER.info("Added non-protein XML %s" % etree.tostring(element))
        xml_added += 1

    # Provide an update every 100,000 or so so user sees progress
    if last_added != xml_added or (last_skipped != (xml_skipped // 100000)):
        LOGGER.info("Xml added %d, Skipped: %d" % (xml_added, xml_skipped))
        last_added = xml_added
        last_skipped = xml_skipped // 100000

LOGGER.info("Closing output file: %s", xml_output_file)
xml_output.write('</interpromatch>\n')

LOGGER.info("Ending successfully with Xml added %d, Skipped: %d", xml_added, xml_skipped)
