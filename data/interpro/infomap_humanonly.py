#!/usr/bin/env python2.7
# Chris Moth 2017-Nov 30
# Extract ONLY human entries from the ~150GB match_complete.xml
# The resulting match_human.xml file in the same directory (specify in -c configfile)
# is used to help generate domain graphics for isoform-specific uniprot ACs in the pipeline
# output reports

import argparse,ConfigParser
import traceback
import sys,os,csv,time,pdb,glob,gzip,shutil
import subprocess as sp
from multiprocessing import cpu_count
from lib import PDBMapIO,PDBMapParser,PDBMapStructure,PDBMapProtein
from lib import PDBMapAlignment,PDBMapData,PDBMapTranscript
from lib import PDBMapIntersect,PDBMapModel

import logging
from logging.handlers import RotatingFileHandler
from logging import handlers
logging.basicConfig(format='%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
    datefmt='%d-%m-%Y:%H:%M:%S',)

root  = logging.getLogger()
root.setLevel(logging.INFO)

fh = RotatingFileHandler("infomap_humanonly.log", maxBytes=(1048576*5), backupCount=7)
formatter = logging.Formatter('%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s')
fh.setFormatter(formatter)
fh.setLevel(logging.WARNING)
root.addHandler(fh)


from lxml import etree
import sys

# Setup the Config File Parser
conf_parser = argparse.ArgumentParser(add_help=False)
conf_parser.add_argument("-c", "--conf_file",
help="Specify config file", metavar="FILE", required=True)
args, remaining_argv = conf_parser.parse_known_args()
conf_file = args.conf_file
print 'Reading config file: %s'%conf_file
config = ConfigParser.SafeConfigParser()
config.read([conf_file])
defaults = {}
defaults.update(dict(config.items("Genome_PDB_Mapper")))
print 'Interpro directory: %s'%defaults['interpro_dir']

xml_input_file = defaults['interpro_dir'] + 'match_complete.xml'
xml_output_file = defaults['interpro_dir'] + 'match_humanonly.xml'

xml_output = open(xml_output_file,"w+")
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
def iterate_xml():
    doc = etree.iterparse(xml_input_file, events=('start', 'end'))
    _, root = next(doc)
    start_tag = None
    for event, element in doc:
        # This logic is not clear - but it ensures
        # that only "protein" and other depth=1
        # xml sub-trees are "yieleded" below
        # Importantly, this approach prevents memory leaks
        if event == 'start' and start_tag is None:
            start_tag = element.tag
        if event == 'end' and element.tag == start_tag:
            yield element
            start_tag = None
            root.clear()

for element in iterate_xml():
  if element.tag == 'protein':
    # import pdb; pdb.set_trace()
    unp = element.get('id')
    # If the base part of the uniprot identifier is in the KB based on our load of idmapping, then this is human
    isHumanUnp = PDBMapProtein.unp2uniprotKB(unp.split('-')[0])
    if isHumanUnp:
      xml_output.write(etree.tostring(element))
      print "Added human unp %s to xml"%unp
      xml_added += 1
    else:
      xml_skipped += 1 
  else:
    xml_output.write(etree.tostring(element))
    print "Added non-protein XML %s"%etree.tostring(element)
    xml_added += 1

  if last_added != xml_added or last_skipped != (xml_skipped / 100000):
    print "Xml added %d, Skipped: %d"%(xml_added,xml_skipped)
    last_added = xml_added
    last_skipped = xml_skipped / 100000

xml_output.write('</interpromatch>\n')

