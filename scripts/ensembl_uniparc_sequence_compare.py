#!/usr/bin/env python
"""
An ongoing challenge with the pipeline is that many uniprot transcripts have no
ENST transcript ID, OR, that the transcript returns a different sequence
fomr the genome than the uniprot UniParc based transcript

This program reports these disconnects.
"""
import logging
from logging import handlers

from lxml import etree, objectify
from collections import defaultdict
import sys, os, glob, re
import gzip
from typing import List, Tuple
import time


# from logging import handlers
sh = logging.StreamHandler()
LOGGER = logging.getLogger()
LOGGER.addHandler(sh)



log_format_string = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
date_format_string = '%H:%M:%S'
log_formatter = logging.Formatter(log_format_string, date_format_string)

LOGGER.setLevel(logging.DEBUG)
sh.setLevel(logging.DEBUG)
sh.setFormatter(log_formatter)


rootdir_log_filename = "ensembl_vs_uniparc.log"
needRoll = os.path.isfile(rootdir_log_filename)
rootdir_fh = logging.handlers.RotatingFileHandler(rootdir_log_filename, backupCount=4)
formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)20s:%(lineno)d] %(message)s', datefmt="%H:%M:%S")
rootdir_fh.setFormatter(formatter)
rootdir_fh.setLevel(logging.INFO)
LOGGER.addHandler(rootdir_fh)

LOGGER.info("sys.argv=%s", str(sys.argv))

if needRoll:
    rootdir_fh.doRollover()

sys.stderr.write("Root (case) directory log file is %s\n" % rootdir_log_filename)


# Parse config file for database parameters
import argparse, configparser
import pprint

# For now, just get the conf_file paramter
# We return later for more
conf_parser = argparse.ArgumentParser(add_help=False)
conf_parser.add_argument("-c", "--conf_file",
                         help="Specify database config file", metavar="FILE")

args, remaining_argv = conf_parser.parse_known_args()
defaults = {
    "dbhost": None,
    "dbuser": None,
    "dbpass": None,
    "dbname": None,
    "xmldir": None
}

config_dict = {}
if args.conf_file:
    config = configparser.ConfigParser()
    config.read([args.conf_file])
    config_dict = dict(config.items("Genome_PDB_Mapper"))
    defaults.update(config_dict)


conf_file = args.conf_file
# Check for command line argument for the XML directory
parser = argparse.ArgumentParser(parents=[conf_parser],
                                 description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.set_defaults(**defaults)
args = parser.parse_args(remaining_argv)
args.conf_file = conf_file

# Get a list of EVERY uniprot ID mentioned in the curated human list
# We can run queries to sifts to figure out which PDBs are alighed to these
def load_idmapping_uniprot_ids() -> List[Tuple]:
    return sql_select_where("SELECT DISTINCT unp FROM Idmapping where ID_type = 'UniParc'")


from pdbmap import PDBMapProtein
from lib import PDBMapTranscriptUniprot
from lib import PDBMapTranscriptEnsembl

PDBMapProtein.load_idmapping(defaults['idmapping'])
PDBMapProtein.load_sec2prim(defaults['sec2prim'])

import requests
import json


class APIError(Exception):
    """An API Error Exception"""

    def __init__(self, status):
        self.status = status

    def __str__(self):
        return "APIError: status={}".format(self.status)





# isoform_json = sifts_get_all_isoforms('6CET')
# pp = pprint.PrettyPrinter()
# pp.pprint(isoform_json)

# Connect to the databsae
import MySQLdb

con = None
def reconnect_sql():
    global con
    if con is not None:
        try:
            con.close()
        except Exception as ex:
            pass

    con = MySQLdb.connect(host=args.dbhost, user=args.dbuser,
                      passwd=args.dbpass, db=args.dbname)

reconnect_sql()

import string

unp_enst_seq_mismatch = []
unp_enst_len_mismatch = []

with open("/tmp/unp_enst_mismatch.txt") as f:
    translator = str.maketrans('', '', "() '\n")
    for line in f.readlines():
        print(line)
        tuple_break = line.split(',')
        unp = tuple_break[0].translate(translator)
        enst = tuple_break[1].translate(translator)
        print("[%s] [%s]" % (unp, enst))
        enst_transcript = PDBMapTranscriptEnsembl(enst)
        unp_transcript = PDBMapTranscriptUniprot(unp)
        if enst_transcript.aa_seq == unp_transcript.aa_seq:
            sys.exit("HUGE FAIL - stop %s" % unp)
        try:

            if len(enst_transcript.aa_seq) != len(unp_transcript.aa_seq):
                unp_enst_len_mismatch.append(str((unp,enst)) + "%d vs %d len" % (
                    len(unp_transcript.aa_seq),  len(enst_transcript.aa_seq)
                ))
            else:
                unp_enst_seq_mismatch.append(str((unp,enst)) + str([(n+1,unp_transcript.aa_seq[n],enst_transcript.aa_seq[n]) for n in range(len(unp_transcript.aa_seq)) \
                                                                if unp_transcript.aa_seq[n] != enst_transcript.aa_seq[n]]))
        except:
            LOGGER.info("FAILURE - CAREFUL on %s %s", unp,enst)
            pass


with open('/tmp/unp_enst_len_mismatch.txt','w') as f:
    for unp_enst in unp_enst_len_mismatch:
        f.write("%s\n" % str(unp_enst))
with open('/tmp/unp_enst_seq_mismatch.txt','w') as f:
    for unp_enst in unp_enst_seq_mismatch:
        f.write("%s\n" % str(unp_enst))
sys.exit(0)


# Increase maximum packet size for this connection
# Oops - not allowed under Redhat 7 new server!
# c = con.cursor()
# c.execute("SET GLOBAL max_allowed_packet=512000000")
# c.close()

def sql_select_where(sql_query) -> List[Tuple]:
    cursor = con.cursor()

    number_of_rows = 0;
    try:
        number_of_rows = cursor.execute(sql_query)
        ids = cursor.fetchall()
        cursor.close()
        return ids
    except:
        LOGGER.exception("Failed to execute: %s", sql_query)

    cursor.close()


uniprot_rows_from_idmapping = load_idmapping_uniprot_ids()
unp_list = [uniprot_id_row[0] for uniprot_id_row in uniprot_rows_from_idmapping]

LOGGER.info("%d unique uniprot IDs", len(unp_list))
# Next, get rid of UNPs that both
# 1) are nonspecific (don't have a dash)
# AND
# 2) Have a specific (isoform specific, dashed UNP)
# that is better




better_unp_list = []
for unp in unp_list:
    if not '-' in unp: # If this is a non-isoformspecific ID
        # Is there a specific one that matches
        best_unp = PDBMapProtein.best_unp(unp)
        # Skip this dash-less unp if there is a better isoformspecific
        if best_unp != unp:
            continue
    better_unp_list.append(unp)

LOGGER.info("Best unps count=%d", len(better_unp_list))

unps_missing_enst = []
unp_enst_mismatches = []
unps_OK = []
for unp in better_unp_list:
    enst_list = PDBMapProtein.unp2enst(unp)
    if not enst_list:
        unps_missing_enst.append(unp)
    else:
        unp_transcript = PDBMapTranscriptUniprot(unp)

        for enst_id in enst_list:
            enst_transcript = PDBMapTranscriptEnsembl(enst_id)
            if unp_transcript.aa_seq == enst_transcript.aa_seq:
                unps_OK.append((unp,enst_id))
                continue
            unp_enst_mismatches.append((unp,enst_id))
    # EARLY break for development testing only
    # if (len(unps_OK) + len(unps_missing_enst) + len(unp_enst_mismatches) > 200):
    #    break

with open('/tmp/unps_OK.txt','w') as f:
    for unp_OK in unps_OK:
        f.write("%s\n" % str(unp_OK))

with open('/tmp/unp_enst_mismatch.txt','w') as f:
    for unp_enst_mismatch in unp_enst_mismatches:
        f.write("%s\n" % str(unp_enst_mismatch))

with open('/tmp/unps_missing_enst.txt', 'w') as f:
    for unp_missing_enst in unps_missing_enst:
        f.write("%s\n" % unp_missing_enst)

logging.info("See outputs in /tmp/unp*.txt")

