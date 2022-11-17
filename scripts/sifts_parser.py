#!/usr/bin/env python
"""Parses a directory of SIFTS XML files
   and uploads data to mysql database.  It needs some ongoing work think.
   Better integration into PDBMap architecture, etc

   2019-09-03  This script is messed up because it also does REST API calls
   to populate two other tables.  Needs to be cleaned up badly

   2019-10-28  Passing legacy_xml gives this script its old behavior which
   saves sifts information for uniprot canonical IDs only"""

from lxml import etree, objectify
from collections import defaultdict
import sys, os, glob, re

import logging
from logging import handlers

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



rootdir_log_filename = "sifts_parser.log"
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
parser.add_argument("xmldir", nargs="?", help="Sifts data XML directory",
                    default=config_dict['sifts'] if 'sifts' in config_dict else None)
parser.add_argument("-x", "--legacy_xml", action='store_true',
                    help="Load sifts .xml to align canonical uniprot IDs for all pdbids ")
parser.add_argument("-b", "--best_isoforms", action='store_true',
                    help="Load 'best' isoforms data from sifts rest api")
parser.add_argument("-a", "--all_isoforms", action='store_true',
                    help="Load alignments for all isoforms from sifts rest api")
parser.add_argument("-p", "--pdb",
                    help="Load alignments for a single pdb for all isoforms in the sifts RESTapi")
args = parser.parse_args(remaining_argv)
args.conf_file = conf_file

assert args.legacy_xml or args.best_isoforms or args.all_isoforms

# Check that all parameters were specified
if not all(vars(args)):
    LOGGER.critical("Must provide database information and XML directory.")
    sys.exit(1)

if not args.xmldir:
    LOGGER.critical("xmldir must be supplied on command line or with sifts entry in -c config file\n");
    sys.exit(1)

from pdbmap import PDBMapProtein

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


def sifts_call_rest_api(REST_request, pdbid):
    try:
        LOGGER.debug("REST request: %s" % REST_request % pdbid)
        resp = requests.get(REST_request % pdbid)

    except requests.exceptions.RequestException as e:  # This is the correct syntax
        print(str(e))
        sys.exit(1)

    if resp.status_code == 404:
        return None

    if resp.status_code != 200:
        # This means something went wrong.
        raise APIError(REST_request + " {}".format(resp.status_code))

    resp_dict = resp.json()
    if not pdbid.lower() in resp_dict:
        raise APIError("fundamental problem in returned json format from sifts {}".format(str(resp)))

    return resp_dict


def sifts_get_best_isoforms(pdbid):
    # Documented at https://www.ebi.ac.uk/pdbe/api/doc/sifts.html
    return sifts_call_rest_api("https://www.ebi.ac.uk/pdbe/api/mappings/isoforms/%s", pdbid)


def sifts_get_all_isoforms(pdbid):
    # Documented at https://www.ebi.ac.uk/pdbe/api/doc/sifts.html
    return sifts_call_rest_api("https://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/%s", pdbid)


# isoform_json = sifts_get_all_isoforms('6CET')
# pp = pprint.PrettyPrinter()
# pp.pprint(isoform_json)

# Connect to the databsae
import MySQLdb

con = MySQLdb.connect(host=args.dbhost, user=args.dbuser,
                      passwd=args.dbpass, db=args.dbname)


# Increase maximum packet size for this connection
# Oops - not allowed under Redhat 7 new server!
# c = con.cursor()
# c.execute("SET GLOBAL max_allowed_packet=512000000")
# c.close()


def load_idmapping_pdbids():
    # Get a list of all pdb files known to the uniprot curated idmapping file
    sql_get_unique_pdbids = "SELECT distinct ID FROM pdbmap_v14.Idmapping where ID_type = 'pdb'"
    cursor = con.cursor()

    number_of_rows = 0;
    try:
        number_of_rows = cursor.execute(sql_get_unique_pdbids)
        pdbids = cursor.fetchall()
        cursor.close()
        return pdbids
    except:
        LOGGER.exception("Failed to execute: %s", sql)


def json_to_INSERTs(pdbid, unp, pdb_unp_dict):
    INSERT_list = []
    for mapping in pdb_unp_dict['mappings']:
        INSERT = {'pdbid': pdbid,
                  'uniprot_acc': unp,
                  'identifier': pdb_unp_dict['identifier'],
                  'name': pdb_unp_dict['name'],
                  'mapping_pdb_chain': mapping['chain_id'],
                  'mapping_entity_id': mapping['entity_id'],

                  'mapping_start_author_residue_number': mapping['start']['author_residue_number'],
                  'mapping_start_author_insertion_code': mapping['start']['author_insertion_code'],
                  'mapping_start_residue_number': mapping['start']['residue_number'],

                  'mapping_end_author_residue_number': mapping['end']['author_residue_number'],
                  'mapping_end_author_insertion_code': mapping['end']['author_insertion_code'],
                  'mapping_end_residue_number': mapping['end']['residue_number'],

                  'mapping_pdb_start': mapping['pdb_start'],
                  'mapping_pdb_end': mapping['pdb_end'],
                  'mapping_unp_start': mapping['unp_start'],
                  'mapping_unp_end': mapping['unp_end'],

                  'mapping_struct_asym_id': mapping['struct_asym_id'],
                  'mapping_seq_identity': mapping['identity']
                  }
        INSERT_list.append(INSERT)
    return INSERT_list


pdbs_processed_count = 0
if args.best_isoforms or args.all_isoforms:
    all_or_best = 'all' if args.all_isoforms else 'best'
    # with open("sifts_retry.txt") as f:
    #    pdbids = [(x.strip(),) for x in f.readlines()]
    #    print("%d pdbids read from sifts_retry.txt"%len(pdbids))

    # pdbs = ['6CES','1M6D']

    cursor = con.cursor()
    pdb_list = []
    if args.pdb:
        pdb_list = args.pdb.split(',')
    else:
        pdbid_rows_from_idmapping = load_idmapping_pdbids()
        pdb_list = [pdbid_row[0] for pdbid_row in pdbid_rows_from_idmapping]

    for pdbid in pdb_list:
        LOGGER.info("Processing pdbid %d/%d %s" % (pdbs_processed_count + 1, len(pdb_list), pdbid))
        isoform_json = sifts_get_all_isoforms(pdbid) if args.all_isoforms else sifts_get_best_isoforms(pdbid)
        if not isoform_json:
            LOGGER.warning("NO SIFTS REPONSE FOR %s ?Retained in uniprot ID Mapping; removed from PDB?" % pdbid)
            continue
        # assert(pdb.lower() in isoform_json)
        # print(isoform_json)
        INSERTs = []
        for pdb_dict in isoform_json:
            # print (isoform_json[pdb_dict])
            for uniprot_dict in isoform_json[pdb_dict]:
                for unp in isoform_json[pdb_dict][uniprot_dict]:
                    # pp = pprint.PrettyPrinter()
                    # pp.pprint(isoform_json[pdb_dict][uniprot_dict][unp])
                    INSERTs.extend(json_to_INSERTs(pdbid, unp, isoform_json[pdb_dict][uniprot_dict][unp]))

        for INSERT in INSERTs:
            placeholders = ', '.join(['%s'] * len(INSERT))
            columns = ', '.join(list(INSERT.keys()))
            sql = "INSERT IGNORE INTO pdbmap_v14.sifts_mappings_pdb_uniprot_%s_isoforms ( %s ) VALUES ( %s )" % (
            all_or_best, columns, placeholders)
            # cursor.executemany(sql, INSERTs)
            cursor.execute(sql, INSERT.values())
            con.commit()
        pdbs_processed_count += 1

    cursor.close()

    print("Processed %d pdbids" % len(pdb_list))
    sys.exit(1)

if args.legacy_xml:
    LOGGER.info("Converting *.xml.gz files to SQL in %s" % args.xmldir)
    # Parse the split XML files
    # for xmlfile in ["%s/1y97.xml.gz"%args.xmldir.rstrip('/')]:
    if args.pdb:
        xmlfile_list = [os.path.join(args.xmldir.rstrip('/'), pdbid.lower() + ".xml.gz") for pdbid in
                        args.pdb.split(',')]
    else:
        xmlfile_list = glob.glob("%s/*.xml.gz" % args.xmldir.rstrip('/'))

    for xmlfile in xmlfile_list:
        print("Parsing %s..." % xmlfile)
        parser = etree.XMLParser(remove_blank_text=True)
        tree = etree.parse(xmlfile, parser)
        root = tree.getroot()

        # Remove XML Namespaces
        for elem in root.getiterator():
            if not hasattr(elem.tag, 'find'): continue
            i = elem.tag.find('}')
            if i >= 0:
                elem.tag = elem.tag[i + 1:]
        objectify.deannotate(root, cleanup_namespaces=True)

        ## Annotate PDB residues
        rlist = []
        # Iterate over PDB chains
        for chain in root.findall("entity"):
            if chain.get("type") == "protein":
                # if chain.get("entityId") == 'S':
                #    import pdb; pdb.set_trace();
                # Iterate over SIFTS annotation segments
                for s in chain.getchildren():
                    # Iterate over segment residues
                    for r in s.find("listResidue"):
                        # Parse residue-level annotations
                        res = {"PDBe_dbResNum": r.get("dbResNum"), "pdbid": "", "pdb_chain": "", "pdb_resnum": None,
                               "pdb_icode": "", "pdb_resname": "",
                               "uniprot_acc": "", "uniprot_resnum": None, "uniprot_resname": "",
                               "ncbi": "", "pfam": [], "cath": [], "scop": [], "interpro": [],
                               "sscode": "", "ssname": ""}
                        for db in r.findall("crossRefDb"):
                            if db.get("dbSource") == "PDB":
                                res["pdbid"] = db.get("dbAccessionId")
                                res["pdb_chain"] = db.get("dbChainId")
                                # Careful, the dbResNum will be "null" in cases where there is no experimental observation
                                pdb_resnum = db.get("dbResNum")
                                res["pdb_icode"] = ""
                                if (not pdb_resnum) or (pdb_resnum == "null"):
                                    pdb_resnum = None
                                elif not pdb_resnum[-1].isdigit():  # If last pos is not a digit, it is insert code
                                    res["pdb_icode"] = pdb_resnum[-1]
                                    pdb_resnum = pdb_resnum[:-1]
                                res['pdb_resnum'] = pdb_resnum
                                res["pdb_resname"] = db.get("dbResName")
                            elif db.get("dbSource") == "UniProt":
                                res["uniprot_acc"] = db.get("dbAccessionId")
                                res["uniprot_resnum"] = db.get("dbResNum")
                                res["uniprot_resname"] = db.get("dbResName")
                            elif db.get("dbSource") == "NCBI":
                                res["ncbi"] = db.get("dbAccessionId")
                            elif db.get("dbSource") == "Pfam":
                                res["pfam"].append(db.get("dbAccessionId"))
                            elif db.get("dbSource") == "CATH":
                                res["cath"].append(db.get("dbAccessionId"))
                            elif db.get("dbSource") == "SCOP":
                                res["scop"].append(db.get("dbAccessionId"))
                            elif db.get("dbSource") == "InterPro":
                                res["interpro"].append(db.get("dbAccessionId"))
                        for rd in r.findall("residueDetail"):
                            if db.get("property") == "codeSecondaryStructure":
                                res["sscode"] = db.text
                            elif db.get("property") == "nameSecondaryStructure":
                                res["ssname"] = db.text
                        # Collapse lists to csv string
                        res["pfam"] = ",".join(res["pfam"])
                        res["cath"] = ','.join(res["cath"])
                        res["scop"] = ','.join(res["scop"])
                        res["interpro"] = ','.join(res["interpro"])
                        # Add residue to list of parsed residues
                        rlist.append(res)

        if len(rlist) < 1:
            LOGGER.warning("No alignments were parsed from %s" % xmlfile)
        else:
            # Sanity check that the pdbid elements match the filename
            assert rlist[0]['pdbid'] in xmlfile, "Crazy - xml file refers to a pdb not the filename %s" % xmlfile
            ## Upload to database
            c = con.cursor()
            sql = """DELETE FROM sifts_legacy_xml WHERE pdbid=%s"""
            try:
                rows_affected = c.execute(sql, (rlist[0]['pdbid'],))
                con.commit()
                if rows_affected > 0:
                    LOGGER.info("PDB %s.  Removed %d stale SIFTS xml rows", rlist[0]['pdbid'], rows_affected)
            except:
                con.rollback()
                LOGGER.exception("Failure on attempt to DELETE old pdb values\n%s", c._last_executed)
                raise

            sql = """INSERT IGNORE INTO sifts_legacy_xml
            (PDBe_dbResNum,pdbid,pdb_chain,pdb_resnum,pdb_icode,pdb_resname,uniprot_acc,uniprot_resnum,
             uniprot_resname,ncbi,pfam,cath,scop,interpro,sscode,ssname)
            values 
            (%(PDBe_dbResNum)s,%(pdbid)s,%(pdb_chain)s,%(pdb_resnum)s,%(pdb_icode)s,%(pdb_resname)s,
             %(uniprot_acc)s,%(uniprot_resnum)s,%(uniprot_resname)s,
             %(ncbi)s,%(pfam)s,%(cath)s,%(scop)s,%(interpro)s,
             %(sscode)s,%(ssname)s);"""

            LOGGER.debug(sql % rlist[0])

            LOGGER.info("%d rows of residue alignments to add via SQL", len(rlist))
            rows_affected = 0
            try:
                rows_affected = c.executemany(sql, rlist)
                con.commit()
                LOGGER.info("Uploaded and committed! %d rows affected" % rows_affected)
            except:
                con.rollback()
                LOGGER.exception("Failed to upload rows.\n%s", c._last_executed)
                raise
            if rows_affected != len(rlist):
                LOGGER.critical("SOMEHOW, in %s not all rows added\n %s to %s",
                                str(rlist[0]['pdbid']),
                                str(rlist[0]),
                                str(rlist[-1]))
        c.close()
