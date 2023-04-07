#!/usr/bin/env python
"""Parses a directory of SIFTS XML files
   and uploads data to mysql database.  It needs some ongoing work think.
   Better integration into PDBMap architecture, etc

   2023-02-13  The uniprot/HUMAN_9606_idmapping_sprot.dat.gz files cross-references
   53,000 PDB structures to SwissProt human uniprot IDs.  HOWEVER, going to sifts directly with
   our curated SwissProt IDs returns over 110,000 transcript-aligned PDB structures.  This,
   it is helpful for us to preload all of these.

   A challenge with the REST-API approach is that it is prone to faults, and crashes.  It takes time.
   A lot of the code is about making sure to "not stop" and gather ALL the data before crashing, then
   Attempting to not repeat loads of previously gathered data.

   2019-09-03  This script is messed up because it also does REST API calls
   to populate two other tables.  Needs to be cleaned up badly

   2019-10-28  Passing legacy_xml gives this script its old behavior which
   saves sifts information for uniprot canonical IDs only"""

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

pdbs_processed_count = 0
pdbs_added_count = 0
pdbs_deleted_count = 0


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
                    help="Use /tmp pdb list from --all_isoforms and Load 'best' isoforms data from sifts rest api")
parser.add_argument("-a", "--all_isoforms", action='store_true',
                    help="Create /tmp pdb list and Load alignments for all isoforms from sifts rest api")
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

def sifts_call_rest_api(REST_request, pdb_or_uniprot_id):
    resp = None
    for attempt in range(20):
        try:
            # LOGGER.debug("REST request: %s" % REST_request % pdb_or_uniprot_id)
            resp = requests.get(REST_request % pdb_or_uniprot_id)
            break  # Success

        except requests.exceptions.RequestException as e:  # This is the correct syntax
            LOGGER.error("Attempt %d failed: %s", attempt + 1, str(e))
            time.sleep((attempt+1) * 1.5)

    if attempt >= 19:
        err_failed = "Failed to reach sifts https after %d tries" % attempt
        LOGGER.critical(err_failed)
        sys.exit(err_failed)

    # These calls return 404 in case the response is _mostly_ empty.
    # This is an odd thing, because even with 404, I think often
    # the content looks like reasonable JSON
    # 404 is NOT an error in this SIFTS API.
    if resp.status_code == 404:
        return None

    if resp.status_code != 200:
        # This means something went wrong.
        raise APIError(REST_request + " {}".format(resp.status_code))

    resp_dict = resp.json()
    # The format of the returned json varies depending on whether we asked for a uniprot ID or a PDBID.
    # For _todaY_, we assume a PDB id if len is <= 4
    if len(pdb_or_uniprot_id) <= 4:
        if not pdb_or_uniprot_id.lower() in resp_dict:
            raise APIError("fundamental problem in returned json format from sifts {}".format(str(resp)))
    else: # It's a uniprot ID - always upper case.
        if not pdb_or_uniprot_id in resp_dict:
            raise APIError("fundamental problem in returned json format from sifts {}".format(str(resp)))

    return resp_dict


def sifts_get_best_isoforms(pdb_id: str):
    # Documented at https://www.ebi.ac.uk/pdbe/api/doc/sifts.html
    return sifts_call_rest_api("https://www.ebi.ac.uk/pdbe/api/mappings/isoforms/%s", pdb_id)


def sifts_get_all_isoforms(pdb_or_uniprot_id: str):
    # Documented at https://www.ebi.ac.uk/pdbe/api/doc/sifts.html
    return sifts_call_rest_api("https://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/%s", pdb_or_uniprot_id)


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

# Get a list of EVERY uniprot ID mentioned in the curated human list
# We can run queries to sifts to figure out which PDBs are alighed to these
def load_idmapping_uniprot_ids() -> List[Tuple]:
    return sql_select_where("SELECT DISTINCT unp FROM Idmapping where ID_type = 'UniParc'")



def load_idmapping_pdbids() -> List[Tuple]:
    # Get a list of all pdb files known to the uniprot curated idmapping file
    return sql_select_where("SELECT DISTINCT ID FROM Idmapping where ID_type = 'pdb'")

def delete_pdb_from_mapping_table(con_cursor: MySQLdb.cursors.Cursor,
                                  pdb_id: str,
                                  sql_mapping_table: str) -> None:
    """
    Delete all entries for a specific pdbid from the given table.
    This uses a cursor created in the calling function
    """

    # CAREFUL: You cannot substitute in a table name as a parameter.
    # So, we're going to "hard code" the table name as a long normal string
    # and just paste in the pdb_id at the concursor.execute() part
    global pdbs_deleted_count
    delete_a_pdb_sql = 'DELETE FROM ' + sql_mapping_table + """ WHERE pdbid=%s"""

    rows_affected = con_cursor.execute(delete_a_pdb_sql, (pdb_id,))
    if rows_affected >= 0:
        LOGGER.info("PDB %s.  DELETEd %d stale SIFTS rows from %s", pdb_id, rows_affected, sql_mapping_table)
        if rows_affected > 0:
            pdbs_deleted_count += 1


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




if args.best_isoforms or args.all_isoforms:
    all_or_best = 'all' if args.all_isoforms else 'best'
    # with open("sifts_retry.txt") as f:
    #    pdbids = [(x.strip(),) for x in f.readlines()]
    #    print("%d pdbids read from sifts_retry.txt"%len(pdbids))

    # pdbs = ['6CES','1M6D']

    save_uniprot_progress = {}

    if all_or_best == 'all': # Then we should pick up additional PDB ids from SIFTS
        # The idea here is that the tmp file will go away periodically...  which is perfect
        # I just want to avoid endless reloading from the REST API which is error prone
        save_uniprot_progress_filename = '/tmp/save_uniprot_progress_all_isoforms.json.gz'
        try:
            with gzip.open(save_uniprot_progress_filename,'rt') as f:
                save_uniprot_progress = json.load(f)

        except:
            save_uniprot_progress = {}

        LOGGER.info('%d previously gathered uniprot IDs loaded from', len(save_uniprot_progress))

        last_saved_len =     len(save_uniprot_progress)

    def save_uniprot_progress_to_tmp():
        LOGGER.info("Writing %d records to %s", len(save_uniprot_progress), save_uniprot_progress_filename)
        with gzip.open(save_uniprot_progress_filename,'wt') as f:
            json.dump(save_uniprot_progress, f, indent=1)


    sql_mapping_table = "sifts_mappings_pdb_uniprot_%s_isoforms" % all_or_best

    pdb_set = set()
    if args.pdb:
        pdb_set = {pdb_id for pdb_id in args.pdb.split(',')}
    else:
        # First let's get PDB IDs for those PDBs associated in the IDmapping file
        pdbid_rows_from_idmapping = load_idmapping_pdbids()
        pdb_set = {pdbid_row[0] for pdbid_row in pdbid_rows_from_idmapping}

        # Let's add in any PDBs that might be associated with our uniprot IDs
        # at SIFTS
        uniprot_rows_from_idmapping = load_idmapping_uniprot_ids()
        unp_list = [uniprot_id_row[0] for uniprot_id_row in uniprot_rows_from_idmapping]

        sifts_pdb_set = set()
        # Let's try to get ALL the PDBs that are aligned by SIFTS
        if all_or_best == 'best':
            # We should mine the PDB IDs from a prior 'all' run of this script.  Read the SQL database for a list
            select_distinct_PDBs = 'select distinct pdbid from sifts_mappings_pdb_uniprot_all_isoforms'
            pdb_list_from_all = sql_select_where(select_distinct_PDBs)
            LOGGER.info("%d pdb ids loaded from sifts_mappings_pdb_uniprot_all_isoforms", len(pdb_list_from_all))
            for pdb_row in pdb_list_from_all:
               sifts_pdb_set.add(pdb_row[0]) # Each list element is a tuple from SQL
            LOGGER.info("%d pdbids loaded sifts_mappings_pdb_uniprot_all_isoforms", len(sifts_pdb_set))
        else:
            assert all_or_best == 'all'
            test_count = 0
            for unp in unp_list:
                pdbs_for_unp = []
                isoform_json = {}

                if unp in save_uniprot_progress:
                    isoform_json = save_uniprot_progress
                else:
                    isoform_json = sifts_get_all_isoforms(unp) if args.all_isoforms else sifts_get_best_isoforms(unp)
                if  isoform_json and 'PDB' in isoform_json[unp]:
                    pdb_dict = isoform_json[unp]['PDB']
                    for pdbid in pdb_dict.keys():
                        pdbs_for_unp.append(pdbid)
                        sifts_pdb_set.add(pdbid)
                    # If we are not retrieving the PDB ilst from the saved isoform json THEN
                    # we should add this restapi information to our saved file.
                    if not isoform_json is save_uniprot_progress:
                        save_uniprot_progress[unp] = isoform_json[unp]
                else: # Next time we save the progress, note that we have no PDBs for this uniprot ID
                    save_uniprot_progress[unp] = {'PDB': {}}
                if len(save_uniprot_progress) - last_saved_len >= 200:
                    save_uniprot_progress_to_tmp()
                    last_saved_len = len(save_uniprot_progress)

                LOGGER.info("%s xrefs to: %s", unp, str(pdbs_for_unp))
                if pdbs_for_unp:
                    test_count += 1
                    # if test_count > 5: # For development it can be helpful to cut this off.
                    #    break

            save_uniprot_progress_to_tmp()
            LOGGER.info("%d pdbs from sifts restAPI calls.  %d pdbs from idmapping", len(sifts_pdb_set), len(pdb_set))

    # Create the larger set of pdbs to round up alignments for
    pdb_set = sifts_pdb_set.union(pdb_set)

    reconnect_sql()

    for pdbid in pdb_set:
        LOGGER.info("Processing pdbid %d/%d %s" % (pdbs_processed_count + 1, len(pdb_set), pdbid))
        isoform_json = sifts_get_all_isoforms(pdbid) if args.all_isoforms else sifts_get_best_isoforms(pdbid)

        try:
            cursor = con.cursor()
        except:
            LOGGER.warning("Reconnecting SQL after inability to get cursor()")
            reconnect_sql()
            cursor = con.cursor()

        delete_pdb_from_mapping_table(
            cursor,
            pdbid,
            sql_mapping_table
        )


        if not isoform_json:
            LOGGER.warning("NO SIFTS REPONSE FOR %s ?Retained in uniprot ID Mapping; removed from PDB?", pdbid)
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
            sql = "INSERT IGNORE INTO %s ( %s ) VALUES ( %s )" % (
            sql_mapping_table, columns, placeholders)
            # cursor.executemany(sql, INSERTs)
            cursor.execute(sql, INSERT.values())

        LOGGER.info("%s: %d rows added", pdbid, len(INSERTs))

        if len(INSERTs) > 0:
            pdbs_added_count += 1
        con.commit()
        pdbs_processed_count += 1

        cursor.close()

    LOGGER.info("PDBs Added: %d  PDBs Deleted: %d  PDBs Processed: %d", \
                pdbs_added_count, pdbs_deleted_count, len(pdb_set))
    sys.exit(0)

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
            delete_pdb_from_mapping_table(
                c,
                rlist[0]['pdbid'],
                'sifts_legacy_xml'
            )

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
            except Exception as ex:
                con.rollback()
                LOGGER.exception("%s: Failed to upload rows.\n", ex)
                raise
            if rows_affected != len(rlist):
                LOGGER.critical("SOMEHOW, in %s not all rows added\n %s to %s",
                                str(rlist[0]['pdbid']),
                                str(rlist[0]),
                                str(rlist[-1]))
            pdbs_added_count += 1
        c.close()
