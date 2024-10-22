#!/usr/bin/env python3
# Make RESTAPI calls to fetch UNIPARC sequences that are needed by uniprot entries in the Idmapping table
# but which are not already in the Uuniparc or UniparcHuman tables

from collections import defaultdict
import sys, os, glob, re, gzip
import requests

# Parse config file for database parameters
import argparse, configparser

cmdline_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
cmdline_parser.add_argument("uniparc_filename", nargs="?", help="Large uniparc filename",
                            default='/dors/capra_lab/data/uniprot/current/uniparc_active.fasta.gz')
cmdline_parser.add_argument("no_per_insert", nargs="?", help="Count of uniparcs to pass to INSERT(IGNORE)",
                            default=5000, type=int)
cmdline_parser.add_argument("-c", "--conf_file", required=True, help="Specify database config file", metavar="FILE")
cmdline_parser.add_argument("--skip_until_id", required=False, help="To speed restarts, specify a 'known uploaded through' UNIPARC id", type=str, metavar="str")
cmdline_parser.add_argument("--humanonly", required=False, help="Add only human uniparc IDs to UniparcHuman table", action='store_true')
args = cmdline_parser.parse_args()


config = configparser.ConfigParser()
config.read([args.conf_file])

configdict = dict(config.items("Genome_PDB_Mapper"))
uniparc_table_name = "UniparcHuman" if args.humanonly else "Uniparc"

import MySQLdb

con = MySQLdb.connect(host=configdict['dbhost'], user=configdict['dbuser'],
                      passwd=configdict['dbpass'], db=configdict['dbname'])

import hashlib

def call_rest_api(REST_request):
    resp = None
    for attempt in range(20):
        try:
            resp = requests.get(REST_request)
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

    resp_text  = resp.text
    # Sanity check first line of fasta file starts with the >
    assert(resp_text[0] == '>')

    # Now convert EOL delimited lines of FASTA format to list of lines
    resp_strings = resp_text.split('\n')
    aa_seq = ""
    for line_no in range(1,len(resp_strings)):
        if len(resp_strings[line_no]) > 0:
            aa_seq += resp_strings[line_no]

    return aa_seq




connection_cursor = con.cursor()

# When we have a null from SQL, that tells us the referenced uniparc row is missing for the 
# Idmapping row
sql_select_missing_uniparc_ids="""\
select DISTINCT Idmapping_uniparc.id FROM  (select * from Idmapping where id_type='Uniparc') AS Idmapping_uniparc 
LEFT OUTER JOIN %s on %s.uniparc=Idmapping_uniparc.id where id_type='Uniparc' and isnull(uniparc);
""" % (uniparc_table_name, uniparc_table_name)

connection_cursor.execute(sql_select_missing_uniparc_ids)
all_rows = connection_cursor.fetchall()

# all_rows = [('UPI000015D2DD',''),('UPI0000F23CBD','')]


if len(all_rows) < 1:
    print("All Uniparc IDs needed in Idmapping table have been already loaded into Uniparc table: %s" % uniparc_table_name)
    sys.exit(0)

print("%d Uniparc IDs must be added to Uniparc table [%s .. %s]" % (len(all_rows),all_rows[0][0],all_rows[-1][0]))

uniparc_dict = {}
insert_list = []
for row_tuple in all_rows:
    uniparc_id = row_tuple[0]
    restapi_uri = "https://rest.uniprot.org/uniparc/%s.fasta" % uniparc_id
    print("%s" % restapi_uri)

    aa_seq = call_rest_api(restapi_uri)
    insert_list.append(
        (uniparc_id, hashlib.md5(aa_seq.encode('UTF-8')).hexdigest(), str(aa_seq))
        )
    if len(insert_list) >= 1000:
        break

## Upload to database
c = con.cursor()
sql = "INSERT IGNORE INTO pdbmap_v14.%s (uniparc,md5sum,fasta) VALUES (%%s, %%s, %%s)" % (
    uniparc_table_name,)

try:
    print("Uploading %d items ranging %s to %s" % (len(insert_list), insert_list[0][0], insert_list[-1][0]),
          end='   ')
    rows_affected = c.executemany(sql, insert_list)
    con.commit()
    print("Uploaded %d!" % rows_affected)
except:
    con.rollback()
    print("Failed to upload rows.")
    print(c._last_executed)
    raise
c.close()

print("Please run this again until 0 uniparc IDs are added to the Uniparc table")


sys.exit(0)
