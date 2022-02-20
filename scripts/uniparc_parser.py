#!/usr/bin/env python
# Parses a (large!) uniparc cross reference file
# and uploads data to mysql database

from collections import defaultdict
import sys, os, glob, re, gzip

# Parse config file for database parameters
import argparse, configparser

cmdline_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
cmdline_parser.add_argument("uniparc_filename", nargs="?", help="Large uniparc filename",
                            default='/dors/capra_lab/data/uniprot/current/uniparc_active.fasta.gz')
cmdline_parser.add_argument("no_per_insert", nargs="?", help="Count of uniparcs to pass to INSERT(IGNORE)",
                            default=5000, type=int)
cmdline_parser.add_argument("-c", "--conf_file", required=True, help="Specify database config file", metavar="FILE")
cmdline_parser.add_argument("--skip_until_id", required=False, help="To speed restarts, specify a 'known uploaded through' UNIPARC id", type=str, metavar="str")
args = cmdline_parser.parse_args()

config = configparser.ConfigParser()
config.read([args.conf_file])

configdict = dict(config.items("Genome_PDB_Mapper"))

# Connect to the databsae
import MySQLdb

con = MySQLdb.connect(host=configdict['dbhost'], user=configdict['dbuser'],
                      passwd=configdict['dbpass'], db=configdict['dbname'])

import hashlib


def flush_uniparcs(uniparc_dict):
    if not uniparc_dict:
        return
    insert_list = []
    first_flag = True
    for uniparc_id in uniparc_dict:
        insert_list.append(
            (uniparc_id, hashlib.md5(uniparc_dict[uniparc_id].encode('UTF-8')).hexdigest(), uniparc_dict[uniparc_id]))

    ## Upload to database
    c = con.cursor()
    sql = "INSERT IGNORE INTO pdbmap_v14.Uniparc (uniparc,md5sum,fasta) VALUES (%s, %s, %s)"

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


with gzip.open(args.uniparc_filename, "rt") as f:
    print("File %s opened successfully" % args.uniparc_filename)
    cur_uniparc = ''  # The Uniparc UPI... identifier we are parsing sequence for
    cur_fasta = ''    # The amino acid letters associated with the Uniparc Id
    cur_is_active = False # True for most entries, where "status=active"
    uniparc_dict = {}
    skipping = True if args.skip_until_id else False
    print("Skipping until uniparc ID %s is encountered" % args.skip_until_id)
    for line in f:
        # When we get to the next >UPI... string
        if len(line) > 10 and line[0] == '>':
            if cur_fasta and cur_uniparc and cur_is_active:
                # Add the UPI... we've been working on to our dict
                if skipping:
                    if args.skip_until_id == cur_uniparc:
                        skipping = False

                if not skipping: # This is usual case - add the uniparc ID and sequence to SQL
                    uniparc_dict[cur_uniparc] = cur_fasta;
                
                    # Flush to SQL if we reach the buffer size limit
                    if len(uniparc_dict) >= args.no_per_insert:
                        flush_uniparcs(uniparc_dict)
                        uniparc_dict = {}

            assert (line[1:4] == "UPI")
            cur_uniparc = line[1:14]  # UPI + 10 character unique ID make the identifier
            cur_fasta = ' '
            cur_is_active = 'status=active' in line[15:]
        else:
            # Keep adding the amino acid letters to our growing
            # uniparc entry
            cur_fasta += line.rstrip()

    if len(uniparc_dict):
        flush_uniparcs(uniparc_dict)
    uniparc_dict = {}

sys.exit(0)
