#!/usr/bin/env python2.7

import csv,sys
import MySQLdb, MySQLdb.cursors
from warnings import filterwarnings,resetwarnings
filterwarnings('ignore', category = MySQLdb.Warning)

# Hard code database credentials and connect for intermediate tables
con = MySQLdb.connect(host='gwar-dev',user='mike',
                            passwd='cheezburger',db='pdbmap_v10',
                            cursorclass=MySQLdb.cursors.Cursor)
c = con.cursor()
query  = "CREATE TABLE IF NOT EXISTS pfam (pdbid VARCHAR(100), "
query += "chain VARCHAR(50), seqstart INT, "
query += "seqend INT, acc VARCHAR(50), "
query += "name VARCHAR(50), description TEXT, evalue DOUBLE, "
query += "PRIMARY KEY(pdbid,chain,seqstart,seqend,acc), "
query += "KEY(acc), KEY(name), KEY(evalue))"
print "Creating pfam table if necessary...",
c.execute(query)
c.close()
print "done"

c = con.cursor()
query  = "LOAD DATA LOCAL INFILE '../data/pfam/pdb_pfam_mapping.txt' "
query += "INTO TABLE pfam "
query += "FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n' "
query += "IGNORE 1 LINES"
print "Loading pfam into database...",
c.execute(query)
c.close()
print "done"