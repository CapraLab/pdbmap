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
query  = "CREATE TABLE IF NOT EXISTS sifts (pdbid VARCHAR(100), "
query += "chain VARCHAR(50), sp VARCHAR(100), "
query += "pdb_seqid INT, sp_seqid INT, "
query += "PRIMARY KEY(pdbid,chain,pdb_seqid), "
query += "KEY(sp,sp_seqid))"
c.execute(query)
c.close()

c = con.cursor()
with open('../data/sifts/pdb_chain_uniprot.tsv','rb') as fin:
  fin.readline(); fin.readline() # skip first two lines
  reader = csv.reader(fin,delimiter='\t')
  for row in reader:
    q = "INSERT IGNORE INTO sifts (pdbid,chain,sp,pdb_seqid,sp_seqid) VALUES "
    v = []
    pdbid,chain,sp = row[0:3]
    pdbid = pdbid.strip()
    try:
      res_beg,res_end,pdb_beg,pdb_end,sp_beg,sp_end = [int(x) for x in row[3:]]
    except:
      sys.stderr.write("icode error: %s,%s,%s\n"%(pdbid,chain,sp))
      continue
    res_range = range(res_beg,res_end+1)
    pdb_range = range(pdb_beg,pdb_end+1)
    sp_range  = range(sp_beg,sp_end+1)
    if len(res_range) != len(pdb_range) or len(pdb_range) != len(sp_range):
      sys.stderr.write("range mismatch: %s,%s,%s\n"%(pdbid,chain,sp))
      continue
    for i,seqid in enumerate(pdb_range):
      v.append("('%s','%s','%s',%d,%d)"%(pdbid,chain,sp,seqid,sp_range[i]))
    # Upload after each row
    q = q+','.join(v)
    try:
      c.execute(q)
    except:
      print row
      print q
      print c._last_executed
      raise
c.close()
con.close()