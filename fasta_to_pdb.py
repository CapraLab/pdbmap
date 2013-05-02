#!/usr/bin/python

import os,sys,sqlite3,MySQLdb
from warnings import filterwarnings	# Disable MySQL warnings
filterwarnings('ignore',category=MySQLdb.Warning)

# Define an in-memory temporary sqlit database
con = sqlite3.connect(':memory:')
c = con.cursor()
c.execute("CREATE TABLE fasta_pep (unp VARCHAR(20),pepindex INT,col INT,aa1 VARCHAR(1))")
c.execute("CREATE TABLE sdp_scores (unp VARCHAR(20),pepindex INT,sdp_score DOUBLE)")

# User provides the fasta and sdpfilenames
fastaname = sys.argv[1]
sdpname = sys.argv[2]
DBNAME = sys.argv[3]
fasta_map = {}

# Build the column->pepindex mapping for each UNP
with open(fastaname,'r') as fin:
	unp = ''
	# For each row
	for line in fin:
		if line[0]=='>':
			unp = line[1:].split('_')[0]
		else:
			# For each column
			pepindex = 0	#pepindex 1 indexed, columns 0 indexed
			for col,aa1 in enumerate(line):
				if aa1 != '-':
					pepindex += 1
					c.execute("INSERT INTO fasta_pep VALUES (?,?,?,?)",(unp,pepindex,col,aa1))

# Build a database mapping UNP+seqres -> SDP score
with open(sdpname,'r') as fin:
	for line in fin:
		if line[0] == '#':
			continue
		col = int(line.split()[0])
		score = line.split()[1]
		if score == 'None':
			score = None
		else:
			score = float(score)
		
		# Get all UNP+Seqres associated with this column:
		c.execute("SELECT unp,pepindex FROM fasta_pep WHERE col=?",(col,))
		results = c.fetchall()
		for res in results:
			unp = str(res[0])
			pepindex = int(res[1])
			# For each UNP+SEQRES, add the score to the database
			c.execute("INSERT INTO sdp_scores VALUES (?,?,?)",(unp,pepindex,score))

# Upload the SDP scores to the MySQL database
MySQLcon = MySQLdb.connect(host='loki.mc.vanderbilt.edu',user='sivleyrm',passwd='19jcpkxw88',db=DBNAME)
MySQLc = MySQLcon.cursor()



c.execute("SELECT unp,pepindex,sdp_score FROM sdp_scores")
results = c.fetchall()
with open('sdp_scores.csv','w') as fout:
	for res in results:
		fout.write("%s,%d,%s\n"%res)
MySQLc.execute("LOAD DATA LOCAL INFILE 'sdp_scores.csv' INTO TABLE sdp_scores FIELDS TERMINATED BY ',' IGNORE 1 LINES")
MySQLcon.commit()

con.close()
MySQLcon.close()