#!/usr/bin/python

aa_code_map = {"ala" : "A",
				"arg" : "R",
				"asn" : "N",
				"asp" : "D",
				"asx" : "B",
				"cys" : "C",
				"glu" : "E",
				"gln" : "Q",
				"glx" : "Z",
				"gly" : "G",
				"his" : "H",
				"ile" : "I",
				"leu" : "L",
				"lys" : "K",
				"met" : "M",
				"phe" : "F",
				"pro" : "P",
				"ser" : "S",
				"thr" : "T",
				"trp" : "W",
				"tyr" : "Y",
				"val" : "V"}

species_map = {"PIG" : "sus_scrofa",
				"BAKER'S YEAST" : "saccharomyces_cerevisiae",
				"YEAST" : "saccharomyces_cerevisiae"}


def read_species(fin):
	for line in fin:
		if line[0:6] == "SOURCE" and line[11:26] == "ORGANISM_COMMON":
			species = line[28:].split(';')[0]
			if species in species_map:
				species = species_map[species]
			return species


def read_dbref(fin,c,pdb_id,species):
	flag = 0
	for line in fin:
		if line[0:5] == "DBREF" and line[26:32].strip() == "UNP":
			flag = 1

			chain = line[12]
			unp_id = line[33:41].strip()
			pdbstart = int(line[14:18])
			dbstart = int(line[55:60])
			offset = dbstart - pdbstart
			c.execute("INSERT INTO chains VALUES (?,?,?,?,?,?)",(chain,unp_id,pdb_id,pdbstart,offset,species))
			con.commit()
		elif flag:
			return

def read_seqres(fin,c):
	flag = 0
	index = 0
	for line in fin:
		if line[0:6] == "SEQRES":
			flag = 1

			chain = line[11]
			pep_subseq = line[19:70].split()
			for aa in pep_subseq:
				c.execute("INSERT INTO peptide VALUES (?,?,?)",(chain,index,aa))
				con.commit()
				index += 1
		elif flag:
			return 

def read_atom(fin,c):
	chainatom_coord_map = {}
	for line in fin:
		if line[0:4] == "ATOM":

			chain = line[21]
			serialnum = int(line[6:11])
			seqres = int(line[22:26])
			aminoacid = line[17:20]
			aminoacid_oneletter = aa_code_map[aminoacid.lower()]
			x = float(line[30:38])
			y = float(line[38:46])
			z = float(line[46:54])
			c.execute("INSERT INTO coords VALUES (?,?,?,?,?,?,?,?)",(chain,serialnum,seqres,aminoacid,aminoacid_oneletter,x,y,z))
			con.commit()

import os,sys,csv,sqlite3,MySQLdb
from warnings import filterwarnings	# Disable MySQL warnings
filterwarnings('ignore',category=MySQLdb.Warning)

con = sqlite3.connect(':memory:')
#con = sqlite3.connect('pdb.db')
c = con.cursor()
c.execute("CREATE TABLE chains (chain VARCHAR(10),unp VARCHAR(10),pdb VARCHAR(20),pdbstart INT,offset INT,species VARCHAR(20),PRIMARY KEY(chain))")
#c.execute("CREATE TABLE peptide (chain VARCHAR(10),seqres INT,aminoacid VARCHAR(5))")
c.execute("CREATE TABLE coords (chain VARCHAR(10),serialnum INT,seqres INT,aa3 VARCHAR(3), aa1 VARCHAR(1),x DOUBLE,y DOUBLE,z DOUBLE,PRIMARY KEY(chain,serialnum))")
con.commit()

# Pull the PDB ID from the filename
pdb_id = os.path.basename(sys.argv[1][0:-4])
DBNAME = sys.argv[2]

# Parse the PDB file
fin = open(sys.argv[1],'r')
species = read_species(fin)
unp_chain_map = read_dbref(fin,c,pdb_id,species)
#hain_pep_map = read_seqres(fin,c)
chainatom_coord_map = read_atom(fin,c)
fin.close()
# Desired information from PDB has now been extracted

# Average the 3D coordinates for each amino acid
c.execute("CREATE TABLE avgcoords AS SELECT chain,seqres,aa3,aa1,AVG(x) as x,AVG(y) as y,AVG(z) as z FROM coords GROUP BY chain,seqres")
con.commit()

# Output the derived data to tabular
c.execute("SELECT avgcoords.chain,chains.species,chains.unp,chains.pdb,(avgcoords.seqres+chains.offset) as seqres,aa3,aa1,x,y,z FROM chains,avgcoords WHERE avgcoords.chain=chains.chain AND avgcoords.seqres >= chains.pdbstart ORDER BY avgcoords.chain,avgcoords.seqres")
fout = open("%s.tab"%pdb_id,'w')
writer = csv.writer(fout,delimiter="\t")
writer.writerow(["chain","species","unp","pdbid","seqres","aa3","aa1","x","y","z"])
writer.writerows(c.fetchall())
fout.close()

# Pull the genomic locations for each chain
with open("GenomicCoords.tab","w") as fout:
	fout.write("transcript\tseqres\taa1\tstart\tend\tchr\tstrand\n")
with open("PDBTranscript.tab","w") as fout:
	fout.write("pdb\tchain\ttranscript\n")
c.execute("SELECT unp,chain FROM chains")
unp_chains = c.fetchall()
for unp_chain in unp_chains:
	unp = str(unp_chain[0])
	chain = str(unp_chain[1])
	print "Loading -> pdb: %s, chain: %s, unp: %s, species: %s"%(pdb_id,chain,unp,species)
	os.system("./protein_to_genomic.pl %s %s %s %s"%(pdb_id,chain,unp,species))
con.close()

# Load all data into Bushlab database
con = MySQLdb.connect(host='loki.mc.vanderbilt.edu',user='sivleyrm',passwd='19jcpkxw88',db=DBNAME)
c = con.cursor()
query = "LOAD DATA LOCAL INFILE 'GenomicCoords.tab' INTO TABLE GenomicCoords FIELDS TERMINATED BY '\t' IGNORE 1 LINES"
c.execute(query)
query = "LOAD DATA LOCAL INFILE 'PDBTranscript.tab' INTO TABLE PDBTranscript FIELDS TERMINATED BY '\t' IGNORE 1 LINES"
c.execute(query)
query = "LOAD DATA LOCAL INFILE '%s.tab' INTO TABLE PDBCoords FIELDS TERMINATED BY '\t' IGNORE 1 LINES (chain,@dummy,@dummy,pdbid,seqres,aa3,aa1,x,y,z)"
c.execute(query%pdb_id)
query = "LOAD DATA LOCAL INFILE '%s.tab' INTO TABLE PDBInfo FIELDS TERMINATED BY '\t' IGNORE 1 LINES (chain,species,unp,pdbid,@dummy,@dummy,@dummy,@dummy,@dummy,@dummy)"
c.execute(query%pdb_id)
# Sanitation query
query = "CALL sanitize_transcripts();"
c.execute(query)
# Pull query?
query = "CALL build_GenomePDB();"
c.execute(query)
con.close()