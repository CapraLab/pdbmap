#!/usr/bin/python27

# This script takes a configuration file as input, specifying the PDB archive,
# the database information in which to store the PDB->Genome mappings,
# and booleans for whether to create a new database or to identify human
# homologues for non-human protein structures.

import argparse,ConfigParser
import os,sys,csv,sqlite3,MySQLdb
from warnings import filterwarnings	# Disable MySQL warnings
filterwarnings('ignore',category=MySQLdb.Warning)

# Setup the Config File Parser
conf_parser = argparse.ArgumentParser(add_help=False)
conf_parser.add_argument("-c", "--conf_file",
help="Specify config file", metavar="FILE")
args, remaining_argv = conf_parser.parse_known_args()
defaults = {
	"dbhost" : "localhost",
	"dbname" : "pdbmap",
	"dbuser" : "user",
	"dbpass" : "pass",
	"pdb_dir" : "./pdb",
	"disable_human_homologue" : False,
	"create_new_db" : False
	}
if args.conf_file:
	config = ConfigParser.SafeConfigParser()
	config.read([args.conf_file])
	defaults = dict(config.items("Genome_PDB_Mapper"))

# Setup the Command Line Parser
parser = argparse.ArgumentParser(parents=[conf_parser],description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
defaults['disable_human_homologue'] = ('True' == defaults['disable_human_homologue'])
defaults['create_new_db'] = ('True' == defaults['create_new_db'])
parser.set_defaults(**defaults)
parser.add_argument("--dbhost", help="Database hostname")
parser.add_argument("--dbuser", help="Database username")
parser.add_argument("--dbpass", help="Database password")
parser.add_argument("--dbname", help="Database name")
parser.add_argument("--pdb_dir", help="Directory containing PDB files")
parser.add_argument("--disable_human_homologue", action='store_true', help="Disable mapping from non-human PDBs to human homologues")
parser.add_argument("--create_new_db", action='store_true', help="Create a new database prior to uploading PDB data")


def main():
	"""Map each PDB to its genomic coordinates and upload to the specified database"""

	# If a new database needs to be created, create one.
	if args.create_new_db:
		create_new_database(args.dbhost,args.dbuser,args.dbpass,args.dbname)

	# Identify all PDB structures in pdb_dir
	if os.path.isdir(args.pdb_dir):
		pdb_files = ['%s/%s'%(args.pdb_dir,pdb_file) for pdb_file in os.listdir(args.pdb_dir)]
	else:
		pdb_files = [args.pdb_dir]
	pdbs   = dict((os.path.basename(pdb_file[0:-4]),pdb_file) for pdb_file in pdb_files)

	

	# Load each PDB structure
	print "%d PDB file(s) found."%len(pdb_files)
	for pdb_id,pdb_file in pdbs.iteritems():
		print "Processing PDB %s..."%pdb_id
		load_pdb(pdb_id,pdb_file)
	

def sqlite_init():
	"""Initializes a temporary SQLite instance for local PDB->Genome Mapping"""
	con = sqlite3.connect(':memory:')
	c = con.cursor()
	c.execute("CREATE TABLE chains (chain VARCHAR(10),unp VARCHAR(10),pdb VARCHAR(20),pdbstart INT,offset INT,species VARCHAR(20),PRIMARY KEY(chain))")
	c.execute("CREATE TABLE coords (chain VARCHAR(10),serialnum INT,seqres INT,aa3 VARCHAR(3), aa1 VARCHAR(1),x DOUBLE,y DOUBLE,z DOUBLE,PRIMARY KEY(chain,serialnum))")
	con.commit()
	return c,con

def load_pdb(pdb_id,pdb_file):
	"""Parse a PDB file and store the important information into SQLite"""

	# Initialize a temporary SQLite database
	c,con = sqlite_init()

	# Load all information from the PDB file
	print "\tParsing info from %s..."%pdb_id
	with open(pdb_file,'r') as fin:
		species = read_species(fin)
		read_dbref(fin,c,pdb_id,species)
		read_atom(fin,c)
		con.commit()

	# Write data to tabular for MySQL upload
	write_pdb_data(pdb_id,c,con)
	c.execute("SELECT unp,chain FROM chains")
	unp_chains = c.fetchall()
	con.close()	# Close the local SQLite connection

	# Create tabulars for the genomic data
	with open("GenomicCoords.tab","w") as fout:
		fout.write("transcript\tseqres\taa1\tstart\tend\tchr\tstrand\n")
	with open("PDBTranscript.tab","w") as fout:
		fout.write("pdb\tchain\ttranscript\thomologue\n")
	for unp_chain in unp_chains:
		unp = str(unp_chain[0])
		chain = str(unp_chain[1])
		if species.lower() == 'human':
			print "\tLoading -> pdb: %s, chain: %s, unp: %s, species: %s"%(pdb_id,chain,unp,species)
			os.system("./protein_to_genomic.pl %s %s %s %s"%(pdb_id,chain,unp,species))
		elif not args.disable_human_homologue:
			ung = unp2ung(unp)
			print "\tLoading any human homologues -> pdb: %s, chain: %s, unp: %s, ung: %s, species: %s"%(pdb_id,chain,unp,ung,species)
			os.system("./unigene_to_homologue_genomic.pl %s %s %s %s"%(pdb_id,chain,ung,species))
		else:
			print "\tSpecies is non-human and human homologue mappning is disabled. Skipping..."

	# Upload to the database
	publish_data(pdb_id,args.dbhost,args.dbuser,args.dbpass,args.dbname)



def read_species(fin):
	"""Parses the species field from a PDB file"""
	for line in fin:
		if line[0:6] == "SOURCE" and line[11:26] == "ORGANISM_COMMON":
			species = line[28:].split(';')[0]
			if species in species_map:
				species = species_map[species]
			return species


def read_dbref(fin,c,pdb_id,species):
	"""Parses the DBREF fields from a PDB file"""
	flag = 0
	for line in fin:
		if line[0:5] == "DBREF":
			flag = 1
			if line[26:32].strip() == "UNP":
				chain = line[12]
				unp_id = line[33:41].strip()
				pdbstart = int(line[14:18])
				dbstart = int(line[55:60])
				offset = dbstart - pdbstart
				c.execute("INSERT INTO chains VALUES (?,?,?,?,?,?)",(chain,unp_id,pdb_id,pdbstart,offset,species))
		elif flag:
			return

def read_seqres(fin,c):
	"""Parses the SEQRES fields from a PDB file"""
	flag = 0
	index = 0
	for line in fin:
		if line[0:6] == "SEQRES":
			flag = 1
			chain = line[11]
			pep_subseq = line[19:70].split()
			for aa in pep_subseq:
				c.execute("INSERT INTO peptide VALUES (?,?,?)",(chain,index,aa))
				index += 1
		elif flag:
			return 

def read_atom(fin,c):
	"""Parses the ATOM fields from a PDB file"""
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

def write_pdb_data(pdb_id,c,con):
	"""Averages the 3D coordinates and outputs the PDB data to tabular for MySQL upload"""
	c.execute("CREATE TABLE avgcoords AS SELECT chain,seqres,aa3,aa1,AVG(x) as x,AVG(y) as y,AVG(z) as z FROM coords GROUP BY chain,seqres")
	con.commit()
	c.execute("SELECT avgcoords.chain,chains.species,chains.unp,chains.pdb,(avgcoords.seqres+chains.offset) as seqres,aa3,aa1,x,y,z FROM chains,avgcoords WHERE avgcoords.chain=chains.chain AND avgcoords.seqres >= chains.pdbstart ORDER BY avgcoords.chain,avgcoords.seqres")
	fout = open("%s.tab"%pdb_id,'w')
	writer = csv.writer(fout,delimiter="\t")
	writer.writerow(["chain","species","unp","pdbid","seqres","aa3","aa1","x","y","z"])
	writer.writerows(c.fetchall())
	fout.close()

def unp2ung(unp):
	import urllib,urllib2
	url = 'http://www.uniprot.org/mapping/'
	params = {
	'from':'ACC+ID',
	'to':'UNIGENE_ID',
	'format':'tab',
	'query': '%s'%unp
	}
	data = urllib.urlencode(params)
	request = urllib2.Request(url, data)
	contact = "mike.sivley@vanderbilt.edu" # Please set your email address here to help us debug in case of problems.
	request.add_header('User-Agent', 'Python %s' % contact)
	response = urllib2.urlopen(request)
	page = response.read(200000)
	return page.split()[-1]

def publish_data(pdb_id,dbhost,dbuser,dbpass,dbname):
	try:
		con = MySQLdb.connect(host=dbhost,user=dbuser,passwd=dbpass,db=dbname)
		c = con.cursor()
	except Exception as e:
		print "There was an error connecting to the database.\n%s"%e
		sys.exit(1)
	query = ["LOAD DATA LOCAL INFILE 'GenomicCoords.tab' INTO TABLE GenomicCoords FIELDS TERMINATED BY '\t' IGNORE 1 LINES"]
	query.append("LOAD DATA LOCAL INFILE 'PDBTranscript.tab' INTO TABLE PDBTranscript FIELDS TERMINATED BY '\t' IGNORE 1 LINES")
	query.append("LOAD DATA LOCAL INFILE '%s.tab' INTO TABLE PDBCoords FIELDS TERMINATED BY '\t' IGNORE 1 LINES (chain,@dummy,@dummy,pdbid,seqres,aa3,aa1,x,y,z)"%pdb_id)
	query.append("LOAD DATA LOCAL INFILE '%s.tab' INTO TABLE PDBInfo FIELDS TERMINATED BY '\t' IGNORE 1 LINES (chain,species,unp,pdbid,@dummy,@dummy,@dummy,@dummy,@dummy,@dummy)"%pdb_id)
	query.append("CALL %s.sanitize_transcripts();"%dbname)
	query.append("CALL %s.build_GenomePDB();"%dbname)
	for q in query:
		c.execute(q)
	con.close()	# Close the remote MySQL connection
	os.system('rm %s.tab'%pdb_id)

def create_new_database(dbhost,dbuser,dbpass,dbname):
	try:
		con = MySQLdb.connect(host=dbhost,user=dbuser,passwd=dbpass)
		c = con.cursor()
	except Exception as e:
		print "There was an error connecting to the database.\n%s"%e
		sys.exit(1)
	query = ["DROP DATABASE IF EXISTS %s"%dbname]
	query.append("CREATE DATABASE IF NOT EXISTS %s"%dbname)
	query.append("USE %s"%dbname)
	query.append("CREATE TABLE %s.GenomicCoords (transcript VARCHAR(20),seqres INT,aa1 VARCHAR(1),start BIGINT,end BIGINT,chr INT,strand INT,PRIMARY KEY(start,end,chr))"%dbname)
	query.append("CREATE TABLE %s.PDBCoords (pdbid VARCHAR(20),chain VARCHAR(1),seqres INT,aa3 VARCHAR(3),aa1 VARCHAR(1),x DOUBLE,y DOUBLE,z DOUBLE,PRIMARY KEY(pdbid,chain,seqres))"%dbname)
	query.append("CREATE TABLE %s.PDBInfo (pdbid VARCHAR(20),species VARCHAR(20),chain VARCHAR(1),unp VARCHAR(20),PRIMARY KEY(pdbid,chain))"%dbname)
	query.append("CREATE TABLE %s.PDBTranscript (pdbid VARCHAR(20),chain VARCHAR(1),transcript VARCHAR(20),homologue BOOLEAN,PRIMARY KEY(pdbid,chain))"%dbname)
	format_dict = {'dbuser':dbuser,'dbhost':dbhost}
	query.append(build_genomePDB_proc%format_dict)
	query.append(sanitize_transcripts_proc%format_dict)
	for q in query:
		c.execute(q)

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

build_genomePDB_proc ="""
CREATE DEFINER=`%(dbuser)s`@`%(dbhost)s` PROCEDURE `build_GenomePDB`()
BEGIN

DROP TABLE IF EXISTS GenomePDB;
DROP TABLE IF EXISTS PDB;

CREATE TABLE PDB AS
(SELECT PDBInfo.pdbid,PDBInfo.species,PDBInfo.chain,PDBInfo.unp,PDBCoords.seqres,PDBCoords.aa1,PDBCoords.x,PDBCoords.y,PDBCoords.z,PDBTranscript.transcript
FROM PDBInfo,PDBCoords,PDBTranscript
WHERE PDBInfo.pdbid=PDBCoords.pdbid AND PDBInfo.chain=PDBCoords.chain AND PDBInfo.pdbid=PDBTranscript.pdbid AND PDBInfo.chain=PDBTranscript.chain AND PDBCoords.pdbid=PDBTranscript.pdbid AND PDBCoords.chain=PDBTranscript.chain);

CREATE TABLE GenomePDB AS
SELECT PDB.pdbid,PDB.species,PDB.chain,PDB.unp,PDB.transcript,PDB.seqres,PDB.aa1,PDB.x,PDB.y,PDB.z,GenomicCoords.start,GenomicCoords.end,GenomicCoords.chr,GenomicCoords.strand
FROM PDB,GenomicCoords
WHERE PDB.transcript=GenomicCoords.transcript AND PDB.seqres=GenomicCoords.seqres
GROUP BY PDB.pdbid,PDB.chain,PDB.seqres
ORDER BY PDB.pdbid,PDB.chain,PDB.seqres;

DROP TABLE PDB;
END
"""

sanitize_transcripts_proc ="""
CREATE DEFINER=`%(dbuser)s`@`%(dbhost)s` PROCEDURE `sanitize_transcripts`()
BEGIN

DROP TABLE IF EXISTS minTrans;
DROP TABLE IF EXISTS maxTrans;

# Genomic transcripts that match the minimum seqres for their PDB Sequence
CREATE TABLE minTrans AS (SELECT PDBTranscript.transcript FROM PDBTranscript,GenomicCoords as GC WHERE
		# PDB seqres must match the GC seqres (no out of range)
		GC.seqres=(SELECT seqres FROM PDBCoords WHERE PDBCoords.pdbid=PDBTranscript.pdbid AND PDBCoords.chain=PDBTranscript.chain ORDER BY seqres LIMIT 1)
		AND
		# PDB amino acid must match the GC amino acid (no mismatched peptides)
		GC.aa1=(SELECT aa1 FROM PDBCoords WHERE PDBCoords.pdbid=PDBTranscript.pdbid AND PDBCoords.chain=PDBTranscript.chain ORDER BY seqres LIMIT 1)
	);

# Genomic transcripts that match the maximum seqres for their PDB Sequence
CREATE TABLE maxTrans AS (SELECT PDBTranscript.transcript FROM PDBTranscript,GenomicCoords as GC WHERE
		# PDB seqres must match the GC seqres (no out of range)
		GC.seqres=(SELECT seqres FROM PDBCoords WHERE PDBCoords.pdbid=PDBTranscript.pdbid AND PDBCoords.chain=PDBTranscript.chain ORDER BY seqres DESC LIMIT 1)
		AND
		# PDB amino acid must match the GC amino acid (no mismatched peptides)
		GC.aa1=(SELECT aa1 FROM PDBCoords WHERE PDBCoords.pdbid=PDBTranscript.pdbid AND PDBCoords.chain=PDBTranscript.chain ORDER BY seqres DESC LIMIT 1)
	);

DELETE FROM PDBTranscript
WHERE
# Delete if the amino acids don't match at the minimum seqres
transcript NOT IN (SELECT transcript FROM minTrans)
OR
# Delete if the amino acids don't match at the maximum seqres
transcript NOT IN (SELECT transcript FROM maxTrans);
DROP TABLE minTrans;
DROP TABLE maxTrans;
END
"""

if __name__=='__main__':
	args = parser.parse_args(remaining_argv)
	args.create_new_db = bool(args.create_new_db)
	args.disable_human_homologue = bool(args.disable_human_homologue)
	main()