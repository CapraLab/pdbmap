#!/usr/bin/python27

# This script takes a configuration file as input, specifying the PDB archive,
# the database information in which to store the PDB->Genome mappings,
# and booleans for whether to create a new database or to identify human
# homologues for non-human protein structures.

import argparse,ConfigParser
import os,sys,csv,time
import sqlite3,MySQLdb
from warnings import filterwarnings	# Disable MySQL warnings
filterwarnings('ignore',category=MySQLdb.Warning)
filterwarnings('ignore',"Unknown table.*")

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
	"create_new_db" : False,
	"force" : False
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
parser.add_argument("-f", "--force", action='store_true', help="Force configuration. No safety checks.")


def main():
	"""Map each PDB to its genomic coordinates and upload to the specified database"""

	# Profiler
	t0 = time.time()

	# If a new database needs to be created, create one.
	if args.create_new_db:
		create_new_db(args.dbhost,args.dbuser,args.dbpass,args.dbname)

	# Identify all PDB structures in pdb_dir
	if os.path.isdir(args.pdb_dir):
		pdb_files = ['%s/%s'%(args.pdb_dir,pdb_file) for pdb_file in os.listdir(args.pdb_dir)]
	else:
		pdb_files = [args.pdb_dir]
	pdbs   = dict((os.path.basename(pdb_file[0:-4]),pdb_file) for pdb_file in pdb_files)

	

	# Load each PDB structure
	pdb_count = 0
	print "%d PDB file(s) found."%len(pdb_files)
	for pdb_id,pdb_file in pdbs.iteritems():
		if pdb_in_db(pdb_id,args.dbhost,args.dbuser,args.dbpass,args.dbname):
			print "PDB %s already included in database. Skipping..."%pdb_id
		else:
			print "Processing PDB %s..."%pdb_id
			pdb_count += 1
			try:
				load_pdb(pdb_id,pdb_file)
			except Exception as e:
				sys.stderr.write("PDB %s could not be processed.\n"%pdb_id)
				sys.stderr.write("%s\n"%e)
				sys.stderr.write("Skipping...\n")
		sys.stdout.flush()	# Flush all output for this PDB file
	
	# Profiler
	t_elapsed = time.time() - t0
	t_average = t_elapsed / pdb_count
	print "Number of PDBs processed: %d"%pdb_count
	print "Total execution time: %f"%t_elapsed
	print "Average execution time per PDB: %f"%t_average
	print "...Complete"

def sqlite_init():
	"""Initializes a temporary SQLite instance for local PDB->Genome Mapping"""
	con = sqlite3.connect(':memory:')
	c = con.cursor()
	c.execute("CREATE TABLE chains (chain VARCHAR(1),unp VARCHAR(10),pdb VARCHAR(20),pdbstart INT,offset INT,species VARCHAR(20),PRIMARY KEY(chain))")
	c.execute("CREATE TABLE coords (chain VARCHAR(1),serialnum INT,seqres INT,aa3 VARCHAR(3), aa1 VARCHAR(1),x DOUBLE,y DOUBLE,z DOUBLE,PRIMARY KEY(chain,serialnum))")
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
		if not species:
			print "PDB contained no species information. Skipping..."
			return
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
		if species.lower() in ['human','homo sapien','homo sapiens']:
			print "\tLoading -> pdb: %s, chain: %s, unp: %s, species: %s"%(pdb_id,chain,unp,species)
			os.system("./protein_to_genomic.pl %s %s %s %s"%(pdb_id,chain,unp,species))
		elif not args.disable_human_homologue:
			ung = unp2ung(unp)
			if not ung:
				sys.stderr.write("No UniGene entry found -> pdb: %s, chain: %s, unp: %s, species: %s"%(pdb_id,chain,unp,species))
				return
			else:
				print "\tSearching for human homologues -> pdb: %s, chain: %s, unp: %s, ung: %s, species: %s"%(pdb_id,chain,unp,ung,species)
				os.system("./unigene_to_homologue_genomic.pl %s %s %s %s"%(pdb_id,chain,ung,species))
		else:
			print "\tSpecies is non-human and human homologue mapping is disabled. Skipping..."

	# Upload to the database
	print "\tUploading data to database..."
	publish_data(pdb_id,args.dbhost,args.dbuser,args.dbpass,args.dbname)


def read_species(fin):
	"""Parses the species field from a PDB file"""
	flag = 0
	common = ''
	scientific = ''
	for line in fin:
		if line[0:6] == "SOURCE":
			flag = 1
			if line[11:26] == "ORGANISM_COMMON":
				common = line[28:].split(';')[0].strip()
				if common in species_map:
					common = species_map[common]
			if line[11:30] == "ORGANISM_SCIENTIFIC":
				scientific = line[32:].split(';')[0].strip()
				if scientific in species_map:
					scientific = species_map[scientific]
		elif flag:
			if common: # Return common if found
				return common
			elif scientific: # Return scientific if not
				return scientific
			else: # Empty return if no species information
				return


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
	except MySQLdb.Error as e:
		print "There was an error connecting to the database.\n%s"%e
		sys.exit(1)
	query = ["LOAD DATA LOCAL INFILE 'GenomicCoords.tab' INTO TABLE GenomicCoords FIELDS TERMINATED BY '\t' IGNORE 1 LINES"]
	query.append("LOAD DATA LOCAL INFILE 'PDBTranscript.tab' INTO TABLE PDBTranscript FIELDS TERMINATED BY '\t' IGNORE 1 LINES")
	query.append("LOAD DATA LOCAL INFILE '%s.tab' INTO TABLE PDBCoords FIELDS TERMINATED BY '\t' IGNORE 1 LINES (chain,@dummy,@dummy,pdbid,seqres,aa3,aa1,x,y,z)"%pdb_id)
	query.append("LOAD DATA LOCAL INFILE '%s.tab' INTO TABLE PDBInfo FIELDS TERMINATED BY '\t' IGNORE 1 LINES (chain,species,unp,pdbid,@dummy,@dummy,@dummy,@dummy,@dummy,@dummy)"%pdb_id)
	query.append("CALL %s.sanitize_transcripts(%s);"%(dbname,pdb_id))
	query.append("CALL %s.update_GenomePDB('%s');"%(dbname,pdb_id))
	end_of_query = ['...','.']
	for q in query:
		tq0  = time.time()
		print("\t%s%s - "%(q[:20],end_of_query[int(len(q)>20))]),
		c.execute(q)
		tqel = time.time()-tq0
		print("%2.2fs"%tqel)
	con.close()	# Close the remote MySQL connection
	os.system('rm -f GenomicCoords.tab')
	os.system('rm -f PDBTranscript.tab')
	os.system('rm -f %s.tab'%pdb_id)

def pdb_in_db(pdb_id,dbhost,dbuser,dbpass,dbname):
	try:
		con = MySQLdb.connect(host=dbhost,user=dbuser,passwd=dbpass,db=dbname)
		c = con.cursor()
	except MySQLdb.Error as e:
		print "There was an error connecting to the database.\n%s"%e
		sys.exit(1)
	c.execute("SELECT * FROM PDBInfo WHERE pdbid=%s",pdb_id)
	res = c.fetchone()
	con.close()
	return res

def create_new_db(dbhost,dbuser,dbpass,dbname):
	try:
		con = MySQLdb.connect(host=dbhost,user=dbuser,passwd=dbpass)
		c = con.cursor()
	except MySQLdb.Error as e:
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
	with open('create_proc_update_GenomePDB.sql','r') as fin:
		query.append(fin.read()%format_dict)
	with open('create_proc_sanitize_transcripts.sql','r') as fin:
		query.append(fin.read()%format_dict)
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


if __name__=='__main__':
	args = parser.parse_args(remaining_argv)
	args.create_new_db = bool(args.create_new_db)
	args.disable_human_homologue = bool(args.disable_human_homologue)
	args.force = bool(args.force)
	if args.create_new_db and not args.force:
		print "You have opted to create a new database."
		print "This will overwrite any existing database by the same name."
		if raw_input("Are you sure you want to do this? (y/n):") == 'n':
			print "Aborting..."
			sys.exit(0)
	main()