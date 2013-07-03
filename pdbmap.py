#!/usr/bin/python27

# This script takes a configuration file as input, specifying the PDB archive,
# the database information in which to store the PDB->Genome mappings,
# and booleans for whether to create a new database or to identify human
# homologues for non-human protein structures.

import argparse,ConfigParser
import os,sys,csv,time,subprocess,gzip
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
	pdbs   = dict((os.path.basename(pdb_file).split('.')[0][-4:].upper(),pdb_file) for pdb_file in pdb_files)

	

	# Load each PDB structure
	pdb_count = 0
	print "%d PDB file(s) found."%len(pdb_files)
	for pdb_id,pdb_file in pdbs.iteritems():
		if not pdb_in_db(pdb_id,args.dbhost,args.dbuser,args.dbpass,args.dbname):
			print "\nProcessing PDB %s..."%pdb_id
			pdb_count += 1
			try:
				load_pdb(pdb_id,pdb_file)
			except (KeyboardInterrupt,SystemExit):
				print("\nExiting...")
				sys.exit(0)
			except Exception as e:
				sys.stderr.write("\tPDB %s could not be processed.\n"%pdb_id)
				sys.stderr.write("\t%s\n"%e)
				sys.stderr.write("\tSkipping...\n")
		sys.stdout.flush()	# Flush all output for this PDB file
	
	# Profiler
	t_elapsed = time.time() - t0
	if pdb_count > 0:
		t_average = t_elapsed / pdb_count
	else:
		t_average = 0
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
	c.execute("CREATE TABLE seqadv (chain VARCHAR(1),seqres INT,aa3 VARCHAR(3),conflict VARCHAR(20),PRIMARY KEY(chain,seqres))")
	con.commit()
	return c,con

def load_pdb(pdb_id,pdb_file):
	"""Parse a PDB file and store the important information into SQLite"""

	# Initialize a temporary SQLite database
	c,con = sqlite_init()

	# Load all information from the PDB file
	# Use gzip if the files appear to be gzipped
	print "\tParsing info from %s..."%pdb_file
	if os.path.splitext(pdb_file)[-1] == '.gz':
		fin = gzip.open(pdb_file,'r')
	else:
		fin = open(pdb_file,'r')
	species = read_species(fin)
	if not species:
		sys.stderr.write("\tPDB contained no species information. Skipping...\n")
		return
	elif (not species.lower() in ['human','homo sapien','homo sapiens']) and args.disable_human_homologue:
		sys.stderr.write("\tSpecies is non-human and human homologue mapping is disabled. Skipping...\n")
		return

	experiment_type = None
	best_model      = None
	cur_model       = None
	for line in fin:
		field = line[0:6].strip()
		code  = line[7:10]
		if field == "REMARK" and code[0] == "2":
			print "REMARK Code: %s"%code
			if not experiment_type:
				experiment_type   = read_experiment_type(line,c)
			elif experiment_type == "NMR" and not best_model:
				best_model = read_best_model(line,c)
		if field == "MODEL":
			cur_model = int(line[10:14])
		elif field == "DBREF":
			read_dbref(line,c,pdb_id,species)
		elif field == "SEQADV":
			read_seqadv(line,c)
		elif field == "ATOM"   and cur_model == best_model:
			read_atom(line,c)
	print("Experiment type: %s"%experiment_type)
	if experiment_type == "NMR":
		print("Best Model: %d"%best_model)
	con.commit()
	fin.close()

	# Remove residues with unwanted conflicts
	print("\tRemoving SEQADV conflicts...")
	c.execute("SELECT chain,seqres FROM seqadv WHERE conflict='EXPRESSION TAG'")
	for chain,seqres in c.fetchall():
		print("\t\tSEQADV %s,%s removed"%(chain,seqres))
		print("\t\tConflict: EXPRESSION TAG")
		print("\t\t#-------------#")
		c.execute("DELETE FROM coords WHERE chain=? AND seqres=?",(chain,seqres))

	# Write data to tabular for MySQL upload
	write_pdb_data(pdb_id,c,con)

	# Pull local processing information
	c.execute("SELECT unp,chain FROM chains")
	unp_chains = c.fetchall()
	c.execute("SELECT chain,seqres FROM seqadv WHERE conflict='ENGINEERED MUTATION'")
	seqadv_protected = c.fetchall()
	con.close()	# Close the local SQLite connection

	# Create tabulars for the genomic data
	with open("GenomicCoords.tab","w") as fout:
		fout.write("transcript\tseqres\taa1\tstart\tend\tchr\tstrand\n")
	with open("PDBTranscript.tab","w") as fout:
		fout.write("pdb\tchain\ttranscript\thomologue\n")
	for unp_chain in unp_chains:
		unp = str(unp_chain[0])
		chain = str(unp_chain[1])
		try:
			if species.lower() in ['human','homo sapien','homo sapiens']:
				print "\tLoading -> pdb: %s, chain: %s, unp: %s, species: %s"%(pdb_id,chain,unp,species)
				exit_code = subprocess.call(["./protein_to_genomic.pl",pdb_id,chain,unp,species])
			else:
				ung = unp2ung(unp)
				if not ung:
					sys.stderr.write("No UniGene entry found -> pdb: %s, chain: %s, unp: %s, species: %s"%(pdb_id,chain,unp,species))
					return
				else:
					print "\tSearching for human homologues -> pdb: %s, chain: %s, unp: %s, ung: %s, species: %s"%(pdb_id,chain,unp,ung,species)
					exit_code = subprocess.call(["./unigene_to_homologue_genomic.pl",pdb_id,chain,ung,species])
		except KeyboardInterrupt:
			if raw_input("\nContinue to next PDB? (y/n):") == 'n':
				raise KeyboardInterrupt
			print("\tSkipping...")
			return
		if exit_code:
			sys.stderr.write("\tPerl script returned a non-zero exit status. Skipping...\n")
			return		

	# Upload to the database
	num_matches = sanitize_data(pdb_id,seqadv_protected)
	print("\t#----#\n\t%d transcripts found.\n\t#----#"%num_matches)
	publish_data(pdb_id,args.dbhost,args.dbuser,args.dbpass,args.dbname,num_matches)


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

def read_experiment_type(line,c):
	print "Checking for experiment type field, field:"
	print line[12:27]
	"""Parses a REMARK 2** field for experiment type"""
	if line[12:27] == "EXPERIMENT TYPE":
		return line[45:].strip()
	else:
		return None

def read_best_model(line,c):
	print "Checking for best conformer field, field:"
	print line[11:57]
	"""Parses a REMARK 210 field for best NMR model"""
	if line[11:57] == "BEST REPRESENTATIVE CONFORMER IN THIS ENSEMBLE":
		best_nmr = line[60:].strip()
		print "best NMR model: %s"%best_nmr
		if best_nmr == "NULL":
			best_nmr = 1
		return int(best_nmr)
	else:
		return None

def read_dbref(line,c,pdb_id,species):
	"""Parses a DBREF field from a PDB file"""
	if line[26:32].strip() == "UNP":
		chain = line[12]
		unp_id = line[33:41].strip()
		pdbstart = int(line[14:18])
		dbstart = int(line[55:60])
		offset = dbstart - pdbstart
		c.execute("INSERT INTO chains VALUES (?,?,?,?,?,?)",(chain,unp_id,pdb_id,pdbstart,offset,species))

def read_seqres(line,c):
	"""Parses a SEQRES field from a PDB file"""
	chain = line[11]
	c.execute("SELECT * FROM chains WHERE chain=?",chain)
	if c.fetchone():
		pep_subseq = line[19:70].split()
		for aa in pep_subseq:
			c.execute("INSERT INTO seqres VALUES (?,?,?)",(chain,read_seqres.index,aa))
			read_seqres.index += 1
read_seqres.index = 0

def read_seqadv(line,c):
	"""Parses a SEQADV field from a PDB file"""
	aa = line[12:15]
	chain = line[16]
	seqres = int(line[18:22])
	conflict = line[49:70].strip()
	c.execute("INSERT INTO seqadv VALUES (?,?,?,?)",(chain,seqres,aa,conflict))

def read_atom(line,c):
	"""Parses a ATOM field from a PDB file"""
	chain = line[21]
	c.execute("SELECT * FROM chains WHERE chain=?",chain)
	if c.fetchone():
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

def sanitize_data(pdb_id,seqadv_protected):
	print("\tSanitizing transcript maps...")
	with open('GenomicCoords.tab','r') as fin:
		reader = csv.reader(fin,delimiter='\t')
		gc = {(row[0],row[1]): row[2] for row in reader}
	with open('PDBTranscript.tab','r') as fin:
		pdb_t_header = fin.next()
		reader = csv.reader(fin,delimiter='\t')
		# store list of tuples (transcript,chain,homologue)
		pdb_t = [(row[2],row[1],row[3]) for row in reader]
	with open('%s.tab'%pdb_id,'r') as fin:
		reader = csv.reader(fin,delimiter='\t')
		pdb_c = {(row[0],row[4]): row[6] for row in reader}
	sanitize = []
	for transcript,chain,homologue in pdb_t:
		# Reduce by first key and remove
		test1 = [key[0] for key in gc if key[0]==transcript]
		test2 = [key[0] for key in pdb_c if key[0]==chain]
		gc_sub    = {int(key[1]):gc[key]    for key in gc    if key[0]==transcript}
		pdb_c_sub = {int(key[1]):pdb_c[key] for key in pdb_c if key[0]==chain}
		# Join by second key and compare values
		for seqres in pdb_c_sub:
			# If a protected SEQADV conflict, skip
			if (chain,seqres) in seqadv_protected:
				pass
			# If the transcript does not contain the seqres, discard
			elif seqres not in gc_sub:
				print("\t\tTranscript %s sanitized. Conflict SEQRES: %s"%(transcript,seqres))
				print("\t\tChain: %s"%chain)
				print("\t\tCause: SEQRES out of transcript range.")
				print("\t\tTranscript range: %d - %d"%(min(gc_sub.keys()),max(gc_sub.keys())))
				print("\t\t#-------------#")
				sanitize.append((transcript,chain,homologue))
				break
			# Or the amino acid is different at that seqres, discard
			elif gc_sub[seqres] != pdb_c_sub[seqres]:
				print("\t\tTranscript %s sanitized. Conflict SEQRES: %s"%(transcript,seqres))
				print("\t\tChain: %s"%chain)
				print("\t\tCause: Amino acid mismatch.")
				print("\t\tPDB.aa: %s, Transcript.aa: %s"%(pdb_c_sub[seqres],gc_sub[seqres]))
				print("\t\t#-------------#")
				sanitize.append((transcript,chain,homologue))
				break

	# Remove the transcript->chain map
	for transcript,chain,homologue in sanitize:
		pdb_t.remove((transcript,chain,homologue))
	# And write the sanitized mappings back to file
	with open('PDBTranscript.tab','w') as fout:
		writer = csv.writer(fout,delimiter='\t')
		for transcript,chain,homologue in pdb_t:
			writer.writerow([pdb_id,chain,transcript,homologue])

	# Return the number of remaining matches
	return len(pdb_t)

def publish_data(pdb_id,dbhost,dbuser,dbpass,dbname,num_matches):
	try:
		con = MySQLdb.connect(host=dbhost,user=dbuser,passwd=dbpass,db=dbname)
		c = con.cursor()
	except MySQLdb.Error as e:
		print "There was an error connecting to the database.\n%s"%e
		sys.exit(1)
	try:
		query = ["LOAD DATA LOCAL INFILE '%s.tab' INTO TABLE PDBInfo FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n' IGNORE 1 LINES (chain,species,unp,pdbid,@dummy,@dummy,@dummy,@dummy,@dummy,@dummy)"%pdb_id]
		query.append("UPDATE PDBInfo SET unp=RTRIM(unp);")
		# Only load the rest if there were any successful matches
		if num_matches > 0:
			query.append("LOAD DATA LOCAL INFILE 'GenomicCoords.tab' INTO TABLE GenomicCoords FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n' IGNORE 1 LINES")
			query.append("LOAD DATA LOCAL INFILE 'PDBTranscript.tab' INTO TABLE PDBTranscript FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n' IGNORE 1 LINES")
			query.append("UPDATE PDBTranscript SET transcript=RTRIM(transcript);")
			query.append("LOAD DATA LOCAL INFILE '%s.tab' INTO TABLE PDBCoords FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n' IGNORE 1 LINES (chain,@dummy,@dummy,pdbid,seqres,aa3,aa1,x,y,z)"%pdb_id)
			query.append("CALL %s.update_GenomePDB('%s');"%(dbname,pdb_id))
		end_of_query = ['.','...']
		print("\tUploading staging files to database...")
		for q in query:
			tq0  = time.time()
			sys.stdout.write("\t\t%s%s - "%(q[:35],end_of_query[int(len(q)>35)]))
			sys.stdout.flush()
			c.execute(q)
			tqel = time.time()-tq0
			print("%2.2fs"%tqel)
	except (KeyboardInterrupt,SystemExit):
		sys.stderr.write("\tUpload to MySQL was interrupted by user. Skipping...")
		raise KeyboardInterrupt
	except Exception as e:
		sys.stderr.write("\tUpload to MySQL failed: %s. Skipping..."%e)
	finally:
		con.close()	# Close the remote MySQL connection
		os.system('rm -f PDBTranscript.tab')
		os.system('rm -f %s.tab'%pdb_id)
		os.system('rm -f GenomicCoords.tab')

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
	query.append("CREATE TABLE %s.GenomicCoords (transcript VARCHAR(20),seqres INT,aa1 VARCHAR(1),start BIGINT,end BIGINT,chr INT,strand INT,PRIMARY KEY(transcript,seqres),KEY(transcript),KEY(start,end,chr))"%dbname)
	query.append("CREATE TABLE %s.PDBCoords (pdbid VARCHAR(20),chain VARCHAR(1),seqres INT,aa3 VARCHAR(3),aa1 VARCHAR(1),x DOUBLE,y DOUBLE,z DOUBLE,PRIMARY KEY(pdbid,chain,seqres))"%dbname)
	query.append("CREATE TABLE %s.PDBInfo (pdbid VARCHAR(20),species VARCHAR(20),chain VARCHAR(1),unp VARCHAR(20),PRIMARY KEY(pdbid,chain))"%dbname)
	query.append("CREATE TABLE %s.PDBTranscript (pdbid VARCHAR(20),chain VARCHAR(1),transcript VARCHAR(20),homologue BOOLEAN,PRIMARY KEY(pdbid,chain),KEY(transcript))"%dbname)
	query.append("""CREATE TABLE `%s`.`GenomePDB` (
  		`pdbid` VARCHAR(20) NOT NULL default '',
  		`species` VARCHAR(20) default NULL,
  		`chain` VARCHAR(1) NOT NULL default '',
  		`unp` VARCHAR(20) default NULL,
  		`transcript` VARCHAR(20) default NULL,
  		`seqres` INT(11) NOT NULL default '0',
  		`aa1` VARCHAR(1) default NULL,
  		`x` DOUBLE default NULL,
  		`y` DOUBLE default NULL,
  		`z` DOUBLE default NULL,
  		`start` BIGINT(20) default NULL,
  		`end` BIGINT(20) default NULL,
  		`chr` INT(11) default NULL,
  		`strand` INT(11) default NULL,
  		PRIMARY KEY  (`pdbid`,`chain`,`seqres`),
  		KEY `pdbid` (`pdbid`),
  		KEY `pos` (`transcript`,`start`,`end`,`chr`)
		)"""%dbname)
	format_dict = {'dbuser':dbuser,'dbhost':dbhost}
	with open('create_proc_update_GenomePDB.sql','r') as fin:
		query.append(fin.read()%format_dict)
	for q in query:
		c.execute(q)
	con.close()

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