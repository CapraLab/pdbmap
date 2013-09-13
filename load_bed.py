#!/usr/bin/python27

import sys,os
import MySQLdb
from warnings import filterwarnings	# Disable MySQL warnings
filterwarnings('ignore',category=MySQLdb.Warning)
filterwarnings('ignore',"Unknown table.*")

# Command Line Arguments:
#	BED File
#	GenomePDB Map File
#	Intersection directory
#	[Optional] Data Name
if len(sys.argv) < 2:
	print "load_bed.py BED MAP intersection_dir [dat_name]"
	sys.exit(0);
bed_file = sys.argv[1]
map_file = sys.argv[2]
intersect_dir = sys.argv[3]
dat_name = "".join(os.path.basename(bed_file).split('.')[:-1])
if len(sys.argv) > 4:
	dat_name = sys.argv[4].replace('.','_')

# Filename Definitions
bed_temp         = bed_file+".temp"
intersect_file   = "%s/genome_pdb_%s.intersect"%(intersect_dir,dat_name)
nointersect_file = "%s/genome_pdb_%s.nointersect"%(intersect_dir,dat_name)

# Convert Chromosome 23+ Encoding
os.system("sed s/^23/X/g %s | sed s/^24/Y/g | sed s/^25/PAR/g | sed s/^26/MT/g > %s"%(bed_file,bed_temp))
os.system("mv %s %s"%(bed_temp,bed_file))

# Find intersections
os.system("/usr/analysis/bin/intersectBed -a %s -b %s -wb > %s"%(map_file,bed_file,intersect_file))

# Database Connection
try:
	con = MySQLdb.connect(host='gwar-dev',user='mike',passwd='cheezburger')
	c = con.cursor()
except MySQLdb.Error as e:
	print "There was an error connecting to the database.\n%s"%e
	sys.exit(1)

# Drop/Create/Load Intersection Table
c.execute("use pdbmap_v5;")
query = ["DROP TABLE IF EXISTS Intersect_Variants_%s;"%dat_name]
query.append("""CREATE TABLE Intersect_Variants_%s (
  `var_chr` varchar(10) NOT NULL default '',
  `var_start` int(11) NOT NULL default '0',
  `var_end` int(11) NOT NULL default '0',
  `var_name` varchar(50) default NULL,
  `chr` varchar(10) NOT NULL default '',
  `start` int(11) NOT NULL default '0',
  `end` int(11) NOT NULL default '0',
  `pdbid` varchar(20) NOT NULL default '',
  `chain` varchar(1) NOT NULL default '',
  `seqres` int(11) default NULL,
  `unp` varchar(20) default NULL,
  `transcript` varchar(20) NOT NULL default '',
  `gene` varchar(20) NOT NULL default '',
  `aa1` varchar(1) default NULL,
  `x` double default NULL,
  `y` double default NULL,
  `z` double default NULL,
  `species` varchar(20) default NULL,
  `strand`INT default 0,
  PRIMARY KEY  (`pdbid`,`chain`,`transcript`,`chr`,`start`,`end`,`var_chr`,`var_start`,`var_end`,`var_name`),
  KEY `var_chr` (`var_chr`,`var_start`,`var_end`),
  KEY `var_name` (`var_name`),
  KEY `bp` (`chr`,`start`,`end`),
  KEY `pdbid` (`pdbid`,`chain`,`seqres`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
"""%dat_name)
query.append("""LOAD DATA LOCAL INFILE '%s' INTO TABLE Intersect_Variants_%s
			(chr, start, end, strand, transcript, pdbid, chain, unp, seqres, 
			aa1, x, y, z, species, gene, var_chr, var_start, var_end, var_name)"""%(intersect_file,dat_name))

# Execute intersection queries
for q in query:
	c.execute(q)
con.close()	# Non-intersection pipeline will cause connection timeout

# Find non-intersections
os.system("/usr/analysis/bin/intersectBed -v -b %s -a %s > %s"%(map_file,bed_file,nointersect_file))

# Annotate the non-intersecting ranges with Ensembl Transcripts and UniProt IDs
os.system("./get_protein.py %s > %s.temp"%(nointersect_file,nointersect_file))
os.system("mv %s.temp %s"%(nointersect_file,nointersect_file))

# Database Connection
try:
	con = MySQLdb.connect(host='gwar-dev',user='mike',passwd='cheezburger')
	c = con.cursor()
except MySQLdb.Error as e:
	print "There was an error connecting to the database.\n%s"%e
	sys.exit(1)

# Drop/Create/Load Non-Intersection Table
c.execute("use pdbmap_v5;")
query = ["DROP TABLE IF EXISTS Non_Intersect_Variants_%s;"%dat_name]
query.append("""CREATE TABLE Non_Intersect_Variants_%s (
  `var_chr` varchar(10) NOT NULL default '',
  `var_start` int(11) NOT NULL default '0',
  `var_end` int(11) NOT NULL default '0',
  `var_name` varchar(50) default NULL,
  `transcript` varchar(50) default NULL,
  `unp` varchar(20) default NULL,
  PRIMARY KEY  (`var_chr`,`var_start`,`var_end`,`transcript`),
  KEY `var_name` (`var_name`),
  KEY `uniprot` (`unp`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
"""%dat_name)
query.append("""LOAD DATA LOCAL INFILE '%s' INTO TABLE
			Non_Intersect_Variants_%s (var_chr, var_start, var_end, var_name,
			transcript, unp)"""%(nointersect_file,dat_name))

# Execute non-intersection queries
for q in query:
	c.execute(q)
con.close()
