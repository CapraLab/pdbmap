#!/usr/bin/python27

import sys,os
import MySQLdb
from warnings import filterwarnings	# Disable MySQL warnings
filterwarnings('ignore',category=MySQLdb.Warning)
filterwarnings('ignore',"Unknown table.*")

# Command Line Arguments:
#	Variant .bed file or .vcf file
#	PDBMap .bed file
#	Intersection directory
#	[Optional] variant source name
if len(sys.argv) < 2:
  print "load_bed.py BED MAP intersection_dir [dat_name]"
  print "Ensure that PDBMAP file is 0-indexed, exclusive end"
  print "  to comply with UCSC BED format."
  sys.exit(0);
var_file = sys.argv[1]
pdbmap_file = sys.argv[2]
intersect_dir = sys.argv[3]
dat_name = sys.argv[4].replace('.','_')
append = False
if len(sys.argv) > 4:
  append = bool(sys.argv[5])

# Filename Definitions
bed_temp         = var_file+".temp"
intersect_file   = "%s/genome_pdb_%s.intersect"%(intersect_dir,dat_name)
intersect_file = intersect_file.replace('-','_')
nointersect_file = "%s/genome_pdb_%s.nointersect"%(intersect_dir,dat_name)
nointersect_file = nointersect_file.replace('-','_')

# Convert Chromosome 23+ Encoding
# UPDATE: Chromosome is now stored as chr1, chr12, chrX, etc.
#       : No conversion is necessary.
#os.system("sed s/^23/X/g %s | sed s/^24/Y/g | sed s/^25/PAR/g | sed s/^26/MT/g > %s"%(var_file,bed_temp))
#os.system("mv %s %s"%(bed_temp,var_file))

# Find intersections
ext = var_file.split('.')[-1].lower()
if ext == 'gz':
  print "Detected gzipped file."
  ext = var_file.split('.')[-2].lower()
if ext == 'bed':
  print "Detected file type: BED"
  os.system("/usr/analysis/bin/intersectBed -a %s -b %s -wa -wb | cut -f 1-21 > %s"%(pdbmap_file,var_file,intersect_file))
  # Adjust the start,end,var_start,var_end positions to 1-indexing
  cmd  = """awk -F"\t" -v OFS="\t" """
  cmd += """'{print $1,$2+1,$3+1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,"""
  cmd += """$16,$17,$18,$19+1,$20+1,$21}' """
  cmd += """%s > %sFIX"""%(intersect_file,intersect_file)
  os.system(cmd)
elif ext == 'vcf':
  print "Detected file type: VCF"
  os.system("/usr/analysis/bin/intersectBed -a %s -b %s -wb | cut -f 1-20 > %s"%(pdbmap_file,var_file,intersect_file))
  # Adjust the start,end positions to 1-indexing
  cmd  = """awk -F"\t" -v OFS="\t" """
  cmd += """'{print $1,$2+1,$3+1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,"""
  cmd += """$16,$17,$18,$19,$20}' """
  cmd += """%s > %sFIX"""%(intersect_file,intersect_file)
  os.system(cmd)
os.system("""mv -f %sFIX %s"""%(intersect_file,intersect_file))

# Database Connection
try:
	con = MySQLdb.connect(host='gwar-dev',user='mike',passwd='cheezburger')
	c = con.cursor()
except MySQLdb.Error as e:
	print "There was an error connecting to the database.\n%s"%e
	sys.exit(1)

# Drop/Create/Load Intersection Table
c.execute("use pdbmap_v7;")
if not append:
  query = ["DROP TABLE IF EXISTS Intersect_Variants_%s;"%dat_name]
  query.append("""CREATE TABLE Intersect_Variants_%s (
    `var_chr` varchar(10) NOT NULL default '',
    `var_start` int(11) NOT NULL default '0',
    `var_end` int(11) NOT NULL default '0',
    `var_name` varchar(50) default NULL,
    `chr` varchar(10) NOT NULL default '',
    `start` int(11) NOT NULL default '0',
    `end` int(11) NOT NULL default '0',
    `strand`INT default 0,
    `gene` varchar(20) NOT NULL default '',
    `transcript` varchar(20) NOT NULL default '',
    `trans_seq` int(11) default NULL,
    `trans_aa1` varchar(1) default NULL,
    `pdbid` varchar(20) NOT NULL default '',
    `chain` varchar(1) NOT NULL default '',
    `species` varchar(20) default NULL,
    `unp` varchar(20) default NULL,
    `chain_seq` int(11) default NULL,
    `chain_aa1` varchar(1) default NULL,
    `x` double default NULL,
    `y` double default NULL,
    `z` double default NULL,
    PRIMARY KEY  (`pdbid`,`chain`,`chain_seq`,`transcript`,`trans_seq`,`chr`,`start`,`end`,`var_chr`,`var_start`,`var_end`,`var_name`),
    KEY `var_chr` (`var_chr`,`var_start`,`var_end`),
    KEY `var_name` (`var_name`),
    KEY `bp` (`chr`,`start`,`end`),
    KEY `pdbid` (`pdbid`,`chain`,`chain_seq`)
    ) ENGINE=MyISAM DEFAULT CHARSET=latin1;
  """%dat_name)
# query.append("""LOAD DATA LOCAL INFILE '%s' INTO TABLE Intersect_Variants_%s
# 			(chr, start, end, strand, transcript, pdbid, chain, unp, seqres, 
# 			aa1, x, y, z, species, gene, var_chr, var_start, var_end, var_name)"""%(intersect_file,dat_name))
if ext == 'bed':
  var_file_cols = "var_chr,var_start,var_end,var_name"
elif ext == 'vcf':
  var_file_cols = "var_chr,var_start,var_name"
query.append("""LOAD DATA LOCAL INFILE '%s' INTO TABLE Intersect_Variants_%s 
                FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n'
                (chr,start,end,strand,gene,transcript,trans_seq,trans_aa1,
                 pdbid,chain,species,unp,chain_seq,chain_aa1,x,y,z,
                 %s)"""%(intersect_file,dat_name,var_file_cols))
if ext == 'vcf':
  query.append("""UPDATE Intersect_Variants_%s SET var_end=var_start+1"""%dat_name)

# Execute intersection queries
for q in query:
	c.execute(q)
con.close()	# Non-intersection pipeline will cause connection timeout

print "Intersection complete. You may now query the intersect table."
print "Annotating non-intersecting variants..."

# Find non-intersections
os.system("/usr/analysis/bin/intersectBed -v -b %s -a %s > %s"%(pdbmap_file,var_file,nointersect_file))

# Adjust the var_start,var_end positions to 1-indexing
os.system("""awk -F"\t" -v OFS="\t" '{print $1,$2+1,$3+1,$4}' %s > %sFIX"""%(nointersect_file,nointersect_file))
os.system("""mv -f %sFIX %s"""%(nointersect_file,nointersect_file))

# Annotate the non-intersecting ranges with Ensembl Transcripts and UniProt IDs
os.system("lib/get_protein.py %s > %s.temp"%(nointersect_file,nointersect_file))
os.system("mv -f %s.temp %s"%(nointersect_file,nointersect_file))

# Database Connection
try:
	con = MySQLdb.connect(host='gwar-dev',user='mike',passwd='cheezburger')
	c = con.cursor()
except MySQLdb.Error as e:
	print "There was an error connecting to the database.\n%s"%e
	sys.exit(1)

# Drop/Create/Load Non-Intersection Table
c.execute("use pdbmap_v7;")
if not append:
  query = ["DROP TABLE IF EXISTS Non_Intersect_Variants_%s;"%dat_name]
  query.append("""CREATE TABLE Non_Intersect_Variants_%s (
    `var_chr` varchar(10) NOT NULL default '',
    `var_start` int(11) NOT NULL default '0',
    `var_end` int(11) NOT NULL default '0',
    `var_name` varchar(50) default NULL,
    `gene` varchar(50) default NULL,
    `transcript` varchar(50) default NULL,
    `unp` varchar(20) default NULL,
    PRIMARY KEY  (`var_chr`,`var_start`,`var_end`,`transcript`),
    KEY `var_name` (`var_name`),
    KEY `uniprot` (`unp`)
    ) ENGINE=MyISAM DEFAULT CHARSET=latin1;
  """%dat_name)
# query.append("""LOAD DATA LOCAL INFILE '%s' INTO TABLE
# 			Non_Intersect_Variants_%s (var_chr, var_start, var_end, var_name,
# 			transcript, unp)"""%(nointersect_file,dat_name))
query.append("""LOAD DATA LOCAL INFILE '%s' INTO TABLE Non_Intersect_Variants_%s FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n'"""%(nointersect_file,dat_name))

# Execute non-intersection queries
for q in query:
	c.execute(q)
con.close()
