#!/usr/bin/env python2.7
#
# Project        : Nucleotide Diversity Analysis
# Filename       : data_setup.py
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2015-10-14
# Description    : Generates data for the nucleotide diversity project. Pulls
#                : structural and genotypic information from PDBMap for a 
#                : given 1000 Genomes super-population. Data produced is ready
#                : for input into pi_analysis.py
#=============================================================================#
# Generates structure files for a 1KG super-population with the columns...
# Structure ID   : PDB structure ID
# Biounit        : Biological assembly #
# Model          : PDB model ID
# Chain          : PDB chain ID
# Seqid          : PDB residue #
# Insertion Code : PDB residue insertion code (or empty string)
# PFAM ACC       : PFAM accession ID
# PFAM Domain    : PFAM domain name
# x, y, z        : PDB residue coordinates (center of mass)
# Chromosome     : Variant chromosome
# Position       : Variant chromosomal position
# Name           : Variant name
# Consequence    : Missense (m), Synonymous SNP (s)
# Ancestral      : Ancestral allele
# Reference      : Reference allele
# Population MAF : Minor allele frequency in the sample population
# Genotypes      : Sample genotypes for this variant (no spaces)
#=============================================================================#
## Package Dependencies ##
import argparse,MySQLdb,csv,subprocess as sp,os,glob,shutil
#=============================================================================#
## Parse Command Line Options ##
desc  = "Generates data for the nucleotide diversity project for a given "
desc += "1000 Genomes super-population. Data generated is ready for input "
desc += "into pi_analyis.py. Requires 1000 Genomes VCF data."
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("population",type=str,
                    help="1000 Genomes super-population code")
parser.add_argument("vcfpath",type=str,
                    help="Path to 1000 Genomes VCF data")
args = parser.parse_args()
print "\n1000 Genomes super-population: %s"%args.population
print "1000 Genomes VCF path: %s"%args.vcfpath
#=============================================================================#
## Data Generation ##

# Create super-population sample ID file
PANELFILE = "integrated_call_samples_v3.20130502.ALL.panel"
pop = args.population
print "Creating %s.sampleid..."%pop
if pop.upper() == "AFR":
  print "The admixed ASW sub-population will be excluded from the AFR super-population."
  cmd = """grep -n %s %s/%s | grep -v ASW | cut -f1 | awk -F":" -vOFS="\t" '{print $2,$1-1}' > %s.sampleids"""%(pop.upper(),args.vcfpath,PANELFILE,pop)
else:
  cmd = """grep -n %s %s/%s | cut -f1 | awk -F":" -vOFS="\t" '{print $2,$1-1}' > %s.sampleids"""%(pop.upper(),args.vcfpath,PANELFILE,pop)
os.system(cmd)
if os.stat("%s.sampleids"%pop).st_size == 0:
  raise Exception("Invalid panel file or population name")

# Query the data from PDBMap
if not os.path.exists(".1kg3_%s.txt"%pop):
  print "Query structural data from PDBMap..."
  # Includes SNPs mapped through multiple transcripts,
  # which are filtered out later by pi_analysis.py
  MYSQLQUERY = "pi_query_template.sql"
  with open(MYSQLQUERY,'rb') as fin:
    sql = fin.read()%pop
  db = MySQLdb.connect(host="chgr2.accre.vanderbilt.edu",user="sivleyrm",
                       passwd="global-trifecta",db="pdbmap_v11")
  c = db.cursor()
  c.execute("SET group_concat_max_len=4096")
  c.close()
  c = db.cursor()
  with open(".1kg3_%s.txt"%pop,'wb') as fout:
    writer = csv.writer(fout,delimiter='\t')
    c.execute(sql)
    for row in c.fetchall():
      writer.writerow(row)
  db.close()

# Query super-population genotypes from 1000 Genomes
print "Query super-population genotypes from %s..."%args.vcfpath
GENOSCRIPT = "./query_genotypes.py"
with open(".1kg3_%s_genos.txt"%pop,'wb') as fout:
  cmd = [GENOSCRIPT,".1kg3_%s.txt"%pop,"%s.sampleids"%pop,args.vcfpath]
  p = sp.Popen(cmd,stdout=fout)
  p.wait()

# Split by structure
print "Separating by structure..."
if os.path.exists(pop):
  shutil.rmtree(pop)
os.makedirs(pop)
cmd  = "tail -n +2 .1kg3_%s_genos.txt | "%pop
cmd += """awk '{print>>"%s/"$1"_"$2"_%s_1kg3.txt"}'"""%(pop,pop)
os.system(cmd)

# Create super-population structfile
print "Creating %s.structfile..."%pop
with open("%s.structfile"%pop,'wb') as fout:
  for f in sorted(glob.glob("%s/*"%pop)):
    fout.write("%s\n"%os.path.abspath(f))

# Clean up the intermediate files
os.remove(".1kg3_%s.txt"%pop)
os.remove(".1kg3_%s_genos.txt"%pop)