#!/usr/bin/env python27
#
# Project        : PDBmap
# Filename       : PDBMapIntersect.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-27
# Description    : Service class designed to intersect structural and 
#                : genomic, and sequence datasets. Does not alter original 
#                : datasets. Calculates intersection with intersectBed and 
#                : stores the crossreference in the PDBMap.IntersectGenome 
#                : database.
#=============================================================================#

# See main check for cmd line parsing
from lib import PDBMapIO
import csv

class PDBMapIntersect():
  def __init__(self,pdbmapio):
    """ Initialization requires a PDBMapIO object """
    self.io = pdbmapio

  @classmethod
  def intersect(slabel,dlabel,dtype='Genomic'):
    # Intersects a structure set with a dataset using intersectBed
    # dtype options: Genomic, Protein, Structural

    # Query and write PDBMap ranges to temp file
    query  = "SELECT chr,start,end,pdbid,chain,seqid FROM "
    query += "Transcript as a "
    query += "INNER JOIN Alignment as b "
    query += "ON a.transcript=b.transcript AND a.seqid=b.trans_seqid "
    query += "WHERE b.pdbid=%s"
    stemp  = "working/%d.TEMP"%multidigit_rand(10)
    with open(temp1,'wb') as fout:
      writer = csv.writer(fout)
      pdbmap = self.io.secure_query(query,(slabel,))
      for row in pdbmap:
        writer.writerow(row)
    
    # Query and write data ranges to temp file
    if dtype == 'Genomic':
      query  = "SELECT chr,start,end,name FROM "
      query += "GenomicData WHERE label=%s"
    elif dtype == 'Protein':
      msg = "ERROR: (PDBMapIntersect) Protein intersection not implemented."
      raise Exception(msg)
    elif dtype == 'Structural':
      msg = "ERROR: (PDBMapIntersect) Structural intersection not implemented."
      raise Exception(msg)
    else:
      msg = "ERROR: (PDBMapIntersect) %s intersection is not a valid option."
      raise Exception(msg%dtype)
    temp2 = "working/%d.TEMP"%multidigit_rand(10)
    with open(temp2,'wb') as fout:
      writer = csv.writer(fout)
      data = self.io.secure_query(query,(dlabel,))
      for row in data:
        writer.writerow(row)

    ## Temp files written. Beginning Intersection and upload. ##

    # Intersect with intersectBed and upload output to PDBMap.IntersectGenome
    cmd = ["intersectBed","-wao","-a",temp1,"-b",temp2]
    p = sp.Popen(cmd,stdout=sp.PIPE)
    parser = process_parser(p)
    self.io.upload_intersection(parse_intersection(parser))
   
    ## Intersection completed and uploaded. Cleaning up temp files ##    

    # Remove temp files
    os.system("rm -f %s"%temp1)
    os.system("rm -f %s"%temp2)
    

## Copied from biolearn
def multidigit_rand(digits):
  randlist = [random.randint(1,10) for i in xrange(digits)]
  multidigit_rand = int(''.join([str(x) for x in randlist]))
  return multidigit_rand

## Copied from snp_stats
def process_parser(p):
  while True:
    line = p.stdout.readline()
    if not line: break
    yield line

## Adapted from snp_stats
def parse_intersection(parser):
  for line in parser:
    line = line.strip()
    if not line or line[0] == "#": continue
    row = line.split('\t')
    yield row[0:10] # original columns, intersected

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
