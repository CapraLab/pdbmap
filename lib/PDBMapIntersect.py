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
import sys,os,csv,random
import subprocess as sp

class PDBMapIntersect():
  def __init__(self,pdbmapio):
    """ Initialization requires a PDBMapIO object """
    self.io = pdbmapio

  def intersect(self,dlabel,slabel=None,dtype='Genomic'):
    # Intersects a structure set with a dataset using intersectBed
    # Slow, but intersection pointers are stored permanently (1-time op)
    # dtype options: Genomic, [Protein, Structural]
    
    # Define the two temp filenames
    temp1  = "temp/%d.TEMP"%multidigit_rand(10)
    temp2 = "temp/%d.TEMP"%multidigit_rand(10)
    try: # Ensure proper file cleanup
      # Query and write data ranges to temp file
      if dtype == 'Genomic':
        query  = "SELECT a.chr,a.start-1,a.end-1,a.gc_id,a.transcript FROM "
        query += "GenomicConsequence as a INNER JOIN GenomicData as b "
        query += "ON a.chr=b.chr AND a.start=b.start AND a.end=b.end "
        query += "WHERE b.label=%s"
      elif dtype == 'Protein':
        msg = "ERROR (PDBMapIntersect) Protein intersection not implemented."
        raise Exception(msg)
      elif dtype == 'Structural':
        msg = "ERROR (PDBMapIntersect) Structural intersection not implemented."
        raise Exception(msg)
      else:
        msg = "ERROR (PDBMapIntersect) %s intersection is not a valid option."
        raise Exception(msg%dtype)
      with open(temp1,'wb') as fout:
        writer = csv.writer(fout,delimiter='\t')
        data = self.io.secure_query(query,(dlabel,),cursorclass='SSCursor')
        print " # Fetching PDBMap::GenomicConsequence,GenomicData #"
        for row in data:
          writer.writerow(row)

      # Query and write PDBMap ranges to temp file, adjusted for UCSC indexing
      query  = "SELECT chr,start-1,end-1,pdbid,chain,chain_seqid,a.transcript FROM "
      query += "Transcript as a "
      query += "INNER JOIN Alignment as b "
      query += "ON a.transcript=b.transcript AND a.seqid=b.trans_seqid "
      if slabel:
        query += "WHERE b.pdbid=%s"
      with open(temp2,'wb') as fout:
        writer = csv.writer(fout,delimiter='\t')
        if slabel:
          structures = self.io.secure_query(query,(slabel,),cursorclass='SSCursor')
        else:
          structures = self.io.secure_query(query,cursorclass='SSCursor')
        print " # Fetching PDBMap::Alignment,Transcript #"
        for row in structures:
          writer.writerow(row)

      # IntersectBed is sensitive, ensure unix encoding
      os.system("dos2unix %s"%temp1)
      os.system("dos2unix %s"%temp2)

      ## Temp files written. Beginning Intersection and upload. ##

      # Intersect with intersectBed and upload output to PDBMap.IntersectGenome
      print " # Performing intersection #"
      cmd = ["intersectBed","-wb","-a",temp1,"-b",temp2]
      p = sp.Popen(cmd,stdout=sp.PIPE)
      parser = process_parser(p)
      nrows = self.io.upload_intersection(parse_intersection(parser))
     
      ## Intersection completed and uploaded. Cleaning up temp files ##    
    except: raise
    finally:
      # Remove temp files
      os.system("rm -f %s"%temp1)
      os.system("rm -f %s"%temp2)
      pass
    return(nrows) # Return the number of intersections
    
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
    # t_chr,t_start,t_end,pdbid,chain,seqid, \
    #   d_chr,d_start,d_end,gc_id = line.split('\t')
    d_chr,d_start,d_end,gc_id,d_trans,t_chr,t_start,t_end, \
    pdbid,chain,seqid,t_trans = line.split('\t')
    if d_trans != t_trans: continue # Not the same transcript
    seqid = int(seqid)
    gc_id = int(gc_id)
    # Return the direct reference
    yield (pdbid,chain,seqid,gc_id)

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
