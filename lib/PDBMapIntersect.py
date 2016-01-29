#!/usr/bin/env python2.7
#
# Project        : PDBMap
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

  def quick_intersect(self,dlabel,slabel=None,dtype='Genomic'):
    # Intersects small datasets using a mysql join
    if dtype == "Genomic":
      pass
    elif dtype == 'Protein':
      msg = "ERROR (PDBMapIntersect) Protein intersection not implemented."
      raise Exception(msg)
    elif dtype == 'Structural':
      msg = "ERROR (PDBMapIntersect) Structural intersection not implemented."
      raise Exception(msg)
    else:
      msg = "ERROR (PDBMapIntersect) %s intersection is not a valid option."
      raise Exception(msg%dtype)
    query  = "INSERT IGNORE INTO GenomicIntersection "
    query += "SELECT a.label,c.label,c.structid,c.chain,c.chain_seqid,a.gc_id "
    query += "FROM GenomicConsequence as a "
    query += "INNER JOIN Transcript as b "
    query += "ON (a.transcript='' OR a.transcript=b.transcript) AND a.chr=b.chr "
    query += "AND a.start >= b.start AND a.end <= b.end "
    query += "INNER JOIN Alignment as c "
    query += "ON b.transcript=c.transcript AND b.seqid=c.trans_seqid "
    query += "WHERE a.label=%s "
    if slabel:
      query += "AND b.label=%s AND c.label=%s"
    if slabel:
      self.io.secure_command(query,(dlabel,slabel,slabel))
    else:
      self.io.secure_command(query,(dlabel,))
    return(-1) # No row feedback for quick-intersect

  def intersect(self,dlabel=None,slabel=None,dtype='Genomic',buffer_size=1):
    # Intersects a structure set with a dataset using intersectBed
    # Slow, but intersection pointers are stored permanently (1-time op)
    # dtype options: Genomic, [Protein, Structural]
    
    # Define the two temp filenames
    temp1  = "temp/%d.TEMP"%multidigit_rand(10)
    temp2  = "temp/%d.TEMP"%multidigit_rand(10)
    try: # Ensure proper file cleanup
      # Query and write data ranges to temp file
      if dtype == 'Genomic':
        query  = "SELECT a.chr,a.start-1,a.end-1,a.gc_id,a.transcript FROM "
        query += "GenomicConsequence as a INNER JOIN GenomicData as b "
        query += "ON a.label=b.label AND a.chr=b.chr AND a.start=b.start AND a.end=b.end "
        query += "AND a.transcript!='' "
        if dlabel:
          query += "WHERE a.label=%s "
        #query += "WHERE (a.consequence LIKE '%%missense_variant' OR a.consequence LIKE '%%synonymous_variant')"
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
        if dlabel:
          data = self.io.secure_query(query,(dlabel,),cursorclass='SSCursor')
        else:
          data = self.io.secure_query(query,cursorclass='SSCursor')
        print " # Fetching PDBMap::GenomicConsequence,GenomicData #"
        for row in data:
          writer.writerow(row)

      # Query and write PDBMap ranges to temp file, adjusted for UCSC indexing
      query  = "SELECT chr,start-1,end-1,structid,chain,chain_seqid,a.transcript FROM "
      query += "Transcript as a "
      query += "INNER JOIN Alignment as b "
      query += "ON a.label=b.label AND a.transcript=b.transcript AND a.seqid=b.trans_seqid "
      if slabel:
        query += "WHERE b.label=%s "
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
      sp.check_call(["dos2unix",temp1])
      sp.check_call(["dos2unix",temp2])

      ## Temp files written. Beginning Intersection and upload. ##

      # Intersect with intersectBed and upload output to PDBMap.IntersectGenome
      print " # Performing intersection #"
      cmd = ["intersectBed","-wb","-a",temp1,"-b",temp2]
      p = sp.Popen(cmd,stdout=sp.PIPE)
      parser = process_parser(p)
      nrows = self.io.upload_intersection(parse_intersection(parser),
                      dlabel=dlabel,slabel=slabel,buffer_size=buffer_size)

      # Remove temp files only if no exception
      #sp.check_call(["rm","-f",temp1])
      #sp.check_call(["rm","-f",temp2])
        
    except Exception as e: 
      msg  = "ERROR (PDBMapIntersect) Exception during "
      msg += "%s Intersection of %s and %s: %s"%(dtype,dlabel,slabel,str(e))
      sys.stderr.write(msg+'\n')
      raise
      # raise Exception(msg)
    finally:
      # # Remove temp files
      # sp.check_call(["rm","-f",temp1])
      # sp.check_call(["rm","-f",temp2])
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
    row = line.split('\t')
    if len(row) < 12 or not row[4]:
      # No transcript info from GenomicConsequence
      d_chr,d_start,d_end,gc_id,t_chr,t_start,t_end, \
      pdbid,chain,seqid,t_trans = row
    else:
      # Transcript checks only if transcript info available
      d_chr,d_start,d_end,gc_id,d_trans,t_chr,t_start,t_end, \
      pdbid,chain,seqid,t_trans = row
      if d_trans != t_trans: continue # Not the same transcript
    seqid = int(seqid)
    gc_id = int(gc_id)
    # Return the direct reference
    yield (pdbid,chain,seqid,gc_id)

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
