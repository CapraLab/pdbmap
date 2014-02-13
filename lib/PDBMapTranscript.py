#!/usr/bin/env python27
#
# Project        : PDBMap
# Filename       : PDBMapTranscript.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-12
# Description    : Defines the PDBMapTranscript class for representing a
#                : genome-mapped Ensembl transcript.
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv,commands

class PDBMapTranscript():
	
  def __init__(self,transid=None,unpid=None):
    if transid:
      self.from_trans(transid)
    elif unpid:
      self.from_unp(unpid)
    else:
    self.transcript = None
    self.gene       = None
    self.seqres     = None
    self.aa1        = None
    self.start      = None
    self.end        = None
    self.chr        = None
    self.strand     = None
    self.sequence   = None

  def from_unp(unpid):
    # Use UniParc to map UNP ID to Ensembl Transcript ID
    pass

  def from_trans(transid):
    # Use Ensembl Transcript ID to load transcript information
    cmd = "lib/transcript_to_genomic.pl %s"%transid
    status, output = commands.getstatusoutput(cmd)
    for line in output.split('\n'):
      t = line.split('\t')
      if line.startswith('#'): continue
      try:
        loc = (t[0], t[1], t[2])
      except IndexError:
        print "Index Error in run_intersectBed"
        print "%s\n"%line
        os.system("cp %s %s.debug"%(bed1,bed1))
        raise Exception("Malformed bed file for intersectBed.")

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
