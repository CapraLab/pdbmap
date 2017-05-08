#!/usr/bin/env python2.7
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
#                : genome-mapped Ensembl transcript. Depends on PDBMapProtein.
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv,commands
from PDBMapProtein import PDBMapProtein

class PDBMapTranscript():
	
  # Transcript query cache, keyed on transcript IDs
  trans_cache = {}
  CACHE_ACCESS_MIN = 25

  def __init__(self,transcript,protein,gene,sequence):
    # Define transcript, gene, and sequence
    self.transcript = transcript
    self.protein    = protein
    self.gene       = gene
    # Sequence | key:   seqid
    # Sequence | value: (rescode,chr,start,end,strand)
    self.sequence   = sequence

  @classmethod
  def cache_transcript(cls,transid,transcript):
    """ Caches a transcript result from a transid query """
    # Cache the transid->transcript query
    PDBMapTranscript.trans_cache[transid] = (transcript,0)
    # Update access counts
    for tid,(tobj,count) in PDBMapTranscript.trans_cache.iteritems():
      PDBMapTranscript.trans_cache[tid] = (tobj,count+1)
    # Remove infrequently accessed entries
    PDBMapTranscript.trans_cache = dict([(x,(y,z)) for x,(y,z) in \
                PDBMapTranscript.trans_cache.iteritems() \
                if z < PDBMapTranscript.CACHE_ACCESS_MIN])

  @classmethod
  def query_from_unp(cls,unpid):
    """ Use UniProt to map UniProt ID to Ensembl Transcript ID """
    if PDBMapProtein.check_loaded() and unpid in PDBMapProtein._sec2prim:
      msg  = "  WARNING (UniProt) %s is a secondary UniProt AC. "%unpid
      unpid = PDBMapProtein._sec2prim[unpid]
      msg += "Using primary AC: %s\n"%unpid
      sys.stderr.write(msg)
    transids = PDBMapProtein.unp2enst(unpid)
    if len(transids) < 1:
      msg = "\n  WARNING (query_from_unp) No transcripts match %s\n"%unpid
      sys.stderr.write(msg)
    # Query all transcript candidates and return
    res = []
    for transid in transids:
      trans = PDBMapTranscript.query_from_trans(transid)
      if trans:
        res.append(trans)
    # Only report an error if NO transcript matches were identified.
    if not res:
      msg = "   WARNING (PDBMapTranscript) No valid transcripts identified for %s"%unpid
      sys.stderr.write(msg)
    return res

  @classmethod
  def query_from_trans(cls,transid):
    """ Use Ensembl Transcript ID to load transcript information """
    # Check for cached transcript query result
    if transid in PDBMapTranscript.trans_cache:
      trans = PDBMapTranscript.trans_cache[transid][0] # exclude the count
      return trans
    # Query the Ensembl API for the transcript
    cmd = "perl lib/transcript_to_genomic.pl %s"%transid
    status, output = commands.getstatusoutput(cmd)
    if status > 0:
      # msg = "   WARNING (transcript_to_genomic.pl) %s: %s\n"%(transid,output)
      # sys.stderr.write(msg)
      PDBMapTranscript.cache_transcript(transid,None)
      return None
    sequence = {} # store sequence keyed on seqid
    for line in output.split('\n'):
      if line.startswith('#'): continue
      fields = line.split('\t')
      transcript = fields[0]
      protein    = PDBMapProtein.ensp2unp(fields[1])[0]
      gene       = fields[2]
      seqid      = int(fields[3])
      rescode    = fields[4].upper()
      if rescode not in aa_code_map.values():
        # msg = "  Non-standard amino acid replacement: %s[%d]->%s with %s\n"%(transid,seqid,rescode,'S')
        # sys.stderr.write(msg)
        rescode  = 'S' # replace non-standard amino acids with Serine
      start      = int(fields[5])
      end        = int(fields[6])
      chrom      = fields[7]
      # If transcript on a haplotype chromosome, ignore
      if chrom not in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8',
                      'chr9','chr10','chr11','chr12','chr13','chr14','chr15',
                      'chr16','chr17','chr18','chr19','chr20','chr21','chr22',
                      'chrX','chrY','chrMT']:
        msg = "  WARNING (transcript_to_genomic.pl) Ignoring non-standard chromosome %s\n"%chrom
        sys.stderr.write(msg)
        return None
      strand     = int(fields[8])
      sequence[seqid] = (rescode,chrom,start,end,strand)
    # Return a new PDBMapTranscript object
    trans = PDBMapTranscript(transcript,protein,gene,sequence)
    PDBMapTranscript.cache_transcript(transid,trans)
    return trans

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

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.\n")
  sys.exit(1)
