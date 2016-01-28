#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : PDBMapNetwork.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-06-27
# Description    : Service class for network analysis of structures
#                : and models in PDBmap.
#=============================================================================#

import subprocess as sp
import numpy as np
import sys,csv,os,re

class PDBMapNetwork():
  def __init__(self,reduce_loc,probe_loc):
    """ Initialization requires REDUCE and PROBE executable locations """
    self.reduce = reduce_loc
    self.probe  = probe_loc

  def generate_rin(self,structid,coord_file,backbone=False,pdb_bonds=False):
    """ Generates a residue interaction network using REDUCE and PROBE
        to identify non-covalent interactions. Inclusion of backbone 
        bonds and bonds explicitly defined in the coordinate file can
        be specified by the user """

    # Use REDUCE to create a temporary file with explicit hydrogens
    tempf1 = "temp/%d.TEMP"%multidigit_rand(10)
    tempf2 = "temp/%d.TEMP"%multidigit_rand(10)
    cmd    = ["zcat",coord_file,">",tempf1,"2>","/dev/null"]
    os.system(' '.join(cmd))
    cmd    = [self.reduce,"-FLIP",tempf1,">",tempf2,"2>","/dev/null"]
    os.system(' '.join(cmd))
    cmd    = ["rm","-f",tempf1]
    sp.Popen(cmd) # Remove temporary file

    # Use PROBE to identify non-covalent interactions
    cmd     = [self.probe,"-self","all",tempf2,"-unformated","-condense"]
    p       = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE)
    out,err = p.communicate()
    cmd     = ["rm","-f",tempf2]
    sp.check_call(cmd) # Remove temporary file

    # Process output into unique set of directed edges
    rows = [l.split() for l in out.split('\n') if l]
    rows = set([(r[1],int(r[2]),r[6],int(r[7])) for r in rows])
    # Remove self-referential edges
    rows = [list(r) for r in rows if not (r[0]==r[2] and r[1]==r[3])]
    
    # Add backbone connections and bonds defined in the coordinate file
    if backbone or pdb_bonds:
      if coord_file.split('.')[-1] == 'gz':
        import gzip
        fin = gzip.open(coord_file,'rb')
      else:
        fin = open(coord_file,'rb')
      for line in fin:
        if backbone and line.startswith('DBREF'):
          chain     = line[12]
          min_seqid = int(line[14:18])
          max_seqid = int(line[20:24])
          for i in range(min_seqid,max_seqid):
            # Link between residue and next residue in sequence
            rows.append([chain,i,chain,i+1])
            # Bidirectional
            rows.append([chain,i+1,chain,i])
        elif pdb_bonds and line.startswith('LINK'):
          link = [line[21],int(line[22:26])]
          link.extend([line[51],int(line[52:56])])
          rows.append(link)
        elif pdb_bonds and line.startswith('SSBOND'):
          ssbond = [line[15],int(line[17:21])]
          ssbond.extend([line[29],int(line[31:35])])
          rows.append(ssbond)
      fin.close()
    return rows


## Copied from biolearn
def multidigit_rand(digits):
  import random
  randlist = [random.randint(1,10) for i in xrange(digits)]
  multidigit_rand = int(''.join([str(x) for x in randlist]))
  return multidigit_rand