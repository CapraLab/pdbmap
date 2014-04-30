#!/usr/bin/env python27
#
# Project        : PDBMap
# Filename       : PDBMapVisualize.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-27
# Description    : Visualization tool for PDBMap variation.
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv,time

class PDBMapVisualize():

  def __init__(self,io):
    """ Initialize the PDBMapVisualize object. """
    self.io = io

  def load_structure(self,pdbid,data_label='1kg'):
    """ Loads the structure from the PDBMap database """
    query = PDBMapVisualize.structure_query
    q = self.io.secure_query(query,qvars=(data_label,pdbid),cursorclass='DictCursor')
    res = {}
    for row in q:
      if not res:
        res = dict([(key,[val]) for key,val in row.iteritems()])
      else:
        for key,val in row.iteritems():
          res[key].append(val)
    return res

  def load_model(self,modelid,data_label='1kg'):
    """ Loads the structure from the PDBMap database """
    query = PDBMapVisualize.model_query
    q = self.io.secure_query(query,qvars=(data_label,modelid),cursorclass='DictCursor')
    res = {}
    for row in q:
      if not res:
        res = dict([(key,[val]) for key,val in row.iteritems()])
      else:
        for key,val in row.iteritems():
          res[key].append(val)
    return res

  def load_unp(self,unpid,data_label='1kg'):
    """ Identifies all associated structures, then pulls those structures. """
    query = PDBMapVisualize.unp_query
    q = self.io.secure_query(query,qvars=(data_label,unpid),cursorclass='Cursor')
    res = []
    for row in q:
      structid = row[0]
      if structid[0:4] == 'ENSP':
        res.append(load_model(structid,data_label))
      else:
        res.append(load_structure(structid,data_label))
    return res

  def visualize_structure(self,pdbid,data_label='1kg',annotation='maf',spectrum_range=None):
    """ Visualizes the data by annotation """
    res  = self.load_structure(pdbid,data_label)
    cols = ['pdbid','chain','seqid',annotation]
    out  = [res[col] for col in cols]   # Extract columns as rows
    out  = [list(i) for i in zip(*out)] # Transpose back to columns
    tempf = "temp/%d.TEMP"%multidigit_rand(10)
    # Write the mappings to a temp file
    with open(tempf,'w') as fout:
      writer = csv.writer(fout,delimiter='\t')
      for row in out:
        writer.writerow(row)
    # Visualize with PyMol
    res_dir = 'results/pdbmap_%s_%s/'%(pdbid,str(time.strftime("%Y%m%d-%H")))
    os.system('mkdir %s'%res_dir)
    cmd  = "run lib/overwrite_pymol.py; "
    cmd += "load /scratch/sivleyrm/pdb/structures/all/pdb/pdb%(pdbid)s.ent,%(pdbid)s; "
    cmd += "overwrite_bfactors('%(pdbid)s','%(tempf)s',col=4,var_spheres=True,spec_range=%(spec_range)s); "
    cmd += "orient;mset 1 x360; movie.roll 1,180,1,axis=y; movie.roll 181,360,1,axis=x; "
    cmd += "save %(res_dir)s/%(pdbid)s_vars_%(anno)s.pdb; "
    cmd += "save %(res_dir)s/%(pdbid)s_vars_%(anno)s.pse; "
    cmd += "png %(res_dir)s/%(pdbid)s_vars_%(anno)s.png, dpi=300, ray=2400,2400;"
    cmd = cmd % {'pdbid':pdbid,'anno':annotation,'tempf':tempf,'res_dir':res_dir,'spec_range':str(spectrum_range)}
    os.system('/usr/bin/pymol -c -d "%s"'%cmd)
    # os.system('rm -f %s'%tempf) # clean up temp file

  # Query definitions
  structure_query = """SELECT
    f.pdbid,g.chain,a.seqid,d.*,c.*
    FROM Residue as a
    INNER JOIN GenomicIntersection as b
    ON a.pdbid=b.pdbid AND a.chain=b.chain AND a.seqid=b.seqid
    INNER JOIN GenomicConsequence as c
    ON b.gc_id=c.gc_id
    INNER JOIN GenomicData as d
    ON c.label=d.label AND c.chr=d.chr AND c.start=d.start AND c.end=d.end AND c.name=d.name
    INNER JOIN Structure as f
    ON a.pdbid=f.pdbid
    INNER JOIN Chain as g
    ON a.pdbid=g.pdbid AND a.chain=g.chain
    WHERE a.label='uniprot-pdb' AND c.label=%s AND a.pdbid=%s;"""
  model_query = """SELECT
    f.modelid,g.chain,a.seqid,d.*,c.*
    FROM Residue as a
    INNER JOIN GenomicIntersection as b
    ON a.pdbid=b.pdbid AND a.chain=b.chain AND a.seqid=b.seqid
    INNER JOIN GenomicConsequence as c
    ON b.gc_id=c.gc_id
    INNER JOIN GenomicData as d
    ON c.label=d.label AND c.chr=d.chr AND c.start=d.start AND c.end=d.end AND c.name=d.name
    INNER JOIN Model as f
    ON a.pdbid=f.modelid
    INNER JOIN Chain as g
    ON a.pdbid=g.pdbid AND a.chain=g.chain
    WHERE a.label='uniprot-pdb' AND c.label=%s AND a.modelid=%s;"""
  unp_query = """SELECT b.pdbid FROM 
    GenomicIntersection as a
    INNER JOIN Residue as b
    ON a.pdbid=b.pdbid AND a.chain=b.chain
    INNER JOIN GenomicConsequence as c
    ON a.gc_id=c.gc_id
    WHERE c.label=%s AND b.unp=%s;"""

## Copied from biolearn
def multidigit_rand(digits):
  import random
  randlist = [random.randint(1,10) for i in xrange(digits)]
  multidigit_rand = int(''.join([str(x) for x in randlist]))
  return multidigit_rand