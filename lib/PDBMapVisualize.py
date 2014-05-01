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

  def __init__(self,io,pdb_dir='data/rcsb',modbase_dir='data/modbase'):
    """ Initialize the PDBMapVisualize object. """
    self.io = io
    self.pdb_dir = pdb_dir
    self.modbase_dir = modbase_dir

  def visualize_structure(self,pdbid,annotation='maf',spectrum_range=None):
    """ Visualize the annotated dataset within a structure """
    res  = self.io.load_structure(pdbid)
    cols = ['pdbid','chain','seqid',annotation]
    out  = [res[col] for col in cols]   # Extract columns as rows
    out  = [list(i) for i in zip(*out)] # Transpose back to columns
    tempf = "temp/%d.TEMP"%multidigit_rand(10)
    # Write the mappings to a temp file
    with open(tempf,'w') as fout:
      writer = csv.writer(fout,delimiter='\t')
      for row in out:
        writer.writerow(row)
    params = {'structid':pdbid,'anno':annotation,'tempf':tempf,
              'spec_range':str(spectrum_range),
              'struct_loc':"%s/pdb%s.ent.gz"%(self.pdb_dir,pdbid)}
    self.visualize(params)

  def visualize_model(self,modelid,annotation='maf',spectrum_range=None):
    """ Visualize the annotated dataset within a model """
    res  = self.io.load_model(modelid)
    cols = ['modelid','chain','seqid',annotation]
    out  = [res[col] for col in cols]   # Extract columns as rows
    out  = [list(i) for i in zip(*out)] # Transpose back to columns
    tempf = "temp/%d.TEMP"%multidigit_rand(10)
    # Write the mappings to a temp file
    with open(tempf,'w') as fout:
      writer = csv.writer(fout,delimiter='\t')
      for row in out:
        writer.writerow(row)
    params = {'structid':modelid,'anno':annotation,'tempf':tempf,
              'spec_range':str(spectrum_range),
              'struct_loc':"%s/models/model/%s.pdb.gz"%(self.modbase_dir,modelid)}
    self.visualize(params)

    def visualize_unp(self,unpid,annotation='maf',spectrum_range=None):
      """ Visualize the annotated dataset associated with a protein """
      res_list  = self.io.load_unp(unpid)
      for entity_type,entity in res_list:
        if entity_type == 'structure':
          self.visualize_structure(entity)
        elif entity_type == 'model':
          self.visualize_model(entity)
        else:
          msg = "ERROR (PDBMapVisualize) Invalid entity_type for %s: %s"%(entity,entity_type)
          raise Exception(msg)

  def visualize(self,params):
    # Visualize with PyMol
    res_dir = 'results/pdbmap_%s_%s'%(params['structid'],str(time.strftime("%Y%m%d-%H")))
    params['res_dir'] = res_dir
    if not os.path.exists(res_dir):
      os.system('mkdir %s'%res_dir)
    cmd  = "run lib/overwrite_pymol.py; "
    cmd += "load %(struct_loc)s,%(structid)s; "
    cmd += "overwrite_bfactors('%(structid)s','%(tempf)s',col=4,var_spheres=True,spec_range=%(spec_range)s); "
    cmd += "orient;mset 1 x360; movie.roll 1,180,1,axis=y; movie.roll 181,360,1,axis=x; "
    cmd += "save %(res_dir)s/%(structid)s_vars_%(anno)s.pdb; "
    cmd += "save %(res_dir)s/%(structid)s_vars_%(anno)s.pse; "
    cmd += "png %(res_dir)s/%(structid)s_vars_%(anno)s.png, dpi=300, ray=2400,2400;"
    cmd = cmd % params
    print cmd
    os.system('/usr/bin/pymol -c -d "%s"'%cmd)
    # os.system('rm -f %s'%params['tempf']) # clean up temp file

## Copied from biolearn
def multidigit_rand(digits):
  import random
  randlist = [random.randint(1,10) for i in xrange(digits)]
  multidigit_rand = int(''.join([str(x) for x in randlist]))
  return multidigit_rand