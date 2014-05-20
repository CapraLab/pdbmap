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
import sys,os,csv,time,math
# from pymol import cmd
max_dist = -999 # Used by show_density and dist

class PDBMapVisualize():

  dens_mat = None # Class variable for show_desnity
  
  def __init__(self,io,pdb_dir='data/rcsb',modbase_dir='data/modbase'):
    """ Initialize the PDBMapVisualize object. """
    self.io = io
    self.pdb_dir = pdb_dir
    self.modbase_dir = modbase_dir

  def visualize_structure(self,pdbid,biounit=0,anno_list=['maf'],spectrum_range=None,group=None):
    """ Visualize the annotated dataset within a structure """
    pdbid = pdbid.lower()
    res  = self.io.load_structure(pdbid,biounit)
    # Output all annotations to results file
    cols = ['model','seqid','chain']
    cols.extend(anno_list)
    timestamp = str(time.strftime("%Y%m%d-%H"))
    params = {'structid':pdbid,'biounit':biounit,'annos':'-'.join(anno_list)}
    res_dir = 'results/pdbmap_%s_%s'
    if group:
      res_dir = res_dir%(group,timestamp)
    else:
      res_dir = res_dir%(params['structid'],timestamp)
    params['res_dir'] = res_dir
    if not os.path.exists(res_dir):
      os.system('mkdir -p %(res_dir)s'%params)
    with open("%(res_dir)s/%(structid)s_biounit%(biounit)d_vars_%(annos)s.txt"%params,'wb') as fout:
      writer = csv.writer(fout,delimiter='\t')
      writer.writerow(cols)
      out = [res[col] for col in cols]    # Extract columns as rows
      out  = [list(i) for i in zip(*out)] # Transpose back to columns
      for row in out:
        writer.writerow(row)

    # Correct submodel ID for undivided structures
    maxm = len(set(res['model']))
    if maxm < 2:
      res['model'] = [0 for m in res['model']]

    # Visualize individual annotations
    for anno in anno_list:
      if anno not in res:
        msg = "ERROR (PDBMapVisualize) Unknown feature %s\n"%anno
        sys.stderr.write(msg)
        continue # Move to next annotation
      cols = ['model','seqid','chain',anno]
      out  = [res[col] for col in cols if res[anno]]   # Extract columns as rows
      out  = [list(i) for i in zip(*out)] # Transpose back to columns
      minval,maxval = (999,-999) if not spectrum_range else spectrum_range
      tempf = "temp/%d.TEMP"%multidigit_rand(10)
      # Write the mappings to a temp file
      with open(tempf,'w') as fout:
        fout.write("#%s\n"%'\t'.join(cols))
        fout.write("attribute: %s\n"%anno)
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for row in out:
          # #0.model:resi.chain value [value ...]
          value = -1 if row[-1] is None else float(row[-1])
          fout.write("\t#0.%d:%d.%s\t%s\n"%(tuple(row)))
          if -1 < value < minval: minval=value
          if value > maxval: maxval=value
      minval,maxval = spectrum_range if spectrum_range else (minval,maxval)
      # Locate the asymmetric unit or biological assembly
      if biounit < 1:
        struct_loc = "%s/structures/all/pdb/pdb%s.ent.gz"%(self.pdb_dir,pdbid)
      else:
        struct_loc = "%s/biounit/coordinates/all/%s.pdb%d.gz"%(self.pdb_dir,pdbid,biounit)
      params = {'structid':pdbid,'biounit':biounit,'anno':anno,'tempf':tempf,
                'minval':minval,'maxval':maxval,'resis':out,'struct_loc':struct_loc}
      self.visualize(params,group=group)

  def visualize_model(self,modelid,biounit=0,anno_list=['maf'],spectrum_range=None,group=None):
    """ Visualize the annotated dataset within a model """
    res  = self.io.load_model(modelid)

    # Output all annotations to results file
    cols = ['model','seqid','chain']
    cols.extend(anno_list)
    timestamp = str(time.strftime("%Y%m%d-%H"))
    params = {'structid':modelid,'biounit':biounit}
    res_dir = 'results/pdbmap_%s_%s'
    if group:
      res_dir = res_dir%(group,timestamp)
    else:
      res_dir = res_dir%(params['structid'],timestamp)
    params['res_dir'] = res_dir
    if not os.path.exists(res_dir):
      os.system('mkdir -p %(res_dir)s'%params)
    with open("%(res_dir)s/%(structid)s_biounit%(biounit)d_vars_annotations.txt"%params,'wb') as fout:
      writer = csv.writer(fout,delimiter='\t')
      writer.writerow(cols)
      out = [res[col] for col in cols]    # Extract columns as rows
      out  = [list(i) for i in zip(*out)] # Transpose back to columns
      for row in out:
        writer.writerow(row)

    # Correct submodel ID for undivided structures
    maxm = len(set(res['model']))
    if maxm < 2:
      res['model'] = [0 for m in res['model']]

    # Visualize individual annotations
    for anno in anno_list:
      if anno not in res:
        msg = "ERROR (PDBMapVisualize) Unknown feature %s\n"%anno
        sys.stderr.write(msg)
        continue # Move to next annotation
      cols = ['model','seqid','chain',anno]
      out  = [res[col] for col in cols if res[anno]]   # Extract columns as rows
      out  = [list(i) for i in zip(*out)] # Transpose back to columns
      tempf = "temp/%d.TEMP"%multidigit_rand(10)
      # Write the mappings to a temp file
      with open(tempf,'w') as fout:
        fout.write("#%s\n"%'\t'.join(cols))
        fout.write("attribute: %s\n"%anno)
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for row in out:
          # #0.model:resi.chain value [value ...]
          fout.write("\t#0.%d:%d.%s\t%s\n"%(tuple(row)))
          if float(row[-1]) < minval: minval=float(row[-1])
          if float(row[-1]) > maxval: maxval=float(row[-1])
      minval,maxval = spectrum_range if spectrum_range else (minval,maxval)
      params = {'structid':modelid,'biounit':biounit,'anno':anno,'tempf':tempf,
                'minval':minval,'maxval':maxval,'resis':out,
                'struct_loc':"%s/models/model/%s.pdb.gz"%(self.modbase_dir,modelid)}
      self.visualize(params,group=group)

  def visualize_unp(self,unpid,anno_list=['maf'],spectrum_range=None):
    """ Visualize the annotated dataset associated with a protein """
    res_list  = self.io.load_unp(unpid)
    # Query all biological assemblies
    query = "SELECT DISTINCT biounit FROM Chain WHERE pdbid=%s"
    res   = self.io.secure_query(query,(entity),cursorclass=Cursor)
    biounits = [r[0] for r in res]
    # Visualize for each biological assembly
    for entity_type,entity in res_list:
      if entity_type == 'structure':
        for biounit in biounits:
          self.visualize_structure(entity,biounit,anno_list,spectrum_range,group=unpid)
      elif entity_type == 'model':
        for biounit in biounits:
          self.visualize_model(entity,biounit,anno_list,spectrum_range,group=unpid)
      else:
        msg = "ERROR (PDBMapVisualize) Invalid entity_type for %s: %s"%(entity,entity_type)
        raise Exception(msg)

  def visualize(self,params,group=None):
    # Visualize with PyMol
    timestamp = str(time.strftime("%Y%m%d-%H"))
    res_dir = 'results/pdbmap_%s_%s'
    if group:
      res_dir = res_dir%(group,timestamp)
    else:
      res_dir = res_dir%(params['structid'],timestamp)
    params['res_dir'] = res_dir
    if not os.path.exists(res_dir):
      os.system('mkdir -p %s'%res_dir)

    # Visualize with Chimera
    cmd  = "chimera --silent --script 'lib/PDBMapVisualize.py "
    keys = ['structid','biounit','anno','tempf','minval','maxval','struct_loc','res_dir']
    cmd += "%s'"%' '.join([str(params[key]) for key in keys])
    print cmd
    os.system(cmd)
    # replyobj.status("Chimera: Visualizing %(structid)s_biounit%(biounit)d_vars_%(anno)s"%params)
    # # Assign the annotation to relevant residues
    # rc("open %(struct_loc)s"%params)
    # rc("defattr %(tempf)s"%params)
    # # Color by chain
    # rc("rainbow chain")
    # # Identify all annotated residues as spheres
    # for resi in params['resis']:
    #   rc("shape sphere center #0.%d:%d.%s@CA radius 1 color grey"%tuple(resi[:3]))
    #   # If annotation is binary, color grey/red
    #   if (params['minval'],params['maxval']) == (0,1):
    #     rc("color red #0.%d:%d.%s@CA"%tuple(resi[:3]))
    # # If annotation is continuous, color spectrum
    # if not (params['minval'],params['maxval']) == (0,1):
    #   rc("rangecolor %(anno)s %(minval)s blue %(maxval)s red"%params)
    # # Export the image
    # rc("export POV-Ray %(res_dir)s/%(structid)s_biounit%(biounit)d_vars_%(anno)s.png"%params)

    # Visualize with PyMol
    # cmd  = "run lib/PDBMapVisualize.py; "
    # cmd += "load %(struct_loc)s,%(structid)s; "
    # cmd += "PDBMapVisualize.overwrite_bfactors('%(structid)s','%(tempf)s',var_spheres=True,spec_range=%(spec_range)s); "
    # cmd += "orient;mset 1 x360; movie.roll 1,180,1,axis=y; movie.roll 181,360,1,axis=x; "
    # cmd += "save %(res_dir)s/%(structid)s_vars_%(anno)s.pdb; "
    # cmd += "save %(res_dir)s/%(structid)s_vars_%(anno)s.pse; "
    # cmd += "png %(res_dir)s/%(structid)s_vars_%(anno)s.png, dpi=300, ray=2400,2400;"
    # cmd = cmd % params
    # print cmd
    # os.system('/usr/bin/pymol -c -d "%s"'%cmd)

    # os.system('rm -f %s'%params['tempf']) # clean up temp file

  # Executed in PyMol
  @classmethod
  def overwrite_bfactor(cls,pdbid,chain,resi,value):
    """ Overwrites the b-factors for all atoms in a residue """
    pdbid = pdbid.lower()
    if chain == 'A':
      # Some modbase models do not have chains, and have been dummy coded to A
      selection = "(%s and (chain %s or chain '') and resi %d)"%(pdbid,chain,resi)
    else:
      selection = "(%s and chain %s and resi %d)"%(pdbid,chain,resi)
    exp = "b=%f"%value
    command = "alter %s, %s"%(selection,exp)
    cmd.alter(selection,exp)
    return command

  # Executed in PyMol
  @classmethod
  def reset_bfactor(cls,pdbid,val=0.0):
    """ Sets b-factor for all atoms to the same value """
    exp = "b=%d"%val
    cmd.alter(pdbid,exp)

  # Called directly from PyMol
  @classmethod
  def overwrite_bfactors(cls,pdbid,scores,resis=None,discrete=False,
                         binary=False,displaytype='cartoon',surface=False,
                         var_spheres=False,spec_range=None):
    """ Visualizes variant characteristics or scores in a structure """
    col = 3 # specify the score column
    pdbid = pdbid.lower()
    if isinstance(scores,str): # if filename, read file
      fin = open(scores,'rU')
      reader = csv.reader(fin,delimiter='\t')
      if resis:
        scores = [row for row in reader if row[0].lower()==pdbid and int(row[2]) in resis]
      else:
        scores = [row for row in reader if row[0].lower()==pdbid]
      fin.close()
    for score in scores:
      if not score[col] or score[col] == 'NULL':
        score[col] = -1
    uniq_scores = [float(score[col]) for score in scores if float(score[col]) > 0]
    if spec_range and uniq_scores:
      min_score = spec_range[0]
      max_score = spec_range[1]
    elif uniq_scores:
      min_score = min(uniq_scores) - 0.0001
      max_score = max(uniq_scores)
    else:
      msg = "ERROR: %s contains no valid scores."%pdbid
      raise Exception(msg)

    PDBMapVisualize.reset_bfactor(pdbid,min_score-10)
    commands = [PDBMapVisualize.overwrite_bfactor(row[0],row[1],int(row[2]),float(row[col])) for row in scores]
    include = "br. b > %f"%min_score
    exclude = "br. b < %f"%min_score
    print("Minimum score: %s"%min_score)
    print("Maximum score: %s"%max_score)

    # Generate models
    if surface:
      cmd.create("surface_obj",pdbid)
      cmd.color("grey","surface_obj")
    cmd.hide("everything","all")
    cmd.show(representation=displaytype)
    cmd.set(name="cartoon_discrete_colors",value="on")
    if var_spheres:
        cmd.show(selection="name CA and br. b > -1.1",representation="spheres")
        cmd.set(name="sphere_scale",value=1)
    if surface:
      cmd.show(selection="surface_obj",representation="surface")
      cmd.set(selection="surface_obj",name="transparency",value=0.5)

    # Assign colors
    cmd.bg_color("white")
    cmd.color("grey","all")
    if binary:
      cmd.color("red", selection="br. b=1")
    elif discrete:
      discrete_colors = ["","firebrick","density","forest","tv_yellow","deeppurple","orange","orange","deepteal","white","black","grey"]
      for i in range(1,int(max_score)+1):
        cmd.color(discrete_colors[i],selection="br. b=%d"%i)
    else:
      cmd.spectrum("b","blue_red",selection=include,minimum=min_score,maximum=max_score)
    cmd.color("grey",exclude)
    cmd.orient("all")
    return commands

  # Executed in PyMol
  @classmethod
  def show_density(cls,pdbid,pdbmap_file,variants_file,radius):
    """ Calculates and visualizes variant density in a structure """
    radius = float(radius)
    if not PDBMapVisualize.dens_mat:
      # columns: chain,seqres,x,y,z
      fin = open(pdbmap_file,'rU')
      fin.readline() # burn the header
      reader = csv.reader(fin,delimiter='\t')
      residues = [row[1:] for row in reader if row[0].lower()==pdbid]
      fin = open(variants_file,'rU')
      fin.readline() # burn the header
      reader = csv.reader(fin,delimiter='\t')
      variants = [row[1:] for row in reader if row[0].lower()==pdbid]

      # Compute the distance "matrix"
      dist_mat = {}
      for res_chain,res_seqres,x,y,z in residues:
        res_pos = (x,y,z)
        for var_chain,var_seqres,i,j,k,chrom,start,end,var_name in variants:
          var_pos = (i,j,k)
          if (res_chain,res_seqres) not in dist_mat:
            dist_mat[(res_chain,res_seqres)] = []
          dist_mat[(res_chain,res_seqres)].append(dist((x,y,z),(i,j,k)))

      # Convert the distance matrix into a density "matrix"
      PDBMapVisualize.dens_mat = []
      for chain,seqres in dist_mat:
        distances = dist_mat[(chain,seqres)]
        for k in range(int(max_dist)):
          num_neighbors = len([x for x in distances if x <= k])
          PDBMapVisualize.dens_mat.append((chain,seqres,k,num_neighbors)) 
            
    # Reduce to requested k
    k_dense_map = [[pdbid,row[0],row[1],row[3]] for row in PDBMapVisualize.dens_mat if row[2]<=radius]
    PDBMapVisualize.overwrite_bfactors(pdbid,k_dense_map,displaytype='mesh')

  # # Once computed, store permanently, extracting different k values
  # PDBMapVisualize.dens_mat = None


def dist(ipos,jpos):
  global max_dist
  dist_x = abs(float(ipos[0]) - float(jpos[0])) ** 2
  dist_y = abs(float(ipos[1]) - float(jpos[1])) ** 2
  dist_z = abs(float(ipos[2]) - float(jpos[2])) ** 2
  xyz_dist = math.sqrt(dist_x+dist_y+dist_z)
  max_dist = max(max_dist,xyz_dist)
  return xyz_dist

## Copied from biolearn
def multidigit_rand(digits):
  import random
  randlist = [random.randint(1,10) for i in xrange(digits)]
  multidigit_rand = int(''.join([str(x) for x in randlist]))
  return multidigit_rand

# Add overwrite_bfactors to PyMol scope
# cmd.extend("PDBMapVisualize.overwrite_bfactors",PDBMapVisualize.overwrite_bfactors)
# cmd.extend("PDBMapVisualize.show_density",PDBMapVisualize.show_density)

# Chimera executable
if __name__ == '__main__':
  from chimera import runCommand as rc
  from chimera import replyobj
  # args = [sys.argv.split(' ')]
  args = sys.argv
  params = {'structid':args[1],'biounit':int(args[2]),'anno':args[3],'tempf':args[4],
            'minval':float(args[5]),'maxval':float(args[6]),'struct_loc':args[7],
            'res_dir':args[8],'resis':[]}
  # Read the residues back into memory
  with open(params['tempf']) as fin:
    for line in fin:
      if line.startswith('\t'):
        params['resis'].append(line.strip().split('\t'))

  # Visualize with Chimera
  replyobj.status("Chimera: Visualizing %(structid)s_biounit%(biounit)d_vars_%(anno)s"%params)
  # Assign the annotation to relevant residues
  rc("open %(struct_loc)s"%params)
  rc("defattr %(tempf)s"%params)
  # Initial colors
  rc("background solid white")
  rc("ribcolor dim grey")
  rc("ribbackbone")
  rc("~disp")
  rc("ribspline card smooth strand")
  # Identify all annotated residues as spheres
  for resi in params['resis']:
    rc("disp %s@ca"%resi[0])
    # # If annotation is binary, color grey/red
    # if (params['minval'],params['maxval']) == (0,1):
    #   rc("color red,a %s@ca"%resi[0])
  # Define variant representation
  rc("represent bs")
  rc("setattr m autochain 0")
  rc("setattr m ballScale .6")
  rc("rangecolor %(anno)s,a %(minval)s blue %(maxval)s red"%params)
  # Orient the image
  rc("define plane name p1 @ca:/%(anno)s>=%(minval)s"%params)
  rc("align p1")
  rc("~define")
  # Export the image
  rc("copy file %(res_dir)s/%(structid)s_biounit%(biounit)d_vars_%(anno)s.png width 800 height 800 units points dpi 72"%params)
  # Export the scene
  rc("save %(res_dir)s/%(structid)s_biounit%(biounit)d_vars_%(anno)s.py"%params)
