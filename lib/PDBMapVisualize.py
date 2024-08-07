#!/usr/bin/env python
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
import sys,os,csv,time,math,platform

# from pymol import cmd
max_dist = -999 # Used by show_density and dist

class PDBMapVisualize():

  dens_mat = None # Class variable for show_desnity
  
  def __init__(self,io,pdb_dir='data/rcsb'):
    """ Initialize the PDBMapVisualize object. """
    self.io = io
    self.pdb_dir = pdb_dir

  def visualize_structure(self,structid,biounit=0,anno_list=['maf'],eps=None,mins=None,spectrum_range=[],group=None,colors=[],permute=False,syn=False):
    """ Visualize the annotated dataset within a structure """
    biounit = 0 if biounit<0 else biounit
    print("Visualizing %s.%s..."%(structid,biounit))
    structid = structid.lower()
    res  = self.io.load_structure(structid,biounit,raw=True,syn=syn)
    if not res:
      msg = "WARNING (PDBMapVisualize) No variants for %s, biounit %s\n"%(structid,biounit)
      sys.stderr.write(msg)
      return
    modelflag = True if res['mpqs'][0] else False

    if anno_list == ['.']:
      anno = self.io.dlabel
      # If anno starts with number, convert to number-name
      try: 
        int(anno[0])
      except:
        pass
      else:
        import inflect
        convA = inflect.engine().number_to_words(anno[0]) + anno[1:]
        anno = convA
        print("\nConverted annotation label from %s to %s"%(anno,convA))
      res[anno] = [1 if l==self.io.dlabel else None for l in res["dlabel"]]
      anno_list = [anno]

    # If synonymous variants are requested, switch the issnp flag
    if syn:
      print("\nSwitching from nonsynonymous to synonymous variants...")
      res['issnp'] = [int('synonymous_variant' in c) if c else 0 for c in res['consequence']]

    # Convert MAF values into DAF values
    for anno in anno_list:
      if anno in ['daf','amr_daf','eas_daf','sas_daf','eur_daf','afr_daf']:
        maf = 'maf' if anno == 'daf' else anno.replace('daf','af')
        res[anno] = [res[maf][i] if not res['anc_allele'][i] or res['anc_allele'][i]=='.' or res['anc_allele'][i]==res['ref_allele'][i] else 1.-res[maf][i] for i in range(len(res[maf]))]

    # Ensure any user-supplied annotations are properly formatted
    anno_list = [a.replace('.','_') for a in anno_list]
    anno_list = [a.replace(' ','_') for a in anno_list]
    if any(a not in res for a in anno_list):
      res = [] # reset the result
      for anno in anno_list:
        # If any specified annotation isn't in the default return
        if anno not in res:
          # Join with the user-supplied annotations
          nres = self.io.load_structure(structid,biounit,useranno=True,raw=True)
          if syn:
            # If synonymous variants are requested, switch the issnp flag
            print("\nSwitching from nonsynonymous to synonymous variants...")
            nres['issnp'] = [int(c.contains('synonymous_variant')) for c in nres['consequence']]
          if not res:
            # Initialize res with the first user-annotated result
            res = nres
          break # all missing annotations filled when the first is found missing

    # If no residues are marked as SNVs, but some have the requested annotation,
    # print a warning and skip the consequence-checks
    if any([a for anno in anno_list for a in res[anno]]) and not any(res["issnp"]):
      msg = "\nWARNING: No residues flagged as SNVs. Visualizing all annotated variants.\n\n"
      sys.stderr.write(msg)
      res['issnp'] = [int(a!=None) for a in res[anno]]

    # Reduce to variable residues
    for key in res:
      if key!="issnp" and any(res['issnp']):
        res[key] = [r for i,r in enumerate(res[key]) if res['issnp'][i]]
      elif key not in ("issnp",anno):
        # Assume we're working with structural data
        res[key] = [r for i,r in enumerate(res[key]) if res[anno][i]!=None]
    if not any(res['issnp']):
      res[anno] = [r for i,r in enumerate(res[anno]) if res[anno][i]!=None]
    del res['issnp']

    # Report the final annotation count
    print("(%d annotated residues)"%len(res['seqid']))
    print("(%d non-duplicate residues)"%len(set(zip(res['chain'],res['seqid']))))
    print("(%d with non-zero/null values)"%len([x for x in res[anno] if x]))

    # Determine the first residue for renumbering
    resrenumber = {}
    for chain in set(res["chain"]):
      resrenumber[chain] = min(r for i,r in enumerate(res["aln_trans_seqid"]) if res["chain"][i]==chain)
    # Correct submodel ID for undivided structures
    maxm = len(set(res['model']))
    if maxm < 2:
      res['model'] = [0 for m in res['model']]

    # If DBSCAN Eps specified
    if eps and mins:
      import numpy as np
      from sklearn.cluster import DBSCAN
      from scipy.spatial.distance import pdist,squareform
      cols = ['model','seqid','chain']
      location = np.array([list(coord) for coord in zip(res['x'],res['y'],res['z'])])
      seen = set()
      seen_add = seen.add
      newres = dict((k,[]) for k in list(res.keys()))
      for i in range(len(res['x'])):
        coord = (res['x'][i],res['y'][i],res['z'][i])
        if not (coord in seen or seen_add(coord)):
          for key in res:
            newres[key].append(res[key][i])
      res = newres
      location = np.array([list(coord) for coord in zip(res['x'],res['y'],res['z'])])
      distance = squareform(pdist(location,'euclidean'))
      clusters = DBSCAN(eps=eps,min_samples=mins,metric='precomputed').fit(distance).labels_
      noise    = int(min(clusters) >= 0)
      clusters = [c+1 for c in clusters]
      res['clusters_%d_%d'%(eps,mins)] = clusters
      cls =  ['grey','red','blue','green','orange','yellow']
      cls += ['brown','purple','black','white','pink','magenta']
      cls += ['sienna','firebrick','khaki','turquoise']
      cls += ['cyan','slate','chartreuse']
      cls *= len(set(clusters)) / len(cls) + 1
      colors.append(cls[noise:len(set(clusters))+noise])
      anno_list.append('clusters_%d_%d'%(eps,mins))

    # Output all annotations to results file
    cols = ['model','seqid','chain']
    cols.extend(anno_list)
    timestamp = str(time.strftime("%Y%m%d-%H"))
    params = {'structid':structid,'biounit':biounit,
              'annos':'-'.join(anno_list),
              'vtype':['nonsynonymous','synonymous'][int(syn)]}
    res_dir = 'results/pdbmap_%s_%s_%s'
    if group:
      res_dir = res_dir%(self.io.dlabel,group,timestamp)
    else:
      res_dir = res_dir%(self.io.dlabel,params['structid'],timestamp)
    params['res_dir'] = res_dir
    if not os.path.exists(res_dir):
      os.system('mkdir -p %(res_dir)s'%params)
    with open("%(res_dir)s/%(structid)s_biounit%(biounit)s_%(annos)s_%(vtype)s.txt"%params,'wb') as fout:
      writer = csv.writer(fout,delimiter='\t')
      writer.writerow(cols)
      out = [res[col] for col in cols]    # Extract columns as rows
      out  = [list(i) for i in zip(*out)] # Transpose back to columns
      for row in out:
        writer.writerow(row)

    # Visualize individual annotations
    for a,anno in enumerate(anno_list):
      if anno not in res:
        msg = "ERROR (PDBMapVisualize) Unknown feature %s\n"%anno
        sys.stderr.write(msg)
        continue # Move to next annotation
      cols = ['model','seqid','chain',anno]
      out  = [[r for i,r in enumerate(res[col]) if res[anno][i]] for col in cols] # Extract columns as rows
      out  = [list(i) for i in zip(*out)]            # Transpose back to columns
      try:
        # For duplicate residues (multi-transcripts, multi-snp codons, etc)
        # use the highest value of the annotation (and hope its not pvalues)
        redout,prevrow = [],[]
        for row in out:
          if ('maf' in anno or 'daf' in anno) and float(row[-1])>0.5:
            row[-1] = 1-float(row[-1]) # convert ref/alt to major/minor
          if row[:-1] == prevrow[:-1]:
            if float(row[-1]) > float(prevrow[-1]):
              redout[-1] = row
          else:
            redout.append(row)
          prevrow = row
        out = redout
      except:
        pass
      minval,maxval = (999,-999)
      attrf = "%(res_dir)s/%(structid)s_biounit%(biounit)s_%%s_%(vtype)s.attr"%params%anno
      # Write the mappings to a temp file
      with open(attrf,'w') as fout:
        fout.write("#%s\n"%'\t'.join(cols))
        fout.write("attribute: %s\n"%anno)
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for row in out:
          # #0.model:resi.chain value [value ...]
          value = -1 if row[-1] is None else row[-1]
          if not modelflag:
            fout.write("\t#0.%d:%d.%s\t%s\n"%(tuple(row)))
          else: # do not include chain specifier for models
            fout.write("\t#0.%d:%d\t%0.6f\n"%(tuple(row[:2]+[row[-1]])))
          try:
            float(value)
            if -1 < value < minval: minval=value
            if value > maxval: maxval=value
          except:
            minval,maxval = -1,-1
      # If the variant type is synonymous, create a second attribute file
      # with "synonymous" as the attribute name. Useful when comparing
      # one attribute between synonymous and nonsynonymous variants
      if params['vtype'] == 'synonymous':
        attrf = "%(res_dir)s/%(structid)s_biounit%(biounit)s_%(vtype)s.attr"%params
        with open(attrf,'w') as fout:
          fout.write("#%s\n"%'\t'.join(cols))
          fout.write("attribute: synonymous\n")
          fout.write("match mode: 1-to-1\n")
          fout.write("recipient: residues\n")
          for row in out:
            # #0.model:resi.chain value [value ...]
            value = -1 if row[-1] is None else row[-1]
            if not modelflag:
              fout.write("\t#0.%d:%d.%s\t%s\n"%(tuple(row)))
            else: # do not include chain specifier for models
              fout.write("\t#0.%d:%d\t%0.6f\n"%(tuple(row[:2]+[row[-1]])))
      minval,maxval = (minval,maxval) if not spectrum_range else spectrum_range[a]
      colors = None if not colors else colors[a]
      # Locate the asymmetric unit or biological assembly
      if modelflag:
        # Extract the Ensembl protein ID from the ModBase model ID
        from lib import PDBMapModel
        struct_loc = PDBMapModel.get_coord_file(structid.upper())
      elif int(biounit) == 0:
        struct_loc = "%s/structures/pdb%s.ent.gz"%(self.pdb_dir,structid)
      else:
        struct_loc = "%s/biounits/%s.pdb%s.gz"%(self.pdb_dir,structid,biounit)
      params = {'structid':structid,'biounit':biounit,'anno':anno,'attrf':attrf,'colors':colors,'vtype':params['vtype'],
                'minval':minval,'maxval':maxval,'resis':out,'struct_loc':struct_loc,'resrenumber':resrenumber}
      self.visualize(params,group=group)

  def visualize_unp(self,unpid,anno_list=['maf'],eps=None,mins=None,spectrum_range=[],colors=[],syn=False):
    """ Visualize the annotated dataset associated with a protein """
    print("Visualizing protein %s"%(unpid))
    res_list  = self.io.load_unp(unpid)
    # Visualize for each biological assembly
    for entity_type,entity in res_list:
      if entity_type == 'structure':
        if self.io.is_nmr(entity):
          biounits = [-1]
        else:
          # Query all biological assemblies
          query = "SELECT DISTINCT biounit FROM Chain WHERE label=%s AND structid=%s AND biounit>0"
          res   = self.io.secure_query(query,(self.io.slabel,entity),cursorclass='Cursor')
          biounits = [r[0] for r in res]
        for biounit in biounits:
          print("Visualizing structure %s.%s"%(entity,biounit))
          self.visualize_structure(entity,biounit,anno_list,eps,mins,spectrum_range,group=unpid,colors=colors,syn=syn)
      elif entity_type == 'model':
        # Query all biological assemblies
        biounits = [-1]
        for biounit in biounits:
          print("Visualizing model %s.%s"%(entity,biounit))
          self.visualize_structure(entity,biounit,anno_list,eps,mins,spectrum_range,group=unpid,colors=colors,syn=syn)
      else:
        msg = "ERROR (PDBMapVisualize) Invalid entity_type for %s: %s"%(entity,entity_type)
        raise Exception(msg)

  def visualize_all(self,anno_list=['maf'],eps=None,mins=None,spectrum_range=[],colors=[],syn=False):
    """ Visualize all structures and models for the annotated dataset """
    print("Visualizing dataset %s"%self.io.dlabel)
    query = "SELECT DISTINCT structid FROM GenomicIntersection WHERE dlabel=%s"
    res   = [r for r in self.io.secure_query(query,(self.io.dlabel,),cursorclass='Cursor')]
    structures = [r[0] for r in res if self.io.detect_entity_type(r[0]) == 'structure']
    for s in structures:
      if self.io.is_nmr(s):
        biounits = [-1]
      else:
        query = "SELECT DISTINCT biounit FROM Chain WHERE structid=%s AND biounit>0"
        bres   = self.io.secure_query(query,(s,),cursorclass='Cursor')
        biounits = [r[0] for r in bres]
      for b in biounits:
        self.visualize_structure(s,b,anno_list,eps,mins,spectrum_range,group='all',colors=colors,syn=syn)
    models = [r[0] for r in res if self.io.detect_entity_type(r[0]) == 'model']
    for m in models:
      biounits = [-1]
      bres   = self.io.secure_query(query,(m,),cursorclass='Cursor')
      biounits = [r[0] for r in bres]
      for b in biounits:
        self.visualize_model(m,b,anno_list,eps,mins,spectrum_range,group='all',colors=colors,syn=syn)

  def visualize(self,params,group=None):
    # Visualize with PyMol
    timestamp = str(time.strftime("%Y%m%d-%H"))
    res_dir = 'results/pdbmap_%s_%s_%s'
    if group:
      res_dir = res_dir%(self.io.dlabel,group,timestamp)
    else:
      res_dir = res_dir%(self.io.dlabel,params['structid'],timestamp)
    params['res_dir'] = res_dir
    if not os.path.exists(res_dir):
      os.system('mkdir -p %s'%res_dir)

    # Visualize with Chimera
    params['minval'] = "%0.6f"%params['minval']
    params['maxval'] = "%0.6f"%params['maxval']
    params['colors'] = '-' if not params['colors'] else ','.join(params['colors'])
    keys   = ['structid','biounit','anno','attrf','minval','maxval','struct_loc','res_dir','colors','resrenumber','vtype']
    script = '"lib/PDBMapVisualize.py %s"'%' '.join([str(params[key]).translate(None,' ') for key in keys])
    cmd    = "TEMP=$PYTHONPATH; unset PYTHONPATH; chimera --nogui --silent --script %s; export PYTHONPATH=$TEMP"%script
    # Allow Mac OSX to use the GUI window
    if platform.system() == 'Darwin':
      cmd  = "TEMP=$PYTHONPATH; unset PYTHONPATH; chimera --silent --script %s; export PYTHONPATH=$TEMP"%script
    try:
      print("Launching Chimera visualization script...")
      status = os.system(cmd)
      if status:
        raise Exception("Chimera process return non-zero exit status.")
    except Exception as e:
      raise

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
    print(("Minimum score: %s"%min_score))
    print(("Maximum score: %s"%max_score))

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
  randlist = [random.randint(1,10) for i in range(digits)]
  multidigit_rand = int(''.join([str(x) for x in randlist]))
  return multidigit_rand

# Chimera executable
if __name__ == '__main__':
  from chimera import runCommand as rc
  from chimera import replyobj
  # args = [sys.argv.split(' ')]
  args = sys.argv
  params = {'structid':args[1],'biounit':int(args[2]),'anno':args[3],'attrf':args[4],
            'minval':float(args[5]),'maxval':float(args[6]),'struct_loc':args[7],
            'res_dir':args[8],'colors':args[9],'resis':[],'resrenumber':args[10],'vtype':args[11]}
  params["resrenumber"] = eval(params["resrenumber"])
  if params['colors'] != '-':
    params['colors'] = params['colors'].split(',')
  else:
    params['colors'] = None
  # Read the residues back into memory
  with open(params['attrf']) as fin:
    for line in fin:
      if line.startswith('\t'):
        params['resis'].append(line.strip().split('\t'))

  # Visualize with Chimera
  replyobj.status("Chimera: Visualizing %(structid)s_biounit%(biounit)d_%(anno)s_%(vtype)s"%params)
  # Assign the annotation to relevant residues
  rc("open %(struct_loc)s"%params)
  rc("defattr %(attrf)s raiseTool false"%params)
  # Initial colors
  rc("background solid white")
  rc("ribcolor dim grey")
  rc("ribbackbone")
  rc("~disp")
  rc("ribspline card smooth strand")
  for resi in params['resis']:
    # Identify all annotated residues as spheres
    rc("disp %s@ca"%resi[0])
    # And show atom/bonds around annotated residues
  # Define variant representation
  rc("represent bs")
  rc("setattr m autochain 0")
  rc("setattr m ballScale .7")
  if params['colors']:
    crange = list(range(int(params['minval']),int(params['maxval'])+1))
    colors = params['colors']
    print(';'.join(["%d->%s"%(val,colors[i]) for i,val in enumerate(crange)]))
    for i,val in enumerate(crange):
      rc("color %s,a :/%s=%d"%(colors[i],params['anno'],val))
  elif len(params['resis']) < 2 or params['minval'] == params['maxval']:
    rc("color red,a :/%s & @ca"%params['anno'])
  elif params['minval'] > params['maxval']:
    rc("rangecolor %(anno)s,a %(maxval)0.6f red %(minval)0.6f blue"%params)
    rc("color green,a :/%(anno)s=%(minval)s"%params)
    rc("color green,a :/%(anno)s>%(minval)s"%params)
    rc("color grey,a :/%(anno)s=%(maxval)s"%params)
    rc("color grey,a :/%(anno)s<%(maxval)s"%params)
  else:
    rc("rangecolor %(anno)s,a %(minval)0.6f blue %(maxval)0.6f red"%params)
    rc("color green,a :/%(anno)s=%(maxval)s"%params)
    rc("color green,a :/%(anno)s>%(maxval)s"%params)
    rc("color gray,a :/%(anno)s=%(minval)s"%params)
    rc("color gray,a :/%(anno)s<%(minval)s"%params)
  rc("transparency 70,r")
  # Export the scene
  rc("save %(res_dir)s/%(structid)s_biounit%(biounit)d_%(anno)s_%(vtype)s.py"%params)
  # Export the image
  rc("copy file %(res_dir)s/%(structid)s_biounit%(biounit)d_struct_%(anno)s_%(vtype)s.png width 5 height 5 units inches dpi 300"%params)
  # remove the structure, leaving only the variants
  rc("bonddisplay off")
  rc("copy file %(res_dir)s/%(structid)s_biounit%(biounit)d_vars_%(anno)s_%(vtype)s.png width 5 height 5 units inches dpi 300"%params)
  rc("close all")
  rc("stop now")
