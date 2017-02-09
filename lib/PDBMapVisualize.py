#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : PDBMapVisualize.py
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2017-02-09
# Description    : Visualization tool for PDBMap.
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv,time,platform

class PDBMapVisualize():
  
  def __init__(self,io,pdb_dir='data/rcsb',modbase_dir='data/modbase'):
    """ Initialize the PDBMapVisualize object. """
    self.io = io
    self.pdb_dir = pdb_dir
    self.modbase_dir = modbase_dir

  def visualize_structure(self,structid,biounit=0,anno_list=['maf'],spectrum_range=[],group=None,colors=[],permute=False,syn=False):
    """ Visualize the annotated dataset within a structure """
    biounit = 0 if biounit<0 else biounit
    print "Visualizing %s.%s..."%(structid,biounit)
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
        print "\nConverted annotation label from %s to %s"%(anno,convA)
      res[anno] = [1 if l==self.io.dlabel else None for l in res["dlabel"]]
      anno_list = [anno]

    # If synonymous variants are requested, switch the issnp flag
    if syn:
      print "\nSwitching from nonsynonymous to synonymous variants..."
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
            print "\nSwitching from nonsynonymous to synonymous variants..."
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
    print "(%d annotated residues)"%len(res['seqid'])
    print "(%d non-duplicate residues)"%len(set(zip(res['chain'],res['seqid'])))
    print "(%d with non-zero/null values)"%len([x for x in res[anno] if x])

    # Determine the first residue for renumbering
    resrenumber = {}
    for chain in set(res["chain"]):
      resrenumber[chain] = min(r for i,r in enumerate(res["aln_trans_seqid"]) if res["chain"][i]==chain)
    # Correct submodel ID for undivided structures
    maxm = len(set(res['model']))
    if maxm < 2:
      res['model'] = [0 for m in res['model']]

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
        struct_loc = "%s/models/model/%s.pdb.gz"%(self.modbase_dir,structid.upper())
      elif int(biounit) == 0:
        struct_loc = "%s/structures/all/pdb/pdb%s.ent.gz"%(self.pdb_dir,structid)
      else:
        struct_loc = "%s/biounit/coordinates/all/%s.pdb%s.gz"%(self.pdb_dir,structid,biounit)
      params = {'structid':structid,'biounit':biounit,'anno':anno,'attrf':attrf,'colors':colors,'vtype':params['vtype'],
                'minval':minval,'maxval':maxval,'resis':out,'struct_loc':struct_loc,'resrenumber':resrenumber}
      self.visualize(params,group=group)

  def visualize_unp(self,unpid,anno_list=['maf'],spectrum_range=[],colors=[],syn=False):
    """ Visualize the annotated dataset in all structures/model of the specified protein """
    print "Visualizing protein %s"%(unpid)
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
          print "Visualizing structure %s.%s"%(entity,biounit)
          self.visualize_structure(entity,biounit,anno_list,spectrum_range,group=unpid,colors=colors,syn=syn)
      elif entity_type == 'model':
        # Query all biological assemblies
        biounits = [-1]
        for biounit in biounits:
          print "Visualizing model %s.%s"%(entity,biounit)
          self.visualize_structure(entity,biounit,anno_list,spectrum_range,group=unpid,colors=colors,syn=syn)
      else:
        msg = "ERROR (PDBMapVisualize) Cannot visualize %s of type %s"%(entity,entity_type)
        raise Exception(msg)

  def visualize_all(self,anno_list=['maf'],spectrum_range=[],colors=[],syn=False):
    """ Visualize all structures and models annotated by the specified dataset """
    print "Visualizing dataset %s"%self.io.dlabel
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
        self.visualize_structure(s,b,anno_list,spectrum_range,group='all',colors=colors,syn=syn)
    models = [r[0] for r in res if self.io.detect_entity_type(r[0]) == 'model']
    for m in models:
      biounits = [-1]
      bres   = self.io.secure_query(query,(m,),cursorclass='Cursor')
      biounits = [r[0] for r in bres]
      for b in biounits:
        self.visualize_model(m,b,anno_list,spectrum_range,group='all',colors=colors,syn=syn)

  def visualize(self,params,group=None):
    """ Run the Chimera visualization with the specified parameters """
    # Setup the visualization output directory
    timestamp = str(time.strftime("%Y%m%d-%H"))
    res_dir = 'results/pdbmap_%s_%s_%s'
    if group:
      res_dir = res_dir%(self.io.dlabel,group,timestamp)
    else:
      res_dir = res_dir%(self.io.dlabel,params['structid'],timestamp)
    params['res_dir'] = res_dir
    if not os.path.exists(res_dir):
      os.system('mkdir -p %s'%res_dir)
    # Setup the chimera visualization parameters
    params['minval'] = "%0.6f"%params['minval']
    params['maxval'] = "%0.6f"%params['maxval']
    params['colors'] = '-' if not params['colors'] else ','.join(params['colors'])
    keys   = ['structid','biounit','anno','attrf','minval','maxval',
                'struct_loc','res_dir','colors','resrenumber','vtype']
    # Generate the subprocess command string
    script = '"lib/PDBMapVisualize.py %s"'%' '.join([str(params[key]).translate(None,' ') for key in keys])
    cmd    = "TEMP=$PYTHONPATH; unset PYTHONPATH; chimera --nogui --silent --script %s; export PYTHONPATH=$TEMP"%script
    if platform.system() == 'Darwin':
      # Allow Mac OSX to use the GUI window
      cmd  = "TEMP=$PYTHONPATH; unset PYTHONPATH; chimera --silent --script %s; export PYTHONPATH=$TEMP"%script
    try:
      status = os.system(cmd)
      if status:
        raise Exception("Chimera process return non-zero exit status.")
    except Exception as e:
      raise

# Chimera visualization executable; called as a subprocess
# by the class methods defined above
if __name__ == '__main__':
  from chimera import runCommand as rc
  from chimera import replyobj
  # Parse the command line arguments
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
  #rc("set bgTransparency")
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
    crange = range(int(params['minval']),int(params['maxval'])+1)
    colors = params['colors']
    print ';'.join(["%d->%s"%(val,colors[i]) for i,val in enumerate(crange)])
    for i,val in enumerate(crange):
      rc("color %s,a :/%s=%d"%(colors[i],params['anno'],val))
  elif len(params['resis']) < 2 or params['minval'] == params['maxval']:
    for resi in params['resis']:
      rc("color red,a %s@ca"%resi[0])
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
  rc("copy file %(res_dir)s/%(structid)s_biounit%(biounit)d_struct_%(anno)s_%(vtype)s.png width 1600 height 1600 units points dpi 300"%params)
  # remove the structure, leaving only the variants
  rc("bonddisplay off")
  rc("copy file %(res_dir)s/%(structid)s_biounit%(biounit)d_vars_%(anno)s_%(vtype)s.png width 1600 height 1600 units points dpi 300"%params)
  rc("close all")
  rc("stop now")