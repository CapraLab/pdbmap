#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : PDBMapModel.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-04-01
# Description    : PDBMapStructure equivalent for ModBase models
#=============================================================================#
# PBS Parameters
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l nodes=vision.mc.vanderbilt.edu
#PBS -l mem=5000mb
#PBS -l walltime=1:00:00:00
#=============================================================================#

# See main check for cmd line parsing
import argparse
import sys,os,csv
from Bio.PDB.Structure import Structure
from lib.PDBMapProtein import PDBMapProtein
from lib.PDBMapTranscript import PDBMapTranscript
from lib.PDBMapAlignment import PDBMapAlignment

class PDBMapModel(Structure):
  
  # Modbase Summary Dictionary
  # Maps Ensembl Protein IDs to Modbase Models
  modbase_dir   = []
  modbase_dict  = {}
  _modelid2info = {}
  _modelid2file = {}
  _info_fields  = ["run_name","mod_seq_id","mod_model_id","modelid","start","end",
                    "sequence_identity","evalue","ga341","mpqs","zdope","template_pdb",
                    "template_chain","template_start","template_end","hit_history",
                    "tsvmod_method","tsvmod_no35","tsvmod_rmsd","filename","unp"]

  def _sfloat(self,obj):
    """ Safe float conversion """
    try: return float(obj)
    except: return "'NULL'" # for MySQL upload

  def _sint(self,obj):
    """ Safe int conversion """
    try: return int(obj)
    except: return "'NULL'" # for MySQL upload

  def __init__(self,s,model_summary,quality=-1):
    # Record the chain information
    chain    = [c for c in s.get_chains()][0]
    oid = chain.id
    chain.id = 'A'
    chain.species = "HUMAN"
    chain.pdbstart = self._sint(model_summary['start'])
    chain.pdbend   = self._sint(model_summary['end'])
    chain.unp      = model_summary['unp']
    chain.offset   = 0 # Assumed. Corrected by alignment.
    chain.hybrid   = 0
    s[0].child_dict['A'] = chain
    del s[0].child_dict[oid]

    # Store the edited structure
    self.structure   = s
    self.quality     = quality
    self.transcripts = []
    self.alignments  = []

    # Store the model summary information
    self.id      = model_summary['modelid']
    self.tvsmod_method = model_summary['tvsmod_method']
    if self.tvsmod_method != 'NA':
      self.tvsmod_no35   = self._sfloat(model_summary['tvsmod_no35'])
      self.tvsmod_rmsd   = self._sfloat(model_summary['tvsmod_rmsd'])
    else:
      self.tvsmod_no35 = 'NULL'
      self.tvsmod_rmsd = 'NULL'
    self.identity = self._sfloat(model_summary['seequence_identity'])
    self.evalue   = self._sfloat(model_summary['evalue'])
    self.ga341    = self._sfloat(model_summary['ga341'])
    self.mpqs     = self._sfloat(model_summary['mpqs'])
    self.zdope    = self._sfloat(model_summary['zdope'])
    self.pdbid    = model_summary['template_pdb']
    self.chain    = model_summary['template_chain']
    self.unp      = model_summary['unp']
    
  def __getattr__(self,attr):
    # Defer appropriate calls to the structure
    if attr in dir(self.structure):
      result = self.structure.__getattribute__(attr)
    else:
      result = self.__getattribute__(attr)
    if callable(result):
      def hooked(*args, **kwargs):
        result = result(*args,**kwargs)
        if result == self.structure:
          return self
        return result
      return hooked
    else:
      return result

  def get_transcripts(self,io=None):
    # io is an used parameter required for polymorphic behavior
    # Retrieve the corresponding transcript for each chain
    if self.transcripts:
      return self.transcripts
    for chain in self.structure[0]:
      # Query all transcripts associated with the chain's UNP ID
      candidate_transcripts = PDBMapTranscript.query_from_unp(self.unp)
      # But only keep the one matching this model's reference ENSP, if specified
      if self.id.startswith("ENSP"):
        candidate_transcripts = [ct for ct in candidate_transcripts if
                               PDBMapProtein.enst2ensp(ct.transcript)==self.id.split('.')[0].split('_')[0]]
      if len(candidate_transcripts) < 1:
        return []
      if len(candidate_transcripts) > 1 and self.id.startswith("ENSP"):
        msg = "WARNING (PDBMapModel) Too many transcripts for %s. Truncating.\n"%self.id
        sys.stderr.write(msg)
        candidate_transcripts = [candidate_transcripts[0]]
      # Align chain to first candidate transcript
      alignments = [PDBMapAlignment(chain,candidate_transcripts[0])]
      # Repeat for remaining chains, select highest scoring alignment
      for trans in candidate_transcripts[1:]:
        new_alignment = PDBMapAlignment(chain,trans)
        alignments.append(new_alignment)
      # Store best transcript alignment as element of chain
      chain.alignments  = alignments
      chain.transcripts = [a.transcript for a in alignments]
    # Return the matched transcripts
    self.transcripts = [t for c in self.structure[0] for t in c.transcripts]
    self.alignments  = [a for c in self.structure[0] for a in c.alignments]
    return self.transcripts

  def get_alignments(self):
    if not self.alignments:
      self.get_transcripts()
    else:
      return self.alignments

  @classmethod
  def get_models(cls):
    """ Returns all recorded ModBase models """
    return [m for v in PDBMapModel.modbase_dict.values() for m in v]

  @classmethod
  def ensp2modbase(cls,ensp):
    """ Maps Ensembl protein IDs to ModBase models """
    models = PDBMapModel.modbase_dict.get(ensp,[])
    return models

  @classmethod
  def unp2modbase(cls,unp):
    """ Maps UniProt protein IDs to ModBase models """
    # Get all matching Ensembl protein IDs
    ensps = PDBMapProtein.unp2ensp(unp)
    models = []
    for ensp in ensps:
      # Get all matching ModBase models
      models.extend(PDBMapModel.ensp2modbase(ensp))
    models = [model for model in models if model]
    return models

  @classmethod
  def get_info(cls,modelid):
    """ Returns the info dictionary for this ModBase model """
    if not PDBMapModel.modbase_dict:
      raise Exception("PDBMapModel.load_modbase must be called before using this method.")
    return PDBMapModel.modbase_dict.get(modelid,None)

  @classmethod
  def get_coord_file(cls,modelid):
    """ Returns the coordinate file location for the ModBase model ID """
    if not PDBMapModel._modelid2file:
      raise Exception("PDBMapModel.load_modbase must be called before using this method.")
    return PDBMapModel._modelid2file.get(modelid,None)

  @classmethod
  def load_modbase(cls,modbase_dir,summary_fname):
    """ Adds a ModBase summary file to the lookup dictionary """
    modbase_dir = modbase_dir.rstrip('/')
    if not PDBMapModel.modbase_dir:
      PDBMapModel.modbase_dir = [modbase_dir]
    else:
      PDBMapModel.modbase_dir.append(modbase_dir)
    if not os.path.exists(summary_fname):
      msg = "ERROR: (PDBMapModel) Cannot load ModBase. %s does not exist."%summary_fname
      raise(Exception(msg))
    fin = open(summary_fname,'rb')
    fin.readline() # burn the header
    reader = csv.reader(fin,delimiter='\t')
    total_count = 0
    mismatch_count = 0
    for row in reader:
      # If row length matches 2013 summary, reconcile columns with 2016
      if len(row) != len(PDBMapModel._info_fields)-2:
        # Insert two dummy rows between columns 1 and 2 to match 2016 summary format
        row.insert(1,None)
        row.insert(1,None)
      # Skip models not constructed from EnsEMBL protein sequences
      if not row[3].startswith("ENSP"): continue
      # Append the coordinate file location to the summary info
      row.append("%s/%s.pdb.gz"%(modbase_dir,row[3]))
      # Extract the EnsEMBL protein identifier
      ensp = row[3].split('.')[0] # Strip the EnsEMBL protein verion number, if present
      ensp = ensp.split('_')[0]   # Strip the ModBase model iteration number
      unps = PDBMapProtein.ensp2unp(ensp)
      unp  = None if not unps else unps[0]
      total_count += 1
      if not unp:
        mismatch_count += 1
      # Skip models without associated UniProt IDs
      if not unp: continue
      # Set UniProt ID as last field in model summary
      row.append(unp)
      # Convert to dictionary
      d = dict(zip(PDBMapModel._info_fields,row))
      # Populate Ensembl protein ID to model summary dictionary
      if ensp in PDBMapModel.modbase_dict:
        PDBMapModel.modbase_dict[ensp].append(d['modelid'])
      else:
        PDBMapModel.modbase_dict[ensp] = [d['modelid']]
      # Populate the ModBase model ID to summary info dictionary
      PDBMapModel._modelid2info[d['modelid']] = d
      # Populate ModBase model ID to coordinate file dictionary
      PDBMapModel._modelid2file[d['modelid']] = d['filename']
  
# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
