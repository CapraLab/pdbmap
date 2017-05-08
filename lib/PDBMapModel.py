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
  modbase_dir  = ''
  modbase_dict = {}

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
    chain.pdbstart = self._sint(model_summary[4])
    chain.pdbend   = self._sint(model_summary[5])
    chain.unp      = model_summary[19]
    chain.offset   = 0 # Assumed. Corrected by alignment.
    chain.hybrid   = 0
    s[0].child_dict['A'] = chain
    del s[0].child_dict[oid]

    # Store the edited structure
    self.structure   = s
    self.quality     = quality
    self.transcripts = []

    # Store the model summary information
    self.id      = model_summary[3]
    self.tvsmod_method = model_summary[16]
    if self.tvsmod_method != 'NA':
      self.tvsmod_no35   = self._sfloat(model_summary[17])
      self.tvsmod_rmsd   = self._sfloat(model_summary[18])
    else:
      self.tvsmod_no35 = 'NULL'
      self.tvsmod_rmsd = 'NULL'
    self.identity = self._sfloat(model_summary[6])
    self.evalue   = self._sfloat(model_summary[7])
    self.ga341    = self._sfloat(model_summary[8])
    self.mpqs     = self._sfloat(model_summary[9])
    self.zdope    = self._sfloat(model_summary[10])
    self.pdbid    = model_summary[11]
    self.chain    = model_summary[12]
    self.unp      = model_summary[17]
    
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
                               PDBMapProtein.enst2ensp(ct.transcript)==self.id.split('.')[0]]
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
  def load_modbase(cls,modbase_dir,summary_fname):
    """ Loads a ModBase summary file into a lookup dictionary """
    PDBMapModel.modbase_dir = modbase_dir
    summary_path = "%s/%s"%(modbase_dir,summary_fname)
    if not os.path.exists(summary_path):
      msg = "ERROR: (PDBMapModel) Cannot load ModBase. %s does not exist."%summary_path
      raise(Exception(msg))
    fin = open(summary_path,'rb')
    fin.readline() # burn the header
    reader = csv.reader(fin,delimiter='\t')
    for row in reader:
      # Skip models not constructed from EnsEMBL protein sequences
      if not row[3].startswith("ENSP"): continue
      # Extract the EnsEMBL protein identifier
      ensp = row[3].split('.')[0] # Strip the EnsEMBL protein verion number
      unps = PDBMapProtein.ensp2unp(ensp)
      unp  = None if not unps else unps[0]
      # Skip models without associated UniProt IDs
      if not unp: continue
      # Set UniProt ID as last field in model summary
      row.append(unp)
      if ensp in PDBMapModel.modbase_dict:
        PDBMapModel.modbase_dict[ensp].append(row)
      else:
        PDBMapModel.modbase_dict[ensp] = [row]
  
# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
