#!/usr/bin/env python27
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

  def __init__(self,s,model_summary):
    # Store the structure
    self.structure = structure
    # Correct some fields
    self.structure[0][' '].id      = 'A'
    self.structure[0]['A'].species = "HUMAN"
    # Store the model summary information
    self.id = model_summary[1]
    self.evalue  = model_summary[4]
    self.ga341   = model_summary[5]
    self.mpqs    = model_summary[6]
    self.zdope   = model_summary[7]
    self.pdbid   = model_summary[8]
    self.chain   = model_summary[9]
    self.tvsmod_method = model_summary[11]
    self.tvsmod_no35   = model_summary[15]
    self.tvsmod_rmsd   = model_summary[16]
    self.unp     = model_summary[17]

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

  def get_transcripts(self):
    # Retrieve the corresponding transcript for each chain
    if self.transcripts:
      return self.transcripts
    for chain in self.structure[0]:
      # Query all transcripts associated with the chain's UNP ID
      candidate_transcripts = PDBMapTranscript.query_from_unp(chain.unp)
      if len(candidate_transcripts) < 1:
        return []
      #UPDATE: Keep all transcript matches
      # Align chain to first candidate transcript
      alignments = [PDBMapAlignment(chain,candidate_transcripts[0])]
      # Repeat for remaining chains, select highest scoring alignment
      for trans in candidate_transcripts[1:]:
        new_alignment = PDBMapAlignment(chain,trans)
        alignments.append(new_alignment)
        # Determine best alignment
        # if new_alignment.score > alignment.score:
        #   alignment = new_alignment
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
    return PDBMapModel.modbase_dict.get(ensp,None)

  @classmethod
  def unp2modbase(cls,unp):
    """ Maps UniProt protein IDs to ModBase models """
    return PDBMapModel.ensp2modbase(PDBMapProtein.unp2ensp(unp))

  @classmethod
  def load_modbase(cls,modbase_dir,summary_fname):
    """ Loads a ModBase summary file into a lookup dictionary """
    PDBMapModel.modbase_dir = modbase_dir
    summary_path = "%s/%s"%(modbase_dir,summary_fname)
    if not os.path.exists(summary_path):
      msg = "ERROR: (PDBMapModel) Cannot load ModBase. %s does not exist."%summary_path
      raise(msg)
    fin = open(summary_path,'rb')
    fin.readline() # burn the header
    reader = csv.reader(fin,delimiter='\t')
    for row in reader:
      row = line.strip().split('\t')
      ensp,i = row[1].split('_')
      unp = PDBMapProtein.ensp2unp(ensp)
      row.append(unp) # set UniProt ID as last field in model summary
      if ensp in PDBMapModel.modbase_dict:
        PDBMapModel.modbase_dict[unp].append(row)
      else:
        PDBMapModel.modbase_dict[unp] = [row]
    #DEBUG
    print "ModBase Summary Dictionary:"
    for key,val in PDBMapModel.modbase_dict.iteritems():
      print key,val
  

# I am just saving these here for the time being, decide where each
# piece belongs once things are worked out.

# Load the ModBase homology index
modbase_dir = '/projects/Bush_eQTL/sivleyrm/data/modbase/H_sapiens_2013/'
homo_dir   = "%s/%s"%(modbase_dir,'ModBase_H_sapiens_2013_GRCh37.70.pep.all')
homo_index = "%s/%s"%(homo_dir,'H_sapiens_2013_GRCh37.70.pep.all.summary.txt')
fin = open(homo_index,'rb')
fin.readline()
homo_dict = {}
for line in fin:
  row    = line.strip().split('\t')
  ensp,i = row[1].split('_')
  if ensp in homo_dict:
    homo_dict[ensp].append(row)
  else:
    homo_dict[ensp] = [row]

# Query an Ensembl protein and get the associated ModBase models
test_pid   = "ENSP00000434723"
model_summaries = homo_dict[test_pid]
model_dir  = "%s/%s"%(homo_dir,'models/model')
for summary in model_summaries:
  id = summary[1]
  model_fname = "%s/%s.pdb"%(model_dir,id)
  if not os.path.exists(model_fname):
    model_fname += ".xz"
  if not os.path.exists(model_fname):
    raise Exception("Model not found.")
# Load the structure in model_fname
# Use get_structure as a template
# Remove the section checking for DBREF
# Add a few lines to read the MODPIPE MODEL ID from REMARK 220,
#   the MODELLER objective function from REMARK 6, and the
#   MODPIPE version number from REMARK 6
# Import the UNP from original query (not ENSP)
# Import the chain ID as 'A'
# Infer start/end from residue values
# Save all of this in a get_model function in PDBParser
# Use it to create an object of type PDBMapModel
# Include metadata/summary info about the ModBase model
# Create new database schema for PDBMapModel
# Upload the models
# Intersect with GenomicConsequence
  

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
