#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : PDBMapSwiss.py
# Author         : Chris Moth
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2017-09-07
# Description    : PDBMapStructure equivalent for SwissModel models
#                  based on Mike Sivley's PDBMapModel.py

import argparse
# import json

import sys,os
import pandas
from Bio.PDB.Structure import Structure
from lib.PDBMapProtein import PDBMapProtein
from lib.PDBMapTranscript import PDBMapTranscript
from lib.PDBMapAlignment import PDBMapAlignment

import logging
logging.basicConfig(Level='INFO')

class PDBMapSwiss(Structure):
  
  # SwissModel Summary Dictionary
  # Maps Uniprot IDs to SwissModel Swisss
  _modelid2info = {}
  # These are the column headings/json field names - for the INDEX_JSON file
  _JSON_fields  = ["uniprot_ac","template","from","to","provider","url","qmean","qmean_norm","coordinate_id"]
  swiss_dir = None

  '''
Fields in modbase that may become more relevant
"run_name","mod_seq_id","mod_model_id","modelid","start","end","sequence_identity","evalue","ga341","mpqs","zdope","template_pdb", "template_chain","template_start","template_end","hit_history", "tsvmod_method","tsvmod_no35","tsvmod_rmsd","filename","unp"]
  '''


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
    self.template   = model_summary['template']
    self.unp      = model_summary['unp']
    self.start = model_summary['start']
    self.end = model_summary['end']
    self.qmean      = model_summary['qmean']
    self.qmean_norm      = model_summary['qmean_norm']
    self.coordinate_id      = model_summary['coordinate_id']
    self.url      = model_summary['url']
    
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
    # import pdb; pdb.set_trace()
    # io is an unused parameter required for polymorphic behavior
    # Retrieve the corresponding transcript for each chain
    if self.transcripts:
      return self.transcripts
    for chain in self.structure[0]:
      # Query all transcripts associated with the chain's UNP ID
      # from the Ensembl mySQL database records
      print "Analyzing unp id", self.unp
      candidate_transcripts = PDBMapTranscript.query_from_unp(self.unp)
      # But only keep the Ensemble transcripts matching this model's reference ENSP, if specified
      # IF the uniporit ID has a DASH (only!) THEN make sure to ONLY keep
      # The transcript which we KNOW from INDEX_JSON and tracing... to be the relevant one
      unp_split = self.unp.split("-",2)
      if (len(unp_split) > 1):
        # This is a sanity check that is similar to the one in PDBMapModel in idea
        # We are making sure that all the candidate_transcripts returned from the 
        # (big) Ensembl SQL query map to an ENSP... protein (...2ensp) that is 
        # also in the PDBMapProtine mapping file.  If not, we have a fundamental problem
        candidate_transcripts = [ct for ct in candidate_transcripts if
                               PDBMapProtein.enst2ensp(ct.transcript) in PDBMapProtein.unp2ensp(self.unp)]
        if len(candidate_transcripts) < 1:
          logging.getLogger(__name__).warning("Unable to  cross reference transcripts from Ensembl for Uniprot AC: %s"%self.unp)
        if len(candidate_transcripts) > 1:
          candidate_transcripts = [candidate_transcripts[0]]
          logging.getLogger(__name__).info( "Too many transcripts for Uniprot isoform AC: %s.  Retained %s only"%(self.unp,candidate_transcripts[0].transcript))

      if len(candidate_transcripts) < 1:
        logging.getLogger(__name__).warning( "No transcripts from Ensembl SQL db found for Uniprot AC: %s"%self.unp)
        return []

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

    return self.alignments

  @classmethod
  def get_swiss_modelids(cls):
    """ Returns all recorded SwissModel models """
    return list(PDBMapSwiss._modelid2info)  

  @classmethod
  def unp2swiss(cls,unp):
    """ Maps UniProt protein IDs to SwissModel models """
    # Get all matching Ensembl protein IDs
    ensps = PDBMapProtein.unp2ensp(unp)
    models = []
    for ensp in ensps:
      # Get all matching SwissModel models
      models.extend(PDBMapSwiss.ensp2swiss(ensp))
    models = [model for model in models if model]
    return models

  @classmethod
  def get_info(cls,modelid):
    """ Returns the info dictionary for this SwissModel model """
    if not PDBMapSwiss._modelid2info:
      raise Exception("PDBMapSwiss.load_swiss_INDEX_JSON must be called before using this method.")
    return PDBMapSwiss._modelid2info.get(modelid,None)

  @classmethod
  def get_coord_file(cls,modelid):
    """ Returns the coordinate file location for a SwissModel model ID """
    if not PDBMapSwiss._modelid2info:
      raise Exception("PDBMapSwiss.load_swiss_INDEX_JSON must be called before using this method.")
    d = PDBMapSwiss._modelid2info[modelid]
    uniprot_directories = d['unp']
    dash_or_end = uniprot_directories.rfind('-')
    # In the SwissModel filename system
    # The Uniprot Isoform is dropped to locate the directory
    # However,the coordinate_id quasi-guid could vary per isoform
    if (dash_or_end == -1):
      dash_or_end = len(uniprot_directories)
    return "%s/%s/%s/%s/swissmodel/%d_%d_%s_%s.pdb"%(
		PDBMapSwiss.swiss_dir,
		d['unp'][0:2],
		d['unp'][2:4],
		d['unp'][4:dash_or_end],
		# d['unp'][4:dash_or_end],
      int(d['start']),
      int(d['end']),
      d['template'],
      d['coordinate_id']   # 
    ) 

  @classmethod
  def get_unp(cls,modelid):
    """ Returns the coordinate file location for the SwissModel model ID """
    if not PDBMapSwiss._modelid2info:
      raise Exception("PDBMapSwiss.load_swiss_INDEX_JSON must be called before using this method.")
    return PDBMapSwiss._modelid2info[modelid]['unp']

  @classmethod
  def concatenate_modelid(cls,unp,start,end,template):
    """ Build a Swiss-model unique key from INDEX file row information """
    return "%s_%d_%d_%s"%(unp,int(start),int(end),template)

  @classmethod
  def load_swiss_INDEX_JSON(cls,swiss_dir,summary_fname):
    """ Adds a SwissModel summary file to the lookup dictionary """
    print "Parameter swiss_dir = " + swiss_dir
    print "Parameter summary_fname = " + summary_fname
    swiss_dir = swiss_dir.rstrip('/')
    if not PDBMapSwiss.swiss_dir:
      PDBMapSwiss.swiss_dir = swiss_dir

    if not os.path.exists(summary_fname):
      msg = "ERROR: (PDBMapSwiss) Cannot load SwissModel. %s does not exist."%summary_fname
      print msg
      raise(Exception(msg))
    print "Opening SwissModel Summary File (JSON format):" + summary_fname
    try:
      df = pandas.read_json(summary_fname) # ,orient='records')
    except ValueError:
      print "Failure to read JSON file %s\nProblem is likely ill-formatted data in the JSON file.  Read error carefully:"%summary_fname
      raise
    except:
      print "Failure attempting read of JSON SwissModel Summary File  Additional details follow:"
      raise

    print "%d rows read from %s successfully."%(len(df),summary_fname)
    total_count = 0
    for index,row in df.iterrows():
      total_count += 1
      try:
        d = dict([(dname,row[dname]) for dname in PDBMapSwiss._JSON_fields])
      except:
        print "Unable to map expected fields on JSON file row %d"%total_count
        print "Read %s for this row:\n"%row
        raise

      # import pdb; pdb.set_trace()
      # For sanity, we need dictionary elements that match the naming convention of our ModBase models
      # Convert columns from/to to start/end  notation like ModBase
      d['start'] = d['from']
      d['end'] = d['to']
      del d['from']
      del d['to']


      # Convert column 'uniprot_ac' to 'unp' for modbase compatability
      d['unp'] = d['uniprot_ac']
      del d['uniprot_ac']
      modelid = PDBMapSwiss.concatenate_modelid(d['unp'],d['start'],d['end'],d['template'])

      d['modelid'] = modelid

      # We have a lot of duplicates - don't worry about that for now (Chris Moth 2017-09-11
      # if (modelid in PDBMapSwiss._modelid2info):
        # print "CAUTION: %s aleady in swiss dictionary\n"%modelid
      # Add this unique modelid to the growing large dictionary      
      PDBMapSwiss._modelid2info[d['modelid']] = d
    logging.getLogger(__name__).info( "%d Swiss models added to in-memory dictionary"%len(PDBMapSwiss._modelid2info))

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
