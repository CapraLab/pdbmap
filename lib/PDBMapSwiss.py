#!/usr/bin/env python
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

import sys,os,re
import pandas
from Bio.PDB.Structure import Structure
from lib.PDBMapProtein import PDBMapProtein
from lib.PDBMapAlignment import PDBMapAlignment
from collections import defaultdict

import logging
LOGGER = logging.getLogger(__name__)

class PDBMapSwiss(Structure):
  # SwissModels are represented in RAM by small dictionaries that are reached individually by modelid
  # or as lists from a Uniprot AC (unp)

  # This class contains only static functions, with the initial call of load_swiss_INDEX_JSON function
  # being a kind of global constructor to enable all the other functions

  # Maps Uniprot IDs to list of SwissModel modelIDs
  _unp2modelids = defaultdict(lambda: [])
  # Map the SwissModel modelIDs to small dictionaries that contain the 
  # key descriptors for each model.
  _modelid2info = {}

  # These are the column headings/json field names - for the INDEX_JSON file
  # We initally parse the small dictionaries with these key names, but then re-arrange the key names
  # to match ModBase nomenclature
  _JSON_fields  = ["uniprot_ac","template","from","to","provider","url","qmean","qmean_norm","coordinate_id","iso_id"]
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
      # print "Analyzing unp id", self.unp
      candidate_transcripts = PDBMapProtein.unp2enst(self.unp)
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
          LOGGER.warning("Unable to  cross reference transcripts from Ensembl for Uniprot AC: %s"%self.unp)
        if len(candidate_transcripts) > 1:
          candidate_transcripts = [candidate_transcripts[0]]
          LOGGER.info( "Too many transcripts for Uniprot isoform AC: %s.  Retained %s only"%(self.unp,candidate_transcripts[0].transcript))

      if len(candidate_transcripts) < 1:
        LOGGER.warning( "No transcripts from Ensembl SQL db found for Uniprot AC: %s"%self.unp)
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
  def ensp2swiss(cls,ensp):
    """ NOT IMPLEMENTED: Maps Ensembl protein IDs to swiss models """
    unps4ensps = PDBMapProtein.ensp2unp(ensp)
    # return PDBMapSwiss._modelid2info.get(PDBMapProtein.,None)
    models = PDBMapSwiss.swiss_dict.get(ensp,[])
    return models

  @classmethod
  def unp2modelids(cls,unp,IsoformOnly = False):
    """ Return all SwissModel ID strings for a Uniprot AC (unp) """
    modelids = []
    if (IsoformOnly):
      for modelid in PDBMapSwiss._unp2modelids[unp]:
        modelids.append(modelid)
    else:
      for modelid in PDBMapSwiss._unp2modelids[unp.split('-')[0]]:
        modelids.append(modelid)
    return modelids
      
  @classmethod
  def unp2swiss(cls,unp,IsoformOnly = False):
    """ Return all SwissModel dictionaries for a Uniprot AC (unp) """
    """ Maps a UniProt protein ID to list of SwissModel models """
    swiss = []
    for modelid in PDBMapSwiss.unp2modelids(unp,IsoformOnly):
      swiss.append(PDBMapSwiss.get_info(modelid))
    return swiss 

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
  def load_REMARK3_metrics(cls,modelid):
    """ Parse the REMARK section of a SwissModel .pdb file to get a variety of quality metrics (eg qmn4, sid etc) """
    metrics = {}
    fname = PDBMapSwiss.get_coord_file(modelid)
    with open(fname,'rt') as f:
      for line in f:
        if line.startswith('ATOM  '): # Then we're done parsing the pdb file.  ATOM comes after all REMARK lines of interest
          break;
        if line.startswith('REMARK   3') and len(line) > 20: # Then it's likely a swiss-specific model REMARK
          splits = re.split(' +',line[:20].rstrip())
          if len(splits) == 3:
            key = splits[2].lower() # Examples: enging, gmqe, qmn4, sid
            value = line[20:].rstrip()
            # 'template' conflicts with the JSON meta file entry of same name, and carries less info
            if (key != 'template'):
              metrics[key] = value
    return metrics


  @classmethod
  def get_unp(cls,modelid):
    """ Returns the coordinate file location for the SwissModel model ID """
    if not PDBMapSwiss._modelid2info:
      raise Exception("PDBMapSwiss.load_swiss_INDEX_JSON must be called before using this method.")
    return PDBMapSwiss._modelid2info[modelid]['unp']

  @classmethod
  def concatenate_modelid(cls,unp,start,end,template):
    """ Build a Swiss-model unique key from INDEX file row information """
    return "%s_%s_%s_%s"%(unp,start,end,template)

  @classmethod
  def load_swiss_INDEX_JSON(cls,swiss_dir,summary_fname):
    """ Adds a SwissModel summary file to the lookup dictionary """
    # print "Parameter swiss_dir = " + swiss_dir
    # print "Parameter summary_fname = " + summary_fname
    swiss_dir = swiss_dir.rstrip('/')
    if not PDBMapSwiss.swiss_dir:
      PDBMapSwiss.swiss_dir = swiss_dir

    if not os.path.exists(summary_fname):
      msg = "ERROR: (PDBMapSwiss) Cannot load SwissModel. %s does not exist."%summary_fname
      logging.exception(msg)
      raise Exception
    logging.getLogger().info("Opening SwissModel Summary File (JSON format):" + summary_fname)
    import json
    with open(summary_fname,'rt') as json_file:
      json_str = json_file.read()
      json_data = json.loads(json_str)
      del json_str

    logging.getLogger().info("%d rows read from %s successfully."%(len(json_data),summary_fname))

    for d in json_data:
      # import pdb; pdb.set_trace()
      # For sanity, we need dictionary elements that match the naming convention of our ModBase models
      # Convert columns from/to to start/end  notation like ModBase
      d['start'] = d.pop('from')
      d['end'] = d.pop('to')

      # Convert column 'uniprot_ac' to 'unp' for modbase compatability
      # First, if there is an isoform identifier which does not match the base iso-less identifier
      # Then be sure to use the new isoform identifier.
      if ('iso_id' in d and
          ('-' not in d['uniprot_ac']) and 
          ('-' in d['iso_id']) and
          (d['uniprot_ac'] != d['iso_id'])):
        d['unp'] = d['iso_id']
      else:
        d['unp'] = d['uniprot_ac']

      del d['uniprot_ac']
      if 'iso_id' in d:
        del d['iso_id']

      d['modelid'] = PDBMapSwiss.concatenate_modelid(d['unp'],d['start'],d['end'],d['template'])

      # We have a lot of duplicates - don't worry about that for now (Chris Moth 2017-09-11
      # if (modelid in PDBMapSwiss._modelid2info):
      # print "CAUTION: %s aleady in swiss dictionary\n"%modelid
      # Add this unique modelid to the growing large dictionary      
      PDBMapSwiss._modelid2info[d['modelid']] = d
      PDBMapSwiss._unp2modelids[d['unp']].append(d['modelid'])
      # If this unp is an isoform, then also add it to the generic (base) unp dictionary
      if (d['unp'].find('-') != -1):
        PDBMapSwiss._unp2modelids[d['unp'].split('-')[0]].append(d['modelid'])

    LOGGER.info( "%d Swiss models added to in-memory dictionary"%len(PDBMapSwiss._modelid2info))

  @classmethod
  def old_load_swiss_INDEX_JSON(cls,swiss_dir,summary_fname):
    """ Adds a SwissModel summary file to the lookup dictionary """
    # print "Parameter swiss_dir = " + swiss_dir
    # print "Parameter summary_fname = " + summary_fname
    swiss_dir = swiss_dir.rstrip('/')
    if not PDBMapSwiss.swiss_dir:
      PDBMapSwiss.swiss_dir = swiss_dir

    if not os.path.exists(summary_fname):
      msg = "ERROR: (PDBMapSwiss) Cannot load SwissModel. %s does not exist."%summary_fname
      logging.exception(msg)
      raise Exception
    logging.getLogger().info("Opening SwissModel Summary File (JSON format):" + summary_fname)
    try:
      df = pandas.read_json(summary_fname) # ,orient='records')
    except ValueError:
      msg = "Failure to read JSON file %s\nProblem is likely ill-formatted data in the JSON file.  Read error carefully:"%summary_fname
      logging.exception(msg)
      raise
    except:
      msg =  "Failure attempting read of JSON SwissModel Summary File  Additional details follow:"
      logging.exception(msg)
      raise

    logging.getLogger().info("%d rows read from %s successfully."%(len(df),summary_fname))
    total_count = 0
    for index,row in df.iterrows():
      total_count += 1
      try:
        d = dict([(dname,row[dname]) for dname in PDBMapSwiss._JSON_fields if pandas.notnull(row[dname])])
      except:
        msg = "Unable to map expected fields on JSON file row %d"%total_count
        msg += "\nRead %s for this row:\n"%row
        logging.exception(msg)
        raise

      # import pdb; pdb.set_trace()
      # For sanity, we need dictionary elements that match the naming convention of our ModBase models
      # Convert columns from/to to start/end  notation like ModBase
      d['start'] = d['from']
      d['end'] = d['to']
      del d['from']
      del d['to']

      # Convert column 'uniprot_ac' to 'unp' for modbase compatability
      # First, if there is an isoform identifier which does not match the base iso-less identifier
      # Then be sure to use the new isoform identifier.
      if ('iso_id' in d and
          ('-' not in d['uniprot_ac']) and 
          ('-' in d['iso_id']) and
          (d['uniprot_ac'] != d['iso_id'])):
        d['unp'] = d['iso_id']
      else:
        d['unp'] = d['uniprot_ac']

      del d['uniprot_ac']
      if 'iso_id' in d:
        del d['iso_id']

      modelid = PDBMapSwiss.concatenate_modelid(d['unp'],d['start'],d['end'],d['template'])

      d['modelid'] = modelid

      # We have a lot of duplicates - don't worry about that for now (Chris Moth 2017-09-11
      # if (modelid in PDBMapSwiss._modelid2info):
      # print "CAUTION: %s aleady in swiss dictionary\n"%modelid
      # Add this unique modelid to the growing large dictionary      
      PDBMapSwiss._modelid2info[d['modelid']] = d
      PDBMapSwiss._unp2modelids[d['unp']].append(d['modelid'])
      # If this unp is an isoform, then also add it to the generic (base) unp dictionary
      if (d['unp'].find('-') != -1):
        PDBMapSwiss._unp2modelids[d['unp'].split('-')[0]].append(d['modelid'])

    LOGGER.info( "%d Swiss models added to in-memory dictionary"%len(PDBMapSwiss._modelid2info))

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
