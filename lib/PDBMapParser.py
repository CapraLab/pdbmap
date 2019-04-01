#!/usr/bin/env python
#
# Project        : PDBMap
# Filename       : PDBMapParser.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-09
# Description    : PDBParser class utilizing Bio.PDB.PDBParser and mysql to
#                : read and upload structures to the PDBMap.Structure database.
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv,collections,gzip,time,random
import logging

import subprocess as sp
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
import Bio.PDB
import numpy as np
from lib.PDBMapModel import PDBMapModel
from lib.PDBMapProtein import PDBMapProtein
from lib.PDBMapStructure import PDBMapStructure
from lib.PDBMapIO import PDBMapIO,aa_code_map
import MySQLdb, MySQLdb.cursors
from warnings import filterwarnings,resetwarnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

LOGGER = logging.getLogger(__name__)

class PDBMapParser(PDBParser):
  def __init__(self,PERMISSIVE=True,get_header=True,
                structure_builder=None,QUIET=True):
    super(PDBMapParser,self).__init__(PERMISSIVE,get_header,
                                      structure_builder,QUIET)

  @classmethod
  def process_structure(cls,s,biounit=0,dbref={},force=False):
    """ The asymmetric unit *must* be processed before any
        biological assemblies, or the additional models
        containing the biological assembly coordinates will
        be removed. """
    # print "   # Processing biological assembly %d"%biounit
    # Process structural elements
    iter_s = [m for m in s]
    if not biounit:
      for m in iter_s[1:]:
        s.detach_child(m.id) # Retain only the first model
    # Record the biological assembly ID for each entity
    for m in s:
      m.biounit = biounit
      for c in m:
        c.biounit = biounit
        for r in m:
          r.biounit = biounit
          for a in r:
            a.biounit = biounit
    iter_s = [m for m in s]
    for m in iter_s:
      iter_m = [c for c in m] # avoid modification during iteration, shallow
      for c in iter_m:
        if dbref:
          if c.id not in dbref:
            m.detach_child(c.id)
            continue # chain not human, or other exclusion criteria
          # Apply DBREF fields from asymmetric unit to biological assembly
          c.unp      = dbref[c.id]['unp']
          c.gene     = dbref[c.id]['gene']
          c.pdbstart = dbref[c.id]['pdbstart']
          c.pdbend   = dbref[c.id]['pdbend']
          c.offset   = dbref[c.id]['offset']
          c.species  = dbref[c.id]['species']
          c.hybrid   = dbref[c.id]['hybrid']
        if 'species' not in dir(c):
          c.species = 'UNKNOWN'
        if c.species != 'HUMAN':
          if force:
            c.species = "HUMAN"
          else:
            # msg = "WARNING (PDBMapIO) Ignoring non-human chain: %s.%s (%s)\n"%(s.id,c.id,c.species)
            # sys.stderr.write(msg)
            m.detach_child(c.id)
            continue
        iter_c = [r for r in c] # avoid modification during iteration, shallow
        for r in iter_c:
          if dbref and r in dbref[c.id]['drop_resis']: # If residue outside dbref range
            c.detach_child(r.id) # Remove residue
          if r.id[0].strip():    # If residue is a heteroatom
            c.detach_child(r.id) # Remove residue
          elif not force and (r.id[1] < c.pdbstart or r.id[1] > c.pdbend): # If residue outside chain boundary
            c.detach_child(r.id)
          elif r.id[1] < 1: # Residue index must be non-negative
            c.detach_child(r.id)
          else:
            # Assign a 1-letter amino acid code
            if 'resname' not in dir(r) or r.resname.lower() not in aa_code_map:
              r.resname='SER' # if unknown, dummy code serine for alignment
            r.rescode = aa_code_map[r.resname.lower()]
            # Compute the center of mass for all residues
            natom = len([a for a in r])
            r.coord  = sum([a.coord for a in r]) / natom
            # Save the structural coordinate independently
            r.x,r.y,r.z = r.coord
            # Save the sequence coordinates independently
            _,r.seqid,r.icode = r.id
        if not len(c): # If chain contained only heteroatoms
          m.detach_child(c.id) # Remove chain
        else:
          # Parse the chain sequence and store as string within the chain
          c.sequence = ''.join([r.rescode for r in c])
      if not len(m): # If the model only contained non-human species
        s.detach_child(m.id)
    if not len(s):
      msg = "ERROR (PDBMapIO) %s contains no human protein chains."%s.id
      raise Exception(msg)
    return s

  @classmethod
  def getBiopythonStructureOrFail(cls,modelid,fname):
    functionNameAsString = sys._getframe().f_code.co_name
    ext = os.path.basename(fname).split('.')[-1]
    try:
      if ext == 'gz':
        fin = gzip.open(fname,'rt')
      elif ext in ['txt','pdb','ent']:
        fin = open(fname,'rt')
      else:
        msg = "   ERROR (PDBMapParser) Unsupported file type: %s.\n"%fname
        raise Exception(msg)
    except Exception as e:
      raise "In %s, could not open file %s"%(functioNameAsString,fname)
    try:
      p = PDBParser()
      filterwarnings('ignore',category=PDBConstructionWarning)
      s = p.get_structure(modelid,fin)
      resetwarnings()
      fin.close()
    except Exception as e:
      msg = "   ERROR (%s) Error while parsing %s from %s: %s"%(functionNameAsString,modelid,fname,str(e).replace('\n',' '))
    return s

  def get_structure(self,pdbid,fname,biounit_fnames=[],quality=-1,io=None):
    s = PDBMapParser.getBiopythonStructureOrFail(pdbid,fname)
    # import pdb; pdb.set_trace()
    try:
      s = PDBMapStructure(s,quality,pdb2pose={})
      # 2018-04-23 It cannot make sense to close a fle handle that is not local to us
      # fin.close()
    except Exception as e:
      msg = "ERROR (PDBMapIO) Error while parsing %s: %s"%(pdbid,str(e).replace('\n',' '))
      raise Exception(msg)

    # Clean up some common header differences
    if 'structure_method' in s.header and 'structure_methods' not in s.header:
      s.header['structure_methods'] = s.header['structure_method']
    if 'expdta' in s.header and 'structure_methods' not in s.header:
      s.header['structure_methods'] = s.header['expdta']
    if 'journal_reference' in s.header and 'journal' not in s.header:
      s.header['journal'] = s.header['journal_reference']
    if 'jrnl' in s.header and 'journal' not in s.header:
      s.header['journal'] = s.header['jrnl']
    if 'keywds' in s.header and 'keywords' not in s.header:
      s.header['keywords'] = s.header['keywds']
    if 'compnd' in s.header and 'compound' not in s.header:
      s.header['compound'] = s.header['compnd']

    # Default any missing fields
    if 'structure_methods' not in s.header:
      s.header['structure_methods'] = ''
    if 'journal' not in s.header:
      s.header['journal'] = ''
    if 'keywords' not in s.header:
      s.header['keywords'] = ''
    if 'compound' not in s.header:
      s.header['compound'] = ''
    if not s.header['resolution']:
      s.header['resolution'] = -1.0
 
    # Load the contents of the PDB file
    ext = os.path.basename(fname).split('.')[-1]
    if ext == 'gz':
      fin = [l for l in gzip.open(fname,'rt')]
    else:
      fin = [l for l in open(fname,'rt')]

    # Extract the SEQADV information and annotate any conflict residues
    for r in s.get_residues():
      if "conflict" not in dir(r):
        r.conflict = None # Initialize all residues to null conflict
    seqadv = [line for line in fin if line.startswith("SEQADV") and
                line[24:28].strip() == "UNP"]
    # Assign each conflict to the associated residue
    for row in seqadv:
      chain  = row[16]
      try:
        seqid  = int(row[18:22])
      except:
        # This SEQADV marks a deletion, which has no PDB position
        continue # Skip this entry and proceed to next
      cnflct = row[49:70].strip()
      if chain in s[0] and seqid in s[0][chain]:
        s[0][chain][seqid].conflict = cnflct

    # Extract the DBREF information or query from SIFTS
    dbref_fields = []
    if io:
      # Query structure-protein alignment information from SIFTS
      q  = "SELECT chain,uniprot_acc,min(resnum),max(resnum), "
      q += "min(uniprot_resnum),max(uniprot_resnum) FROM sifts "
      q += "WHERE pdbid=%s AND uniprot_acc!='' group by chain,uniprot_acc"
      res = io.secure_query(q,(pdbid,),cursorclass='Cursor')
      dbref_fields = [list(r) for r in res]
    # if len(dbref_fields) < 1:
      # # Attempt to parse DBREF fields if SIFTS unavailable
      # dbref_fields = [line for line in fin if line[0:6]=="DBREF " and
      #           line[26:32].strip() == "UNP"]
      # Attempt to read from the pdb file any information missing in SIFTS
      for line in fin:
         # If a chain is present in DBREF but not SIFTS, add it.
        if line[0:6].strip() == "DBREF" and line[26:32].strip() == "UNP":
          if line[12:14].strip() not in [r[0] for r in dbref_fields]:
            dbref_fields.append(line)

    if len(dbref_fields) < 1:
      msg = "   ERROR (PDBMapIO) No DBREF or SIFTS available for %s."%s.id
      raise Exception(msg)
    dbref = {}
    for ref in dbref_fields:
      if type(ref)==str:
        # PDB/ENT file
        chain    = ref[12:14].strip()
        unp      = ref[33:41].strip()
        species  = ref[42:54].strip().split('_')[-1]
        pdbstart = int(ref[14:18])
        pdbend   = int(ref[20:24])
        dbstart  = int(ref[55:60])
        dbend    = int(ref[62:67])
        hybrid   = 0
      else:
        # SIFTS query, species lookup in SwissProt
        chain,unp,pdbstart,pdbend,dbstart,dbend = ref
        if unp in PDBMapProtein.unp2species:
          species = PDBMapProtein.unp2species[unp]
        else:
          species = "UNKNOWN"
        hybrid  = 0
      # Documented offset between PDB and canonical sequence
      offset   = dbstart - pdbstart
      # Handle split chains
      drop_resis = []
      if 'unp' in dir(s[0][chain]):
        hybrid = 1 # Flag this chain as a hybrid
        # If existing chain segment is human, and this one is not
        if s[0][chain].species == 'HUMAN' and species != 'HUMAN':
          # Drop this segment of the chain
          drop = [r for r in s[0][chain] if pdbstart <= r.id[1] <= pdbend]
          for r in drop:
            # s[0][chain].detach_child(r.id)
            drop_resis.append(r.id)
          continue # Do not overwrite existing chain DBREF
        # Else if this chain segment is human, and the exiting is not   
        elif s[0][chain].species != 'HUMAN' and species == 'HUMAN':
          # Drop the existing segment of the chain
          nhstart,nhend = s[0][chain].pdbstart,s[0][chain].pdbend
          drop = [r for r in s[0][chain] if nhstart <= r.id[1] <= nhend]
          for r in drop:
            # s[0][chain].detach_child(r.id)
            drop_resis.append(r.id)
        # NO ACTION NECESSARY FOR THE REMAINING CONDITIONS.
        # If both chains are human, the first will be retained, the 
        # information will be the same, the start offset will default
        # to that specified in the first DBREF, which should precede
        # the new segment, and the new segment will be concatenated.
        # If both chains are non-human, then this chain will either be replaced
        # with a human chain-segment in a future DBREF, or it will be dropped
        # as a non-human chain during preprocessing.
      s[0][chain].unp  = unp
      gene = None if not species == "HUMAN" else PDBMapProtein.unp2hgnc(unp) 
      s[0][chain].gene = gene
      # Do not overwrite an earlier pdbstart (disordered hybrid chains)
      if pdbstart not in dir(s[0][chain]) or pdbstart < s[0][chain].pdbstart:
        s[0][chain].pdbstart = pdbstart
        s[0][chain].offset   = offset
      # Do not overwrite a later pdbend (disordered hybrid chains)
      if pdbend not in dir(s[0][chain]) or pdbend > s[0][chain].pdbend:
        s[0][chain].pdbend   = pdbend
      s[0][chain].species  = species
      s[0][chain].hybrid   = hybrid
      dbref[chain] = {'unp':unp,'gene':gene,'pdbstart':pdbstart,'offset':offset,'pdbend':pdbend,
                      'species':species,'hybrid':hybrid,'drop_resis':[r for r in drop_resis]}

    # Sanitize free text fields
    junk_stripper = str.maketrans('','',"'\"")
    s.header["name"]     = str(s.header["name"]    ).translate(junk_stripper)
    s.header["author"]   = str(s.header["author"]  ).translate(junk_stripper)
    s.header["keywords"] = str(s.header["keywords"]).translate(junk_stripper)
    s.header["compound"] = str(s.header["compound"]).translate(junk_stripper)
    s.header["journal"]  = str(s.header["journal"] ).translate(junk_stripper)
    s.header["structure_method"]    = str(s.header["structure_method"]).translate(junk_stripper)
    s.header["structure_reference"] = str(s.header["structure_reference"]).translate(junk_stripper)

    try:
      # Preprocess the asymmetric and biological units for this protein
      print("   # Processing the asymmetric unit")
      s = PDBMapParser.process_structure(s) # *must* occur before biounits
    except Exception as e:
      # If the structure is invalid (non-human, empty, errors, etc) pass along exception
      raise Exception("   ERROR (PDBMapIO): %s"%str(e))
    m = s[0]
    # Calculate relative solvent accessibility for this model
    ssrsa = {}
    resinfo = PDBMapIO.dssp(m,fname,dssp_executable='dssp')
    for c in m:
      for r in c:
        if (c.id,r.seqid,r.icode) not in resinfo:
          info = dict((i,None) for i in ['ss','rsa','phi','psi','tco','kappa','alpha'])
        else:
          info  = resinfo[(c.id,r.seqid,r.icode)]
        r.ss  = info['ss']
        r.rsa = info['rsa']
        r.phi = info['phi']
        r.psi = info['psi']
        r.tco = info['tco']
        r.k   = info['kappa']
        r.a   = info['alpha']
        ssrsa[(c.id,r.seqid,r.icode)] = (r.ss,r.rsa,r.phi,r.psi,r.tco,r.k,r.a)

    # Process the biological assemblies for this structure
    for biounit_fname in biounit_fnames:
      try:
        biounit = PDBMapParser.getBiopythonStructureOrFail(pdbid,fname)
        if os.path.basename(biounit_fname).split('.')[-1] == 'gz':
          bioid = int(os.path.basename(biounit_fname).split('.')[-2][3:])
        else:
          bioid = int(os.path.basename(biounit_fname).split('.')[-1][3:])
        biounit = PDBMapStructure(biounit,pdb2pose={}) # must pass empty dictionary: python bug
      except Exception as e:
        msg = "   ERROR (PDBMapIO) Error while parsing %s biounit %d: %s"%(pdbid,bioid,str(e).replace('\n',' '))
        LOGGER.exception(msg)
        raise Exception(msg)
      print("   # Processing biounit %d"%bioid)
      biounit = PDBMapParser.process_structure(biounit,biounit=bioid,dbref=dbref)
      if not biounit:
        msg = "   ERROR (PDBMapIO) Biological assembly %s.%d contains no human protein chains.\n"%(pdbid,bioid)
        LOGGER.exception(msg)
        continue
      # Add the models for this biological assembly to the PDBMapStructure
      for m in biounit:
        m.id = "%d.%d"%(m.biounit,m.id)
        for c in m:
          for r in c:
            r.ss,r.rsa,r.phi,r.psi,r.tco,r.k,r.a = ssrsa[(c.id,r.seqid,r.icode)]
            if (c.id,r.id[1]) in seqadv:
              r.conflict = seqadv[(c.id,r.id[1])]
            else:
              r.conflict = None
        s.add(m)
    return s

  @classmethod
  def process_structure_dssp_unp2hgnc(cls, m, model_summary,fname,unp=None):
    unp = unp if unp else model_summary['unp']
    s = PDBMapParser.process_structure(m)
    m = s[0]
    s.unp = unp
    m.unp = unp
    resinfo = PDBMapIO.dssp(m,fname,dssp_executable='dssp')
    n = 0
    for c in m:
      c.unp  = unp
      c.gene = PDBMapProtein.unp2hgnc(unp)
      for r in c:
        info  = resinfo[(c.id,r.seqid,r.icode)]
        r.ss  = info['ss']
        r.rsa = info['rsa']
        r.phi = info['phi']
        r.psi = info['psi']
        r.tco = info['tco']
        r.k   = info['kappa']
        r.a   = info['alpha']
        r.conflict = None
    return s
