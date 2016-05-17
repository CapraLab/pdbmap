#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : PDBMapStructure.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-12
# Description    : Wrapper for Bio.PDB.Structure.Structure with additional
#                : information and functionality pertinent to PDBMap. All
#                : requests for Structure attributes are deferred to the
#                : Structure object.
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv,copy,time,random,tempfile
import subprocess as sp
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB import Superimposer
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
import lib.PDBMapIO
from lib.PDBMapTranscript import PDBMapTranscript
from lib.PDBMapAlignment import PDBMapAlignment
from warnings import filterwarnings,resetwarnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Import helper for PyRosetta
def import_rosetta():
  if not import_rosetta.imported:
    try:
      global rosetta,mutants
      import rosetta
      import lib.mutants as mutants
      rosetta.init()
      import_rosetta.imported = True
    except:
      sys.stderr.write("PyRosetta not found. Rosetta utilities are unavailable.\n")
import_rosetta.imported = False

from multiprocessing import Pool,cpu_count
def unwrap_self_mutate(arg,**kwargs):
  return PDBMapStructure.mutate(*arg,**kwargs)
def unwrap_self_relax(arg,**kwargs):
  return PDBMapStructure.relax(*arg,**kwargs)

class PDBMapStructure(Structure):

  def __init__(self,s,quality=-1,pdb2pose={},refseq=None,alignment={}):
    # Assign the Structure, and quality
    if isinstance(s,PDBMapStructure):
      self = copy.deepcopy(s)
    else:
      self.structure   = s
      self.id          = s.id
      self.quality     = quality
      self.transcripts = []
      self.alignments  = []
      self._pdb2pose   = pdb2pose
      if not self._pdb2pose:
        for m in s:
          self._pdb2pose[m.id] = {}
          for c in m:
            self._pdb2pose[m.id][c.id] = {}
        for i,r in enumerate(s.get_residues()):
          c = r.get_parent()
          mid,cid = c.get_parent().id,c.id
          self._pdb2pose[mid][cid][r.id[1]] = i+1
      else:
        # Revert indexing to original PDB indexing per pdb2pose
        for mid,d1 in self._pdb2pose.iteritems():
          for cid,d2 in d1.iteritems():
            # The pose ID should always be less than the PDB ID
            for rid,pose_id in sorted(d2.iteritems(),reverse=True):
              # Isolate the residue
              res = self.structure[mid][cid][pose_id]
              # Detach from the chain
              self.structure[mid][cid].detach_child(res.id)
              # Update residue's position in the chain
              resid = list(res.id)
              resid[1] = rid
              res.id = tuple(resid)
              # Add residue back to the chain
              self.structure[mid][cid].add(res)
            # Reverse the residue list to correct for backward algorithm
            cl = self.structure[mid][cid].child_list
            self.structure[mid][cid].child_list = cl[::-1]
      # Align to reference sequence if one is provided
      self.refseq = refseq
      if self.refseq:
        print "Aligning %s to reference sequence"%self.id
        self.align2refseq(self.id,refseq)

  def __str__(self):
    return self.structure.id

  def __getstate__(self):
    return self.structure,self.quality,self.transcripts,self._pdb2pose

  def __setstate__(self,state):
    self.structure,self.quality,self.transcripts,self._pdb2pose = state

  def __getattr__(self,attr):
    # Defer appropriate calls to the internal structure
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

  def get_residue(self,chain,seqid,model=0,refpos=False,strict=False):
    if not refpos:
      try:
        return self.structure[model][chain][seqid]
      except:
        return None
    else:
      # Adjust for alignment between reference and structure
      # print "Reference position %d is aligned with PDB position..."%seqid,
      if seqid in self.structure[model][chain].alignment.seq2pdb:
        seqid = self.structure[model][chain].alignment.seq2pdb[seqid]
        # print seqid
      else:
        if strict: raise Exception("PDB cannot map position %d"%seqid)
        # else: print 'NA'; return None
        else: return None
      #FIXME: This is only necessary if the underlying structure *is* a pose
      #     : There's nothing to say that this is an appropriate expectation
      # # Adjust for alignment between structure and pose
      # print "Structure position #%d:%d.%s is aligned with pose position..."%(model,seqid,chain),
      # if seqid in self._pdb2pose[model][chain]:
      #   seqid = self._pdb2pose[model][chain][seqid]
      #   print seqid
      # else:
      #   if strict: raise Exception("Pose cannot map position %d"%seqid)
      #   # else: print 'NA'; return None
      #   else: return None
      # print ""
      try:
        return self.structure[model][chain][seqid]
      except:
        return None

  def align2refseq(self,sid,refseq):
    for c in self.get_chains():
      refdict = dict((i+1,(r,"NA",0,0,0)) for i,r in enumerate(refseq))
      c.transcript = PDBMapTranscript("ref","ref","ref",refdict)
      c.alignment  = PDBMapAlignment(c,c.transcript)
      self.transcripts.append(c.transcript)
      self.alignments.append(c.alignment)

  def pose(self):
    """ Loads the PDBMapStructure as a Rosetta::Pose object """
    import_rosetta()
    io = PDBIO()
    io.set_structure(self.structure)
    with tempfile.NamedTemporaryFile('wrb',suffix='.pdb',delete=False) as tf:
      io.save(tf.name)
    pose = rosetta.Pose()
    rosetta.pose_from_pdb(pose,tf.name)
    os.remove(tf.name)
    return pose

  def from_pose(self,pose):
    """ Updates the PDBMapStructure from a Rosetta::Pose object """
    with tempfile.NamedTemporaryFile('wrb',suffix='.pdb',delete=False) as tf:
      pose.dump_pdb(tf.name)
    p = PDBParser()
    with open(tf.name,'rb') as fin:
      filterwarnings('ignore',category=PDBConstructionWarning)
      s = p.get_structure(self.id,tf.name)
      p = lib.PDBMapIO.PDBMapParser()
      s = p.process_structure(s,force=True)
      resetwarnings()
    os.remove(tf.name)
    return PDBMapStructure(s,pdb2pose={},refseq=self.refseq)

  def get_transcripts(self,io=None):
    # Retrieve the corresponding transcript for each chain
    # Check if transcripts have been previously identified
    error_msg = ""
    if self.transcripts:
      return self.transcripts
    # Identify and align corresponding transcripts
    # prot2chain = {}
    for chain in self.structure[0]:
      print "Getting transcripts for %s.%s"%(self.id,chain.id)
      # # If a chain of the same protein has already been solved, use solution
      # # BUT NOTE: must update Alignment->Chain backreference after deep copy
      # if chain.unp in prot2chain:
      #   chain.alignments = [copy.deepcopy(a) for a in prot2chain[chain.unp]]
      #   # Update the back-reference on the copied Alignment to point to this chain
      #   for a in chain.alignments:
      #     a.chain = chain
      # else:
      # Query all transcripts associated with the chain's UNP ID
      candidate_transcripts = PDBMapTranscript.query_from_unp(chain.unp)
      if len(candidate_transcripts) < 1:
        error_msg += "UniProt indicates no EnsEMBL transcripts for %s.%s (%s)"%(self.id,chain.id,chain.unp)
        # raise Exception("ERROR (PDBMapStructure): %s\n"%error_msg)
      # Align chains candidate transcripts
      alignments = {}
      for trans in candidate_transcripts:
        alignment = PDBMapAlignment(chain,trans,io=io)
        # Exclude alignments with <90% identity, likely bad matches
        if alignment.perc_identity >= 0.9:
          # Setup tuples for primary sort on length, secondary sort on transcript name (for tie-breaking consistency)
          if alignment.transcript.gene not in alignments:
            alignments[alignment.transcript.gene] = [(len(alignment.transcript.sequence),alignment.transcript.transcript,alignment)]
          else:
            alignments[alignment.transcript.gene].append((len(alignment.transcript.sequence),alignment.transcript.transcript,alignment))
        else:
          # Note that at least one transcript was dropped due to low alignment quality
          error_msg += "%s (%s.%s (%s)) dropped due to low alignment quality (%.2f); "%(trans.transcript,self.id,chain.id,chain.unp,alignment.perc_identity)
      # Store canonical transcript for each gene alignment as element of chain
      chain.alignments = []
      # prot2chain[chain.unp] = []
      for gene in alignments:
        alignments[gene].sort() # ascending by transcript length, then name
        if len(alignments[gene]) > 0:
          chain.alignments.append(alignments[gene][-1][-1]) # last alignment (longest) length
      # prot2chain[chain.unp] = chain.alignments
      # Recover transcripts from alignments
      chain.transcripts = [a.transcript for a in chain.alignments]
    # Return the matched transcripts
    try:
      self.transcripts = [t for c in self.structure[0] for t in c.transcripts]
      self.alignments  = [a for c in self.structure[0] for a in c.alignments]
    except:
      self.transcripts = []
      self.alignments  = []
    if not self.transcripts:
      raise Exception("ERROR (PDBMapStructure): %s"%error_msg)
    return self.transcripts

  def get_alignments(self):
    if not self.alignments:
      self.get_transcripts()
    else:
      return self.alignments  

  def permute(self,nperm=1000):
    """ Generator: permutes SNP assignments within the structure """
    # Record the original state of SNPMapper
    ts = copy.deepcopy(self.structure.snpmapper)
    # Resolve the SNPMapper to an alias
    sm = self.structure.snpmapper
    unps = sm.keys()
    # Chains may not cover full length of protein sequence
    # Only permute protein sequence present in the structure
    pr = {} # protein range
    for m in self.structure:
      for c in m:
        seqs = [r.seqid for r in c]
        minseq = np.min(seqs)
        maxseq = np.max(seqs)
        unp = c.unp
        if unp not in pr:
          pr[unp] = [minseq,maxseq]
        elif pr[unp][0]>minseq:
          pr[unp][0] = minseq
        elif pr[unp][1]<maxseq:
          pr[c.unp][1] = maxseq
    try:
      for i in range(nperm):
        for unp in unps:
          # Shuffle between min/max observed residues
          np.random.shuffle(sm[unp][pr[unp][0]:pr[unp][1]])
        # Yield this structure object with each permutation
        yield self
    except: 
      raise Exception('Structure must be initialized with SNP assignments.')
    finally:
      # Return the SNPMapper to its original state
      self.structure.snpmapper = ts
      for m in self.structure:
        for c in m:
          for r in c:
            r.snp = self.snpmapper[c.unp][r.seqid-1]

  def contains(m,c,r):
    return m in structure and c in s.structure[m] and r in s.structure[m][c]

  def snps(self):
    for m in self.structure:
      for c in m:
        for r in c:
          if r.snp[0]:
            yield (c.id,r.seqid,r.snp[1])

  def clean(self):
    p  = PDBParser()
    io = PDBIO()
    with tempfile.NamedTemporaryFile('wrb',suffix='.pdb',delete=False) as tf:
      io.set_structure(self.structure)
      io.save(tf.name)
    cmd  = ['lib/clean_pdb.py',tf.name,'ignorechain','nopdbout']
    print "\n%s"%' '.join(cmd)
    proc = sp.Popen(cmd,stdout=sp.PIPE)
    s = p.get_structure(self.get_id(),proc.stdout)
    p = lib.PDBMapIO.PDBMapParser()
    s = p.process_structure(s,force=True)
    s = PDBMapStructure(s,pdb2pose=self._pdb2pose,refseq=self.refseq)
    os.remove(tf.name)
    return s

  def mutate(self,mut,strict=True):
    """ Point mutation """
    pose = self.pose()
    c,m = mut
    c = c if c else ' '
    a1,r,a2 = m[0],int(m[1:-1]),m[-1]
    print "\nSimulating mutation: %s%s%s"%(a1,r,a2)
    # Adjust for alignment between reference and structure
    if "alignment" in dir(self.structure[0][c]):
      print "Reference position %d is aligned with PDB position..."%r,
      if r in self.structure[0][c].alignment.seq2pdb:
        r = self.structure[0][c].alignment.seq2pdb[r]
      else:
        if strict: raise Exception("PDB cannot map position %d"%r)
        else: print 'NA'; return None
      print r
    # Adjust for alignment between structure and pose
    print "The reference allele is %s"%a1
    print "The observed  allele is %s"%self.structure[0][c][r].rescode
    print "The alternate allele is %s\n"%a2
    # Check that the reference amino acid matches observed amino acid
    if not self.structure[0][c][r].rescode == a1:
      if strict: raise Exception("Reference allele does not match.")
      else: return None
    print "Structure position %d is aligned with pose position..."%r,
    if r in self._pdb2pose[0][c]:
      r = self._pdb2pose[0][c][r]
    else:
      if strict: raise Exception("Pose cannot map position %d"%r)
      else: print 'NA'; return None
    print r
    try:
      print "Inserting the alternate allele...\n"
      pose = mutants.mutate_residue(pose,r,a2,pack_radius=10)
    except:
      if strict: raise
      else: return None
    return self.from_pose(pose)

  def pmutate(self,muts,maxprocs=cpu_count(),strict=True):
    """ Handles parallelization of point mutations """
    if len(muts) < 2:
      print "\n## Only one mutation to model. No workers spawned. ##\n"
      return [self.mutate(muts[0],strict=strict)]
    else:
      print "\n## Modeling mutations in serial. No workers spawned. ##\n"
      return [self.mutate(m,strict=strict) for m in muts]
    # print "Spawning a pool of %d mutate workers"%(min(len(muts),maxprocs))
    # pool = Pool(processes=min(len(muts),maxprocs))
    # return pool.map(unwrap_self_mutate,zip([self]*len(muts),muts,[strict]*len(muts)))

  def prelax(self,iters=1,maxprocs=cpu_count()):
    """ Handles parallelization of relaxation """
    if iters < 2:
      print "\n## Only one iteration of FastRelax. No workers spawned. ##\n"
      sr    = self.relax()
      sc    = sr.score()
      rmsd  = sr.rmsd(self)
      farep = sr.fa_rep()
      return sr,sc,np.nan,np.nan,rmsd,np.nan,np.nan,farep
    print "\n## Spawning a pool of %d relax workers ##\n"%(min(iters,maxprocs))
    pool     = Pool(processes=min(iters,maxprocs))
    ensemble = pool.map(unwrap_self_relax,zip([self]*iters))
    scores = [e.score()    for e in ensemble]
    rmsds  = [e.rmsd(self) for e in ensemble]
    return ensemble[np.argmin(scores)],np.min(scores),np.mean(scores),np.std(scores),np.min(rmsds),np.mean(rmsds),np.std(rmsds),ensemble[np.argmin(scores)].fa_rep()

  def relax(self):
    """ Apply Rosetta:FastRelax """
    import_rosetta()
    pose  = self.pose()
    relax = rosetta.FastRelax()
    relax.constrain_relax_to_start_coords(True)
    relax.set_scorefxn(rosetta.create_score_function('talaris2014'))
    # relax.set_scorefxn(rosetta.create_score_function_ws_patch("standard", "score12"))
    relax.apply(pose)
    return self.from_pose(pose)

  def score(self,norm=True):
    """ Use Rosetta::Score to evaluate a structure """
    import_rosetta()
    sfxn = rosetta.create_score_function('talaris2014')
    # sfxn = rosetta.create_score_function_ws_patch("standard", "score12")
    p    = self.pose()
    sc   = sfxn(p)
    if norm:
      return sc / p.n_residue()
    return sc

  def fa_rep(self,norm=True):
    import_rosetta()
    sfxn = rosetta.create_score_function('talaris2014')
    # sfxn = rosetta.create_score_function_ws_patch("standard", "score12")
    p    = self.pose()
    sfxn(p)
    fa_rep = p.energies().total_energies()[rosetta.fa_rep]
    if norm:
      return fa_rep / p.n_residue()
    return fa_rep

  def rmsd(self,s2,norm100=True):
    import_rosetta()
    p = self.pose()
    rmsd = rosetta.CA_rmsd(p,s2.pose())
    if not norm100:
      return rmsd
    return rmsd * 100. / p.n_residue()

  # DEPRECATED #
  def scwrl(self,altseq):
    """ Repacks sidechains using SCWRL4 and returns a copy """
    io = PDBIO()
    seqfname = "temp/%d.txt"%multidigit_rand(10)
    with open(seqfname,'wb') as seqfile:
      structfile = "temp/%d.pdb"%multidigit_rand(10)
      seqfile.write(altseq)
      scwrlfile = structfile+".scwrl"
      io.set_structure(self.structure)
      io.save(structfile)
    cmd = ["scwrl","-0","-i",structfile,'-s',seqfname,'-o',scwrlfile]
    print "\n%s"%' '.join(cmd)
    sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE).communicate()
    p = PDBParser()
    with open(scwrlfile,'rb') as fin:
      filterwarnings('ignore',category=PDBConstructionWarning)
      s = p.get_structure(self.id,scwrlfile)
      resetwarnings()
    s = PDBMapStructure(s,pdb2pose={},refseq=self.refseq)
    os.remove(structfile)
    os.remove(scwrlfile)
    os.remove(seqfname)
    return s

## Copied from biolearn
def multidigit_rand(digits):
  randlist = [random.randint(1,10) for i in xrange(digits)]
  multidigit_rand = int(''.join([str(x) for x in randlist]))
  return multidigit_rand

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)

# This class was a pain in the ass to write. Thank you to:
# http://stackoverflow.com/questions/1466676/create-a-wrapper-class-
# to-call-a-pre-and-post-function-around-existing-functions
