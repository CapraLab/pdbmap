#!/usr/bin/env python
#
# Project        : PDBMap
# Filename       : PDBMapIO.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-09
# Description    : Manages all database and IO interactions.
#=============================================================================#


# See main check for cmd line parsing
import sys,os,csv,collections,gzip,time,random
import subprocess as sp
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
import Bio.PDB
import tempfile
import numpy as np
from .PDBMapModel import PDBMapModel
from .PDBMapSwiss import PDBMapSwiss
from .PDBMapProtein import PDBMapProtein
from .PDBMapStructure import PDBMapStructure
import MySQLdb, MySQLdb.cursors
# Support for caching sql queries
import hashlib
import pickle as pickle

from warnings import filterwarnings,resetwarnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import logging
LOGGER = logging.getLogger(__name__)

class PDBMapIO(PDBIO):
  def __init__(self,dbhost=None,dbuser=None,dbpass=None,dbname=None,slabel="",dlabel="",createdb=False):
    super(PDBMapIO,self).__init__()
    self.dbhost = dbhost
    self.dbuser = dbuser
    self.dbpass = dbpass
    self.dbname = dbname
    self.slabel = slabel
    self.dlabel = dlabel
    self._cons  = []   # Define the connection pool
    self._con   = None # Define the connection
    self._c     = None # Define the cursor
    self.check_schema()
    # Open the connection only if this wasn't a database creation
    if not createdb:
      self._connect()

  def __del__(self):
    # Guarantee all database connections are closed
    for con in self._cons:
      # LOGGER.info("Calling con.close() to close SQL connection")
      con.close()


  def is_nmr(self,pdbid,label=-1):
    # None is a valid argument to label
    if label == -1:
      label=self.slabel
    self._connect()
    query  = "SELECT * FROM Structure WHERE pdbid=%s and method like '%%nmr%%' "
    if label:
      query += "AND label=%s "
    query += "LIMIT 1"
    if label:
      self._c.execute(query,(pdbid,label))
    else:
      self._c.execute(query,(pdbid,))
    res = True if self._c.fetchone() else False
    self._close()
    return res

  def structure_in_db(self,pdbid,label=-1):
    # None is a valid argument to label
    if label == -1:
      label=self.slabel
    self._connect()
    query  = "SELECT * FROM Structure WHERE pdbid=%s "
    if label:
      query += "AND label=%s "
    query += "LIMIT 1"
    if label:
      self._c.execute(query,(pdbid,label))
    else:
      self._c.execute(query,(pdbid,))
    res = True if self._c.fetchone() else False
    self._close()
    return res

  def model_in_db(self,modelid,label=-1):
    # None is a valid argument to label
    if label == -1:
      label=self.slabel
    self._connect()
    query = "SELECT * FROM Model WHERE modelid=%s "
    if label:
      query += "AND label=%s "
    query += "LIMIT 1"
    if label:
      self._c.execute(query,(modelid,label))
    else:
      self._c.execute(query,(modelid,))
    res = True if self._c.fetchone() else False
    self._close()
    return res

  def swiss_in_db(self,modelid,label=-1):
    # None is a valid argument to label
    if label == -1:
      label=self.slabel
    self._connect()
    query = "SELECT * FROM Swiss WHERE modelid=%s "
    if label:
      query += "AND label=%s "
    query += "LIMIT 1"
    if label:
      self._c.execute(query,(modelid,label))
    else:
      self._c.execute(query,(modelid,))
    res = True if self._c.fetchone() else False
    self._close()
    return res

  def unp_in_db(self,unpid,label=-1):
    # None is a valid argument to label
    if label == -1:
      label=self.slabel
    self._connect()
    query = "SELECT * FROM Chain WHERE (unp=%s OR unp LIKE %s) "
    if label:
      query += "AND label=%s"
    query += "LIMIT 1"

    base_unp = unpid.split('-')[0]
    if label:
      self._c.execute(query,(base_unp,base_unp+"-%%",label))
    else:
      self._c.execute(query,(base_unp,base_unp+"-%%"))
    res = True if self._c.fetchone() else False
    self._close()
    return res

  def gene_in_db(self,gene,label=-1):
    # None is a valid argument to label
    if label == -1:
      label=self.slabel
    self._connect()
    query  = "SELECT unp FROM Chain a WHERE gene=%s "
    if label:
      query += "AND label=%s"
    query += "LIMIT 1"
    if label:
      self._c.execute(query,(gene,label))
    else:
      self._c.execute(query,(gene,))
    res = self._c.fetchone()
    self._close()
    if res:
      return True,res['unp'][0]
    return False,None

  def genomic_datum_in_db(self,name,label=None):
    self._connect()
    query  = "SELECT * FROM Structure WHERE name=%s "
    if label:
      query += "AND label=%s "
    query += "LIMIT 1"
    if label:
      self._c.execute(query,(name,label))
    else:
      self._c.execute(query,name)
    res = True if self._c.fetchone() else False
    self._close()
    return res

  def upload_structure(self,model=False,sfields=None):
    """ Uploads the current structure in PDBMapIO """
    # Verify that structure is not already in database
    if self.structure_in_db(self.structure.id):
      msg =  "WARNING (PDBMapIO) Structure %s "%self.structure.id
      msg += "already in database.\n"
      sys.stderr.write(msg)
      return(1)
    # Upload entire structure (structure, chains, residues)
    queries = []
    s = self.structure

    # Upload the structure if skip not specified
    if not model:
      h = s.header
      squery  = 'INSERT IGNORE INTO Structure '
      squery += '(label,pdbid,method,quality,resolution,`name`,author,deposition,`release`,compound,keywords,reference,structure_reference) '
      squery += 'VALUES ('
      squery += '"%(label)s","%(id)s","%(structure_method)s",%(quality)f,'
      squery += '%(resolution)f,"%(name)s","%(author)s","%(deposition_date)s",'
      squery += '"%(release_date)s","%(compound)s","%(keywords)s",'
      squery += '"%(journal)s","%(structure_reference)s")'
      sfields = dict((key,s.__getattribute__(key)) for key in dir(s) 
                    if isinstance(key,collections.Hashable))
      # Add the header fields
      sfields.update(s.header)
      # Add the internal structure fields
      sfields.update(dict((key,s.structure.__getattribute__(key)) for key in 
                      dir(s.structure) if isinstance(key,collections.Hashable)))
      sfields["label"] = self.slabel
      squery = squery%sfields
      queries.append(squery)

    # Upload the models, chains,and residues
    for m in s:
      mfields = dict((key,m.__getattribute__(key)) for key in dir(m) 
                      if isinstance(key,collections.Hashable))
      if mfields['biounit'] > 0: # remove biounit tag
        mfields['id'] = int(mfields['id'].split('.')[-1])
        mfields['id'] += 1 # correct indexing
      for c in m:
        cquery  = "INSERT IGNORE INTO Chain "
        cquery += "(label,structid,biounit,model,chain,unp,gene,offset,hybrid,sequence) "
        cquery += "VALUES "
        cquery += '("%(label)s","%(id)s",'%sfields # structure id
        cquery += '%(biounit)d,%(id)d,'%mfields # model id
        cquery += '"%(id)s","%(unp)s","%(gene)s",%(offset)d,%(hybrid)d,"%(sequence)s")'
        cfields = dict((key,c.__getattribute__(key)) for key in dir(c) 
                        if isinstance(key,collections.Hashable))
        cfields["label"] = self.slabel
        cquery = cquery%cfields
        queries.append(cquery)
        rquery  = "INSERT IGNORE INTO Residue "
        rquery += "(label,structid,biounit,model,chain,resname,rescode,seqid,icode,x,y,z,ss,rsa,phi,psi,tco,kappa,alpha,conflict) "
        rquery += "VALUES "
        for r in c:
          # If length exceeds 10 thousand characters, start a new query
          if len(rquery) > 10000:
            queries.append(rquery[:-1])
            rquery  = "INSERT IGNORE INTO Residue "
            rquery += "(label,structid,biounit,model,chain,resname,rescode,seqid,icode,x,y,z,ss,rsa,phi,psi,tco,kappa,alpha,conflict) "
            rquery += "VALUES "
          # Continue with the next query
          rquery += '("%(label)s","%(id)s",'%sfields # structure id
          rquery += '%(biounit)d,%(id)d,'%mfields # model id
          rquery += '"%(id)s",'%cfields  # chain id
          rquery += '"%(resname)s","%(rescode)s",%(seqid)d,'
          rquery += '"%(icode)s",%(x)f,%(y)f,%(z)f,'
          rquery += '"%(ss)s",%(rsa)s,%(phi)s,%(psi)s,'
          rquery += '%(tco)s,%(k)s,%(a)s,%(conflict)s),'
          rfields = dict((key,r.__getattribute__(key)) for key in dir(r) 
                          if isinstance(key,collections.Hashable))
          rfields["conflict"] = '"%s"'%rfields["conflict"] if rfields["conflict"] else '\\N'
          for key,val in list(rfields.items()):
            if val is None:
              if key == 'ss':
                rfields[key] = '?'
              else:
                rfields[key] = 'NULL'
          rfields["label"] = self.slabel
          rquery = rquery%rfields
        queries.append(rquery[:-1])

    # Upload the transcripts
    # Inner exception if empty list, outer exception if error during query
    if not len(s.get_transcripts(io=self)):
      msg = "No transcripts for structure %s, proteins: %s"%(s.id,','.join([c.unp for m in s for c in m]))
      LOGGER.critical(msg)
      return (0)

    tquery  = "INSERT IGNORE INTO Transcript "
    tquery += "(label,transcript,protein,gene,seqid,rescode,chr,start,end,strand) "
    tquery += "VALUES "
    seen = set([])
    for t in s.get_transcripts(io=self):
      if t.transcript not in seen or seen.add(t.transcript):
        for seqid,(rescode,chr,start,end,strand) in list(t.sequence.items()):
          # If length exceeds 10 thousand characters, start a new query
          if len(tquery) > 10000:
            queries.append(tquery[:-1])
            tquery  = "INSERT IGNORE INTO Transcript "
            tquery += "(label,transcript,protein,gene,seqid,rescode,chr,start,end,strand) "
            tquery += "VALUES "
          tquery += '("%s","%s","%s","%s",'%(self.slabel,t.transcript,t.protein,t.gene)
          tquery += '%d,"%s",'%(seqid,rescode)
          tquery += '"%s",%d,%d,%d),'%(chr,start,end,strand)
    queries.append(tquery[:-1])
    # except Exception as e:
    #  msg = "ERROR (PDBMapIO) Failed to INSERT transcripts for %s: %s"%(s.id,str(e).rstrip('\n'))
    #  raise Exception(msg)

    # Execute all queries at once to ensure everything completed.
    try:
      self._connect()
      for q in queries[::-1]:
        self._c.execute(q)
      self._con.commit()
    except:
      msg  = "ERROR (PDBMapIO) Query failed for %s: "%s.id
      msg += "%s\n"%self._c._last_executed
      LOGGER.exception(msg)
      self._con.rollback()
      raise
    finally:
      self._close()

    # Upload the alignments
    # These should be done alignment by alignment as there is a master:detail relationship here to maintain
    for a in s.get_alignments():
      aquery  = "INSERT IGNORE INTO AlignChainTranscript "
      aquery += "(label,structid,chain,transcript,score,perc_aligned,perc_identity,alignment) " # chain_seqid,transcript,trans_seqid) "
      aquery += 'VALUES ("%s","%s","%s","%s",%f,%f,%f,"%s")'%(self.slabel,s.id,a.chain.id,a.transcript.transcript,a.score,a.perc_aligned,a.perc_identity,a.aln_str)

      last_al_id = 0
      # Execute all queries at once to ensure everything completed.
      try:
        self._connect()
        self._c.execute(aquery)
        last_al_id = self._c.lastrowid
        # self._con.commit() # Don't commit quite yet!
      except:
        msg  = "ERROR (PDBMapIO) Query failed for %s: "%s.id
        msg += "%s\n"%self._c._last_executed
        LOGGER.exception(msg)
        self._con.rollback()
        raise

      if last_al_id == 0: # This is a drag - but if the INSERT was IGNORED we have to get the last_al_id the hard way
        self._c.execute("SELECT al_id FROM AlignChainTranscript where label = %s and structid = %s and chain = %s and transcript = %s",
                                 (self.slabel,s.id,a.chain.id,a.transcript.transcript))
        for (temp_al_id) in cursor:
          last_al_id = temp_al_id
          break;

      assert (last_al_id != 0),"al_id referential problems in AlignChainTranscript"

      # We'll close after writing all the rewidue mappings below
      # finally:
      #  self._close()

      queries = []
      INSERT_VALUES = "INSERT IGNORE INTO AlignChainTranscriptResidue (chain_res_num, chain_res_icode, trans_seqid, al_id) VALUES "
      aquery  =       INSERT_VALUES

      for c_seqid,t_seqid in list(a.pdb2seq.items()):
        # If length exceeds 10 thousand characters, start a new query
        if len(aquery) > 10000:
          queries.append(aquery[:-1])
          aquery  = INSERT_VALUES
        aquery += '(%d,"%s",%d,%d),'%(c_seqid[1],c_seqid[2],t_seqid,last_al_id)

        queries.append(aquery[:-1])

      # Execute all queries at once to ensure everything completed.
      try:
        for q in queries[::-1]:
          self._c.execute(q)
        self._con.commit()
      except:
        msg  = "ERROR (PDBMapIO) Query failed for %s: "%s.id
        msg += "%s\n"%self._c._last_executed
        sys.stderr.write(msg)
        self._con.rollback()
        raise
      finally:
        self._close()
    # End for a in get_alignments()
    return(0)

  def delete_all_modbase (self):
    self._connect()
    deleteCommands = [
      "DELETE FROM GenomicIntersection WHERE slabel='modbase'",
      "DELETE FROM Alignment WHERE label='modbase'",
      "DELETE FROM AlignmentScore WHERE label='modbase'",
      "DELETE FROM Transcript WHERE label='modbase'",
      "DELETE FROM Residue WHERE label='modbase'",
      "DELETE FROM Chain WHERE label='modbase'",
      "DELETE FROM Model"]

    for mquery in deleteCommands:
      LOGGER.info(mquery)
      try:
        self._c.execute(mquery)
      except (MySQLdb.Error, MySQLdb.Warning) as e:
        LOGGER.critical(e)
        sys.exit(1)
    self._close()

  def upload_model(self):
    """ Uploades the current model in PDBMapIO """
    m = self.structure
    if self.model_in_db(m.id,self.slabel):
      msg =  "WARNING (PDBMapIO) Structure %s "%m.id
      msg += "already in database.\n"
      sys.stderr.write(msg)
      return(1)

    # Prepare the Model summary information (if no errors occurred)
    mquery  = 'INSERT IGNORE INTO Model '
    mquery += '(label,modelid,unp,method,no35,rmsd,mpqs,evalue,ga341,zdope,pdbid,chain,identity)'
    mquery += 'VALUES ('
    mquery += '"%(label)s","%(id)s","%(unp)s","%(tsvmod_method)s",'
    mquery += '%(tsvmod_no35)s,%(tsvmod_rmsd)s,%(mpqs)s,%(evalue)s,%(ga341)f,'
    mquery += '%(zdope)s,"%(pdbid)s","%(chain)s",%(identity)f)'
    mfields = dict((key,m.__getattribute__(key)) for key in dir(m) 
                  if isinstance(key,collections.Hashable))
    mfields["label"] = self.slabel
    mquery = mquery%mfields
    # First, pass the underlying PDBMapStructure to upload_structure
    self.upload_structure(model=True,sfields=mfields)
    # Then execute the Model upload query
    self._connect()
    self._c.execute(mquery)
    self._close()

  def delete_all_swiss(self):
    self._connect()
    deleteCommands = [
      "DELETE FROM GenomicIntersection WHERE slabel='swiss'",
      "DELETE FROM Alignment WHERE label='swiss'",
      "DELETE FROM AlignmentScore WHERE label='swiss'",
      "DELETE FROM Transcript WHERE label='swiss'",
      "DELETE FROM Residue WHERE label='swiss'",
      "DELETE FROM Chain WHERE label='swiss'",
      "DELETE FROM Swiss"]

    for mquery in deleteCommands:
      LOGGER.info(mquery)
      try:
        self._c.execute(mquery)
      except (MySQLdb.Error, MySQLdb.Warning) as e:
        LOGGER.critical(e)
        sys.exit(1)
    self._close()

  def upload_swiss(self,remark3_metrics):
    """ Uploades the current swissmodel in PDBMapIO """
    m = self.structure
    if self.swiss_in_db(m.id,self.slabel):
      msg =  "WARNING (PDBMapIO) Structure %s "%m.id
      msg += "already in swiss database.\n"
      sys.stderr.write(msg)
      return(1)

  # Prepare the Model summary information (if no errors occurred)
    mquery  = 'INSERT IGNORE INTO Swiss '
    mquery += '(label,modelid,unp,start,end,qmean,qmean_norm,template,coordinate_id,url,pdbid,chain,identity,ostat,mthd)'
    mquery += 'VALUES ('
    mquery += '"%(label)s","%(id)s","%(unp)s",'
    mquery += '%(start)d,%(end)d,'
    mquery += '%(qmean)f,%(qmean_norm)f,"%(template)s",'
    mquery += '"%(coordinate_id)s","%(url)s",'
    mquery += '"%(pdbid)s","%(chain)s",%(identity)f,"%(ostat)s","%(mthd)s")'
    mfields = dict((key,m.__getattribute__(key)) for key in dir(m) if isinstance(key,collections.Hashable))
    mfields["label"] = self.slabel
    mfields.update(remark3_metrics)
    mfields["identity"] = float(remark3_metrics["sid"])
    mquery = mquery%mfields
    # First, pass the underlying PDBMapStructure to upload_structure
    self.upload_structure(model=True,sfields=mfields)
    # Then execute the Model upload query
    self._connect()
    self._c.execute(mquery)
    self._close()


  # Prepare the Model summary information (if no errors occurred)
  def upload_genomic_data(self,dstream,dname):
    """ Uploads genomic data via a PDBMapData generator """
    filterwarnings('ignore', category = MySQLdb.Warning)
    i=0 # ensure initialization
    for i,record in enumerate(dstream):
      if not i%1000:
        sys.stdout.write("\rRecords uploaded: %5d"%i)
      # Upload all but the consequences and optional Fst information to GenomicData
      record.INFO['LABEL'] = dname
      query  = "INSERT IGNORE INTO GenomicData "
      query += "(label,chr,start,end,name,variation,vtype,svtype,ref_allele,alt_allele,"
      query += "svlen,quality,avgpost,rsq,erate,theta,ldaf,ac,an,aa,da,maf,amr_af,asn_af,"
      query += "eas_af,sas_af,afr_af,eur_af,ens_gene,hgnc_gene,"
      query += "snpsource,format,gt) VALUES "
      query += "(%(LABEL)s,%(CHROM)s,%(START)s,%(END)s,%(ID)s,%(EXISTING)s,%(VT)s,"
      query += "%(SVTYPE)s,%(REF)s,%(ALT)s,%(SVLEN)s,%(QUAL)s,%(AVGPOST)s,%(RSQ)s,"
      query += "%(ERATE)s,%(THETA)s,%(LDAF)s,%(AC)s,%(AN)s,%(AA)s,%(DA)s,"
      query += "%(AF)s,%(AMR_AF)s,%(ASN_AF)s,%(EAS_AF)s,%(SAS_AF)s,%(AFR_AF)s,%(EUR_AF)s,"
      query += "%(GENE)s,%(HGNC)s,%(SNPSOURCE)s,%(FORMAT)s,%(GT)s)"
      try:
        self._connect(cursorclass=MySQLdb.cursors.Cursor)
        self._c.execute(query,record.INFO)
        self._con.commit()
      except Exception as e:
        if "_last_executed" in dir(self._c):
          msg = self._c._last_executed.replace('\n',';')
          sys.stderr.write("WARNING (PDBMapIO) GenomicData query failed, query: %s\n"%msg)
        else:
          msg  = str(e).replace('\n',';')
          msg += "WARNING (PDBMapIO) GenomicData query failed, exception: %s\n"%msg
          sys.stderr.write(msg)
        self._con.rollback()
        continue # halt upload of this variant
      finally:
        self._close()

      # Upload each consequence to GenomicConsequence
      query  = "INSERT IGNORE INTO GenomicConsequence "
      query += "(label,chr,start,end,name,transcript,protein,uniprot,canonical,allele,"
      query += "consequence,cdna_pos,cds_pos,protein_pos,ref_amino_acid,"
      query += "alt_amino_acid,ref_codon,alt_codon,polyphen,sift,biotype,"
      query += "domains) VALUES "
      query += "(%(LABEL)s,%(CHROM)s,%(START)s,%(END)s,%(ID)s,"
      query += "%(Feature)s,%(ENSP)s,%(SWISSPROT)s,%(CANONICAL)s,%(Allele)s,"
      query += "%(Consequence)s,%(cDNA_position)s,%(CDS_position)s,"
      query += "%(Protein_position)s,%(Ref_AminoAcid)s,%(Alt_AminoAcid)s,"
      query += "%(Ref_Codon)s,%(Alt_Codon)s,%(PolyPhen)s,%(SIFT)s,"
      query += "%(BIOTYPE)s,%(DOMAINS)s)"
      for csq in record.CSQ:
        csq["LABEL"] = self.dlabel
        try: 
          self._connect(cursorclass=MySQLdb.cursors.Cursor)
          self._c.execute(query,csq)
          self._con.commit()
        except Exception as e:
          if "_last_executed" in dir(self._c):
            msg = self._c._last_executed.replace('\n',';')
            sys.stderr.write("WARNING (PDBMapIO) GenomicConsequence query failed, query:: %s\n"%msg)
          else:
            msg  = str(e).replace('\n',';')
            msg += "WARNING (PDBMapIO) GenomicConsequence query failed, exception: %s\n"%msg
            sys.stderr.write(msg)    
          self._con.rollback()
          continue # halt upload of this variant
        finally:
          self._close()
    resetwarnings()
    return i # return the number of uploaded rows

  def upload_intersection(self,dstream,buffer_size=1):
    """ Uploads an intersection via a process parser generator """
    # Query header
    queryh  = "INSERT IGNORE INTO GenomicIntersection "
    queryh += "(dlabel,slabel,structid,chain,seqid,gc_id) VALUES "
    # Query value string
    queryv  = "(%s,%s,%s,%s,%s,%s)" # Direct reference
    # Query value list
    vals = []
    # Query argument list
    args = []
    i=0 # ensure initialization
    for i,row in enumerate(dstream):
      row = list(row) # convert from tuple
      vals.append(queryv)
      args.extend(row)
      if not (i+1) % buffer_size:
        # Construct query string
        query = queryh + ','.join(vals)
        try:
          self._connect()
          self._c.execute(query,args)
          self._con.commit()
          vals,args = [],[]
        except:
          self._con.rollback()
          raise
        finally:
          self._close()
    # Upload any rows left in the buffer
    if args:
      # Construct query string
      query = queryh + ','.join(vals)
      try:
        self._connect()
        self._c.execute(query,args)
        self._con.commit()
      except:
        self._con.rollback()
        raise
      finally:
        self._close()
    return(i) # Return the number of uploaded rows

  def load_structure(self,structid,biounit=0,useranno=False,raw=False,syn=False):
    """ Reconstructs annotated PDBMapStructure from the database """
    biounit = 0 if biounit < 0 else biounit
    query = PDBMapIO.structure_query
    if useranno:
      supp_select = ",z.* "
      supp_table  = "LEFT JOIN pdbmap_supp.%s AS z ON z.chr=g.chr "%self.dlabel
      supp_table += "AND z.start=g.start AND z.name=g.name "
      try:
        # Most custom datasets joined on chromosomal position
        query = query%(supp_select,supp_table) # Insert genomic useranno text
        q = self.secure_query(query,qvars=(self.dlabel,self.slabel,
                        structid,biounit),cursorclass='DictCursor')
        q = [row for row in q] # must throw exception here if the query failed
      except:
        # But some are joined on structural position
        query = PDBMapIO.structure_query
        supp_select = ",z.* "
        supp_table  = "LEFT JOIN pdbmap_supp.%s AS z ON z.pdbid=a.structid "%self.dlabel
        supp_table += "AND z.chain=a.chain AND z.seqid=a.seqid "
        query = query%(supp_select,supp_table) # Insert structural useranno text
        q = self.secure_query(query,qvars=(self.dlabel,self.slabel,
                        structid,biounit),cursorclass='DictCursor')
    else:
      supp_select = ''
      supp_table  = ''
      query = query%(supp_select,supp_table) # Insert useranno text
      q = self.secure_query(query,qvars=(self.dlabel,self.slabel,
                      structid,biounit),cursorclass='DictCursor')
      
    # Return raw result dictionary if specified
    if raw:
      res = {}
      for row in q:
        if not res:
          res = dict([(key,[val]) for key,val in list(row.items())])
        else:
          for key,val in list(row.items()):
            res[key].append(val)
      return res
    # Construct and return as a Bio.PDB.Structure
    s = None
    for row in q:
      if not s:
        s = Bio.PDB.Structure.Structure(row['structid'])
        s.biounit  = row['biounit']
        s.resolution = row['resolution'] if row['resolution'] else None
        s.mpqs  = row['mpqs'] if row['mpqs'] else None
        s.type       = 'model' if row['mpqs'] else 'structure'
        s.method  = row['method']
        s.snpmapper  = {}
      m = row['model']
      if m not in s: s.add(Bio.PDB.Model.Model(m))
      m = s[m]
      c = row['chain']
      if c not in m: m.add(Bio.PDB.Chain.Chain(c))
      c = m[c]
      c.biounit = row['biounit']
      c.unp  = row['unp']
      c.transcript  = row['enst']
      c.gene  = row['ensg']
      c.perc_identity = row['perc_identity']
      c.hybrid  = row['hybrid']
      r = (' ',row['seqid'], ' ')
      c.add(Bio.PDB.Residue.Residue(r,aa_code_map[row['rescode']],' '))
      r = c[r]
      r.biounit = row['biounit']
      r.rescode = row['rescode']
      r.seqid = row['seqid']
      r.coords = (row['x'],row['y'],row['z'])
      r.pfam  = (row['pfamid'],row['pfam_domain'],row['pfam_desc'],row['pfam_evalue'])
      r.ss  = row['ss']
      r.angles = (row['phi'],row['psi'],row['tco'],row['kappa'],row['alpha'])
      if c.unp not in s.snpmapper:
        s.snpmapper[c.unp] = np.array([[None for j in range(19)] for i in range(r.seqid)])
      if not np.any(s.snpmapper[c.unp][r.seqid-1]):
        snp  = [row['issnp'],row['snpid'],row['chr'],row['start'],row['end'],row['hgnc_gene'],row['ens_gene']]
        snp += [row['anc_allele'],row['ref_allele'],row['alt_allele']]
        snp += [row['maf'],row['amr_af'],row['asn_af'],row['eas_af'],row['sas_af'],row['eur_af'],row['afr_af']]
        snp += [row['ref_codon'],row['alt_codon'],row['vep_ref_aa'],row['vep_alt_aa']]
        s.snpmapper[c.unp][r.seqid-1] = snp
      r.snp = s.snpmapper[c.unp][r.seqid-1]
    return PDBMapStructure(s)

  def load_unp(self,unpid):
    """ Identifies all associated structures and models, then pulls those structures. """
    # Queries all structures and models associated with a given UniProt ID
    query = """SELECT DISTINCT structid FROM Chain WHERE label=%s AND (unp='%s' OR unp LIKE '%s-%%');"""
    base_unp = unpid.split('-')[0]
    q = self.secure_query(query,qvars=('pdb',base_unp,base_unp),
                                        cursorclass='Cursor')
    entities = [r[0] for r in q]
    q = self.secure_query(query,qvars=('modbase',base_unp,base_unp),
                                        cursorclass='Cursor')
    entities.extend([r[0] for r in q])
    q = self.secure_query(query,qvars=('swiss',base_unp,base_unp),
                                        cursorclass='Cursor')
    entities.extend([r[0] for r in q])
    LOGGER.info("%s found in %d structures/models."%(unpid,len(entities)))
    res = []
    for entity in entities:
      entity_type = self.detect_entity_type(entity)
      if entity_type == 'structure':
        res.append((entity_type,entity))
      elif entity_type == 'model':
        res.append((entity_type,entity))
      elif entity_type == 'swiss':
        res.append((entity_type,entity))
    return res

  @classmethod
  def dssp(cls,m,fname,dssp_executable='dssp'):
    dssp_executable = '/dors/capra_lab/projects/psb_collab/psb_pipeline/data/dssp/dssp_local.exe'
    # Create temp file containing only ATOM rows
    try:
      temp_fileobject = tempfile.NamedTemporaryFile(delete=False)
      tmpf = temp_fileobject.name
    except Exception as e:
      raise Exception("Unable to create temporary file for dssp processing:\n%s"%str(e));
    temp_fileobject.close();

    cmd  = "grep ^ATOM %s > %s"%(fname,tmpf)
    try:
      if fname.split('.')[-1] == 'gz':
        cmd = 'z'+cmd
      os.system(cmd)
      cmd  = [dssp_executable,tmpf]
      # A problem with this code is that stderr coming from dssp messes up log badly
      p = sp.Popen(cmd,stdout=sp.PIPE)
      resdict = {}
      start_flag = False
      stdout,stderr = p.communicate()
      for i,linebytes in enumerate(stdout.splitlines()):
        if not linebytes: continue
        line = linebytes.decode('latin')
        if not start_flag:
          if line.startswith('  #'):
            start_flag = True
          continue
        resinfo = {}
        if not line[5:10].strip(): continue
        resinfo['seqid'] = int(line[5:10])
        resinfo['icode'] = line[10:11]
        resinfo['chain'] = line[11:12].upper()
        if not resinfo['chain'].strip():
          resinfo['chain'] = 'A'
        resinfo['aa']    = line[13:14].upper()
        resinfo['ss']    = line[16:17].upper()
        resinfo['acc']   = float(line[34:38])
        resinfo['tco']   = float(line[85:91])
        resinfo['kappa'] = float(line[92:97])
        resinfo['alpha'] = float(line[98:103])
        resinfo['phi']   = float(line[103:109])
        resinfo['psi']   = float(line[109:115])
        resinfo['rsa']   = resinfo['acc'] / solv_acc[resinfo['aa']]
        resdict[(resinfo['chain'],resinfo['seqid'],resinfo['icode'])] = resinfo.copy()
    except Exception as e:
      msg = "ERROR (PDBMapIO) Unable to parse output from DSSP invocation: %s\nException:%s"%(str(cmd),str(e))
      raise Exception(msg)
    finally:        
      # Remove the temp ATOM file
      try:
        os.remove(tmpf)
      except Exception as e:
        raise Exception("Unable to delete dssp temporary input file %s\n:Exception:",str(tmpf),str(e))
    return resdict

  def detect_entity_type(self,entity):
    """ Given an entity ID, attempts to detect the entity type """
    if self.model_in_db(entity,label=None):
      return "model"
    elif self.swiss_in_db(entity,label=None):
      return "swiss"
    elif self.structure_in_db(entity,label=None):
      return "structure"
    elif PDBMapProtein.isunp(entity):
      return "unp"
    elif self.unp_in_db(entity,label=None):
      return "unp"
    # If a gene, return the associated UniProt ID
    elif PDBMapProtein.ishgnc(entity):
      return PDBMapProtein.hgnc2unp(entity)
    else:
      isgene,unp = self.gene_in_db(entity,label=None)
      if isgene:
        return unp
    return None

  def load_sifts(self,fdir,conf_file):
    """ Given a SIFTS XML directory, uploads to PDBmap """
    cmd = "scripts/sifts_parser.py -c %s %s"%(conf_file,fdir)
    os.system(cmd)
    return -1
    ## This code used to have to do range-inferences to derive
    ## 1-to-1 mappings from start-end mappings. We're now using
    ## the 1-to-1 XML files from SIFTS, so this preprocessing
    ## is no longer necessary and the data can be loaded directly
    ## from the XML file via a simple sifts XML parser.
    # if sprot:
    #   PDBMapProtein.load_sprot(sprot)
    #   humansp = PDBMapProtein.sprot
    # else:
    #   humansp = None
    # rc = 0
    # with open(fname,'rb') as fin:
    #   fin.readline(); fin.readline() # skip first two lines
    #   reader = csv.reader(fin,delimiter='\t')
    #   for row in reader:
    #     q  = "INSERT IGNORE INTO sifts "
    #     q += "(pdbid,chain,sp,pdb_seqid,sp_seqid) VALUES "
    #     v = []
    #     pdbid,chain,sp = row[0:3]
    #     pdbid = pdbid.strip()
    #     if humansp and row[2] not in humansp:
    #       continue # not a human protein
    #     try:
    #       res_beg,res_end,pdb_beg,pdb_end,sp_beg,sp_end = [int(x) for x in row[3:]]
    #     except:
    #       # msg  = "WARNING (PDBMapIO) SIFTS icode error: "
    #       # msg += "%s,%s,%s\n"%(pdbid,chain,sp)
    #       # sys.stderr.write(msg)
    #       continue
    #     res_range = range(res_beg,res_end+1)
    #     pdb_range = range(pdb_beg,pdb_end+1)
    #     sp_range  = range(sp_beg,sp_end+1)
    #     if len(res_range) != len(pdb_range) or len(pdb_range) != len(sp_range):
    #       # msg  = "WARNING (PDBMapIO) SIFTS range mismatch: "
    #       # msg += "%s,%s,%s\n"%(pdbid,chain,sp)
    #       # sys.stderr.write(msg) 
    #       continue
    #     for i,seqid in enumerate(pdb_range):
    #       v.append("('%s','%s','%s',%d,%d)"%(pdbid,chain,sp,seqid,sp_range[i]))
    #     # Upload after each row
    #     q = q+','.join(v)
    #     rc += self.secure_command(q)
    # return rc

  def load_pfam(self,fname):
    query  = "LOAD DATA LOCAL INFILE %s "
    query += "INTO TABLE pfam "
    query += "FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n' "
    query += "IGNORE 1 LINES"
    rc = self.secure_command(query,(fname,))
    return rc

  def show(self):
    return(self.__dict__['structure'])

  def _alive(self,con):
    try:
      c = con.cursor()
      c.execute("SELECT 1")
      return True
    except:
      return False

  def _connect(self,usedb=True,cursorclass=MySQLdb.cursors.DictCursor):
    # Finally, open a cursor from the active connection
    # Search for existing connection with desired cursorclass
    self._con = None
    for con in self._cons:
      if con.cursorclass == cursorclass:
        if self._alive(con):
          # Set as active connection
          self._con = con
        else:
          # Remove dead connection
          LOGGER.info("Calling con.close() to close SQL connection")
          con.close()
          self._cons.remove(con)
    # If none was found, open a new connection with desired cursorclass
    if not self._con:
      trycount = 0
      nMaxRetries = 100000
      while (True):
        trycount += 1
        LOGGER.info("Connecting to %s with %s"%(self.dbhost,cursorclass.__name__))
        try:
          if usedb:
            if self.dbpass:
              con = MySQLdb.connect(host=self.dbhost,
                    user=self.dbuser,passwd=self.dbpass,
                    db=self.dbname,cursorclass=cursorclass)
            else:
              con = MySQLdb.connect(host=self.dbhost,
                    user=self.dbuser,db=self.dbname,
                    cursorclass = cursorclass)
          else:
            if self.dbpass:
              con = MySQLdb.connect(host=self.dbhost,
                    user=self.dbuser,passwd=self.dbpass,
                    cursorclass = cursorclass)
            else:
              con = MySQLdb.connect(host=self.dbhost,
                  user=self.dbuser,cursorclass=cursorclass)
          # Add to connection pool and set as active connection
          # This point is "success" from MySQLdb.connect
          # LOGGER.info("Successful connect to %s"%self.dbhost)
          self._cons.append(con)
          self._con = con
          # Success... Leave while loop
          break
        except MySQLdb.OperationalError as e:
          if e[0] == 1040:          
            msg = "Unable to connect to connect to mySQL: %s\n"%e
            LOGGER.exception(msg)
            sys.stderr.write(msg)
            # If nMaxRetries, terminate with exception
            if (trycount == nMaxRetries):
              msg = "Max retries of %d reached\n"%nMaxRetries
              logging.exception(msg)
              raise Exception(msg)
            else:
              waitTime = 30*(1+trycount % 5)
              msg = "Try %d failed.  Waiting %d secs and retrying...\n"%(trycount,waitTime)
              LOGGER.exception(msg)
              time.sleep(waitTime) # Wait 30/60/90 etc seconds and retry (return to top of while loop)
            # Loop back through while loop
          else:
            # If nMaxRetries, terminate with exception
            msg  = "\nError database connection unsuccessful: %s after %d retries\n"%(e,trycount)
            msg += "Parameters:\n"
            msg += " DBHOST = %s\n"%self.dbhost
            msg += " DBNAME = %s\n"%self.dbname
            msg += " DBUSER = %s\n\n"%self.dbuser
            logging.exception(msg)
            raise
    # Finally, open a cursor from the active connection
    self._c = self._con.cursor() # and open a new one
    return self._c

  def _close(self):
    # LOGGER.info("Closing SQL connection")
    self._c.close() # Close only the cursor
    return

  def secure_query(self,query,qvars=None,cursorclass='SSDictCursor'):
    """ Executes queries using safe practices """
    logmsg = "With cursorclass %s, performing SQL query:\n"%str(cursorclass)
    if qvars:
      logmsg += query%qvars
    else:
      logmsg += query
    LOGGER.info(logmsg)

    if cursorclass == 'SSDictCursor':
      self._connect(cursorclass=MySQLdb.cursors.SSDictCursor)
    elif cursorclass == 'SSCursor':
      self._connect(cursorclass=MySQLdb.cursors.SSCursor)
    elif cursorclass == 'Cursor':
      self._connect(cursorclass=MySQLdb.cursors.Cursor)
    else:
      self._connect()
    filterwarnings('ignore', category = MySQLdb.Warning)
    try:
      if qvars: self._c.execute(query,qvars)
      else:     self._c.execute(query)
      for row in self._c:
        yield row
    except Exception as e:
      msg  = "ERROR (PDBMapIO) Secure query failed; Exception: %s; "%str(e)
      msg += " Provided query: %s"%query
      msg += " Provided args: %s"%str(qvars)
      if "_last_executed" in dir(self._c):
        msg += "\n Executed Query: \n%s"%self._c._last_executed
      LOGGER.exception(msg)
      raise Exception(msg)
    finally:
      # msg = "Executed Query: \n%s\n"%self._c._last_executed
      # print msg
      resetwarnings()
      self._close()


  def secure_cached_query(self,cache_dir,query,qvars=None,cursorclass='SSDictCursor'):
    """ Executes queries using safe practices, and loads/stores results in a cache_directory """
    logmsg = "With cursorclass %s, performing SQL query:\n"%str(cursorclass)
    query_hash = cursorclass 
    if qvars:
      query_hash += hashlib.sha256((query%qvars).encode('utf-8')).hexdigest()
      logmsg += query%qvars
    else:
      query_hash += hashlib.sha256(query.encode('utf-8')).hexdigest()
      logmsg += query
    
    query_hash_filename = os.path.join(cache_dir,query_hash)
    logmsg += "\nquery_hash = %s"%query_hash_filename
    LOGGER.info(logmsg)
    if os.path.exists(query_hash_filename):
      LOGGER.info("returning query results from cache")
      with open(query_hash_filename,'rb') as cache_file_handle:
        rows = pickle.load(cache_file_handle)
        for row in rows:
          yield row
    else: 
      LOGGER.info("Cache not found.  Executing new query")
      if cursorclass == 'SSDictCursor':
        self._connect(cursorclass=MySQLdb.cursors.SSDictCursor)
      elif cursorclass == 'SSCursor':
        self._connect(cursorclass=MySQLdb.cursors.SSCursor)
      elif cursorclass == 'Cursor':
        self._connect(cursorclass=MySQLdb.cursors.Cursor)
      else:
        self._connect()
      filterwarnings('ignore', category = MySQLdb.Warning)
      rows = []
      try:
        if qvars: self._c.execute(query,qvars)
        else:     self._c.execute(query)
        for row in self._c:
          rows.append(row)
      except Exception as e:
        msg  = "ERROR (PDBMapIO) Secure query failed; Exception: %s; "%str(e)
        msg += " Provided query: %s"%query
        msg += " Provided args: %s"%str(qvars)
        if "_last_executed" in dir(self._c):
          msg += "\n Executed Query: \n%s"%self._c._last_executed
        LOGGER.exception(msg)
        raise Exception(msg)
      finally:
        # msg = "Executed Query: \n%s\n"%self._c._last_executed
        # print msg
        resetwarnings()
        self._close()
      
      with tempfile.NamedTemporaryFile(dir=cache_dir) as cache_file_handle:
        pickle.dump(rows,cache_file_handle)
        try: # Atomic rename can fail due to multi-process contention.  don't halt for that!
          os.link(cache_file_handle.name,query_hash_filename)
        except:
          LOGGER.warn("Unable to rename %s to %s"%(cache_file_handle.name,query_hash_filename))
        else:
          umask = os.umask(0)
          os.umask(umask)
          os.chmod(query_hash_filename, 0o666 & ~umask)
          LOGGER.info("Query results saved to %s"%query_hash_filename)
      for row in rows:
        yield row

  def secure_command(self,query,qvars=None):
    """ Executes commands using safe practices """
    self._connect()
    filterwarnings('ignore', category = MySQLdb.Warning)
    try:
      if qvars: self._c.execute(query,tuple(qvars))
      else:     self._c.execute(query)
      self._con.commit()
      rc = self._c.rowcount

    except Exception as e:
      self._con.rollback()
      msg  = "ERROR (PDBMapIO) Secure command failed; Exception: %s; "%str(e)
      msg += "\nCommand: %s"%self._c._last_executed
      raise Exception(msg)
    finally:
      self._close()
      resetwarnings()
    return rc

  def check_schema(self):
    return  # I have no idea why we'd want to do this with every startup
    self._connect(usedb=False)
    format_dict = {'dbuser':self.dbuser,'dbhost':self.dbhost}
    queries = [ 'lib/create_schema_Idmapping.sql',
                'lib/create_schema_Structure.sql',
                'lib/create_schema_Model.sql',
                'lib/create_schema_Chain.sql',
                'lib/create_schema_Residue.sql',
                'lib/create_schema_Transcript.sql',
                'lib/create_schema_AlignChainTranscript.sql',
                'lib/create_schema_GenomicData.sql',
                'lib/create_schema_GenomicConsequence.sql',
                'lib/create_schema_GenomicIntersection.sql',
                'lib/create_schema_sifts.sql',
                # 'lib/create_schema_pfam.sql',
                'lib/create_procedure_assign_foreign_keys.sql',
                'lib/create_procedure_get_protein.sql',
                'lib/create_procedure_get_structure.sql',
                'lib/create_procedure_get_full_structure.sql']
                # 'lib/create_procedure_get_repr_subset.sql']
    filterwarnings('ignore', category = MySQLdb.Warning)
    try: # Create the database if not exists
      self._c.execute("CREATE DATABASE %s"%self.dbname)
    except: pass # Database exists. Using existing.
    finally: self._c.execute("USE %s"%self.dbname)
    for sql_file in queries:
      with open(sql_file) as sql_fp:
        lines = [line.split('#')[0].rstrip() for line in sql_fp]
        try: # Create the tables and procedures if not exists
          query = ' '.join(lines)
          self._c.execute(query%format_dict)
        except: pass # Table exists. Using existing.
    resetwarnings()
    self._close()

  # Query definitions
 
  # Queries all biological assemblies and models with the given structid
  structure_query = """SELECT
  /*SLABEL*/a.label as slabel,
  /*ID*/a.structid,
  /*METHOD*/IF(c.method IS NULL,d.method,c.method) as method,
  /*STRUCTURE*/c.resolution,
  /*MODEL*/d.mpqs,
  /*CHAIN*/b.biounit,b.model,b.chain,b.unp,b.hybrid,
  /*RESIDUE*/a.seqid,a.icode,a.rescode,a.ss,a.phi,a.psi,a.tco,a.kappa,a.alpha,a.x,a.y,a.z,
  /*PFAM*/p.acc as pfamid,p.name as pfam_domain,p.description as pfam_desc,p.evalue as pfam_evalue,
  /*ALIGN*/h.chain_seqid as aln_chain_seqid,h.trans_seqid as aln_trans_seqid,j.perc_identity,
  /*TRANS*/i.transcript as enst,i.gene as ensg,i.protein as ens_prot,i.seqid as ens_prot_seqid,i.rescode as ens_prot_aa,
  /*VARIANT*/IF(f.consequence LIKE '%%%%missense_variant%%%%',1,0) as issnp,
  /*DLABEL*/g.label as dlabel,
  /*VARIANT*/g.name as snpid,g.chr,g.start,g.end,g.hgnc_gene,g.ens_gene,g.aa as anc_allele,g.ref_allele,g.alt_allele,g.maf,g.amr_af,g.asn_af,g.eas_af,g.sas_af,g.eur_af,g.afr_af,
  /*CONSEQUENCE*/f.gc_id,f.transcript as vep_trans,f.protein as vep_prot,f.uniprot as unp,f.protein_pos as vep_prot_pos,f.ref_codon,f.alt_codon,f.ref_amino_acid as vep_ref_aa,f.alt_amino_acid as vep_alt_aa,
  /*CONSEQUENCE*/f.consequence,f.polyphen,f.sift,f.biotype
  /*USERANNO*/%s
  FROM Residue as a
  LEFT JOIN Chain as b
  ON a.label=b.label AND a.structid=b.structid AND a.biounit=b.biounit AND a.model=b.model AND a.chain=b.chain
  LEFT JOIN Structure as c
  ON b.label=c.label AND b.structid=c.pdbid
  LEFT JOIN Model as d
  ON b.label=d.label AND b.structid=d.modelid
  LEFT JOIN GenomicIntersection as e
  ON a.label=e.slabel AND a.structid=e.structid AND a.chain=e.chain AND a.seqid=e.seqid AND e.dlabel=%%s
  LEFT JOIN GenomicConsequence as f
  ON e.dlabel=f.label AND e.gc_id=f.gc_id
  LEFT JOIN GenomicData as g
  ON f.label=g.label AND f.chr=g.chr AND f.start=g.start AND f.end=g.end AND f.name=g.name
  LEFT JOIN Alignment as h USE INDEX(PRIMARY)
  ON a.label=h.label AND a.structid=h.structid AND a.chain=h.chain AND a.seqid=h.chain_seqid #AND f.transcript=h.transcript
  LEFT JOIN Transcript as i USE INDEX(PRIMARY)
  ON h.label=i.label AND h.transcript=i.transcript AND h.trans_seqid=i.seqid
  LEFT JOIN AlignmentScore as j
  ON h.label=j.label AND h.structid=j.structid AND h.chain=j.chain AND h.transcript=j.transcript
  LEFT JOIN pfam as p
  ON a.structid=p.pdbid AND a.chain=p.chain AND a.seqid BETWEEN p.seqstart AND p.seqend
  /*USERANNO*/%s
  WHERE a.label=%%s
  AND (f.transcript IS NULL OR f.transcript="" OR f.transcript=h.transcript)
  AND a.structid=%%s AND a.biounit=%%s
  ORDER BY b.unp,a.seqid DESC;
  """
  # Removing non-asymmetric unit clause since biounit is always specified explicitly anyway
  # and sometimes the asymmetric unit is wanted. Retaining below:
  # AND (c.method LIKE '%%%%nmr%%%%' OR b.biounit>0 OR NOT ISNULL(d.modelid))


aa_code_map = {"ala" : "A",
        "arg" : "R",
        "asn" : "N",
        "asp" : "D",
        "asx" : "B",
        "cys" : "C",
        "glu" : "E",
        "gln" : "Q",
        "glx" : "Z",
        "gly" : "G",
        "his" : "H",
        "ile" : "I",
        "leu" : "L",
        "lys" : "K",
        "met" : "M",
        "phe" : "F",
        "pro" : "P",
        "ser" : "S",
        "thr" : "T",
        "trp" : "W",
        "tyr" : "Y",
        "val" : "V"}
for key,val in list(aa_code_map.items()):
  aa_code_map[val] = key        
solv_acc = {"A" : 115,
        "R" : 225,
        "N" : 150,
        "D" : 160,
        "B" : 155,
        "C" : 135,
        "E" : 190,
        "Q" : 180,
        "Z" : 185,
        "G" : 75,
        "H" : 195,
        "I" : 175,
        "L" : 170,
        "K" : 200,
        "M" : 185,
        "F" : 210,
        "P" : 145,
        "S" : 115,
        "T" : 140,
        "W" : 255,
        "Y" : 230,
        "V" : 155,
        "X" : 180} # Undetermined amino acid SA taken from Sander 1994

## Copied from biolearn
def multidigit_rand(digits):
  randlist = [random.randint(1,10) for i in range(digits)]
  multidigit_rand = int(''.join([str(x) for x in randlist]))
  return multidigit_rand

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
