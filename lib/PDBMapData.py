#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : PDBMapData.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-18
# Description    : Defines the PDBMapData class. This class is not intended
#                : to store or retain any data. It's purpose is to manage the 
#                : upload of annotation files to the PDBMap.GenomicData 
#                : database.
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv,re,gzip
import numpy as np
import subprocess as sp
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import vcf # PyVCF
from lib import bed # PyVCF emulator for BED files

class PDBMapData():
  
  def __init__(self,vep="",dname=''):
    self.vep   = vep
    self.dname = dname
    if self.vep and not os.path.exists(self.vep):
      msg = "ERROR (PDBMapData) VEP location invalid: %s"%vep
      raise Exception(msg)
    if self.vep:
      # Check for a dbconn file
      cache    = "/dors/capra_lab/data/vep/"
      # cache    = os.path.expanduser(cache) # replace ~ with explicit home directory
      if not os.path.exists(cache):
        msg = "WARNING (PDBMapData) No cache exists. Using network connection.\n"
        sys.stderr.write(msg)
        cache = None
      # Construct the VEP command
      self.vep_cmd = [self.vep,'-i','']
      # Specify the input type
      self.vep_cmd.extend(['--format',''])
      # Use the local VEP cache
      self.vep_cmd.extend(['--cache','--offline','--dir',cache])
      # Disable the progress bars
      self.vep_cmd.extend(['--no_progress'])
      # Increase buffer size to improve runtime (default 5,000)
      self.vep_cmd.extend(['--buffer_size','500000']) # ~5GB per 100,000
      # Annotate with functional info/prediction
      self.vep_cmd.extend(['--sift','s','--polyphen','s'])
      # Annotate with variant, gene, protein, and domain identifiers
      self.vep_cmd.extend(['--check_existing','--symbol','--protein','--uniprot','--domains'])
      # Annotate transcript with canonical bool and biotype
      self.vep_cmd.extend(['--canonical','--biotype','--pubmed'])
      ## We are now allowing Synonymous SNPs to be mapped ##
      # Retain only coding-region variants
      self.vep_cmd.extend(['--coding_only'])
      # Specify output format
      self.vep_cmd.extend(['--vcf'])
      # Send to stdout and don't generate a summary stats file
      self.vep_cmd.extend(['--no_stats'])
      self.vep_cmd.extend(['-o','stdout'])

  def record_parser(self,record,info_headers,csq_headers):

    # Ensure that an end is specified, default to start+1
    if "END" not in record.INFO:
      record.INFO["END"] = int(record.POS) + 1

    # Ensure that each record contains all fields
    for header in info_headers:
      if header not in record.INFO:
        record.INFO[header] = None

    # Process fields which may or may not be included in VEP output
    if 'SNPSOURCE' in record.INFO and record.INFO['SNPSOURCE']:
      record.INFO['SNPSOURCE'] = ';'.join(record.INFO['SNPSOURCE'])
    else: record.INFO['SNPSOURCE'] = None

    # Make any necessary population allele frequency corrections
    # Enforce biallelic assumption
    freqs = ['AMR_AF','AFR_AF','EUR_AF','EAS_AF','SAS_AF','ASN_AF']
    for freq in freqs:
      if freq not in record.INFO: record.INFO[freq] = None # Ensure all fields present
      record.INFO[freq] = 0.0 if not record.INFO[freq] else record.INFO[freq]
      if type(record.INFO[freq]) in [type((None,)),type([])]:
        record.INFO[freq] = float(record.INFO[freq][0])

    # Ensure the ancestral allele is encoded properly
    if 'AA' in record.INFO and record.INFO['AA']:
      # Format: AA|REF|ALT|INDELTYPE
      record.INFO['AA'] = record.INFO['AA'].split('|')[0].upper()
    else: 
      record.INFO['AA'] = None

    # Extract the format and genotypes at this site
    record.INFO['FORMAT'] = record.FORMAT
    record.INFO['GT'] = ','.join([str(s['GT']) for s in record.samples])

    # Enforce biallelic assumption
    # Record only the first alternate allele count
    record.INFO['AC'] = record.INFO['AC'][0] if 'AC' in record.INFO and record.INFO['AC'] else None

    # Ensure necessary fields are present in record.INFO
    if 'AN'      not in record.INFO: record.INFO['AN']      = None
    if 'AF'      not in record.INFO: record.INFO['AF']      = None
    if 'VT'      not in record.INFO: record.INFO['VT']      = None
    if 'SVTYPE'  not in record.INFO: record.INFO['SVTYPE']  = None
    if 'SVLEN'   not in record.INFO: record.INFO['SVLEN']   = None
    if 'AVGPOST' not in record.INFO: record.INFO['AVGPOST'] = None
    if 'RSQ'     not in record.INFO: record.INFO['RSQ']     = None
    if 'ERATE'   not in record.INFO: record.INFO['ERATE']   = None
    if 'THETA'   not in record.INFO: record.INFO['THETA']   = None
    if 'LDAF'    not in record.INFO: record.INFO['LDAF']    = None

    # Reformat the variant-type value
    if type(record.INFO["VT"]) == type(tuple()):
      record.INFO["VT"] = ','.join(vt for vt in list(record.INFO["VT"]))
    
    # Ensure 1000 Genomes Fst fields are populated or None
    popfst  = ['AMREASFST','AMRSASFST','AMREURFST','AMRAFRFST','EASSASFST']
    popfst += ['EASEURFST','EASAFRFST','SASEURFST','SASAFRFST','EURAFRFST']
    popfst += ['AMREASSASEURAFRFST'] # all-populations Fst
    for pop in popfst:
      if pop in record.INFO:
        nhat,dhat,fst = record.INFO[pop]
        pop = pop[:-3] # Remove the FST suffix
        pop = pop if pop != 'AMREASSASEURAFR' else 'ALLPOP'
        if nhat=='nan' or dhat=='nan':
          print "nhat/dhat is string nan:",nhat
        if np.isnan(nhat) or np.isnan(dhat):
          print "nhat/dhat is numpy nan:",nhat,dhat
        record.INFO["%s_Nhat"%pop] = nhat if not np.isnan(nhat) else None
        record.INFO["%s_Dhat"%pop] = dhat if not np.isnan(dhat) else None
        record.INFO["%s_Fst"%pop]  = fst  if not np.isnan(fst)  else None
      else:
        pop = pop[:-3] # Remove the FST suffix
        pop = pop if pop != 'AMREASSASEURAFR' else 'ALLPOP'
        record.INFO["%s_Nhat"%pop] = None
        record.INFO["%s_Dhat"%pop] = None
        record.INFO["%s_Fst"%pop]  = None

    # Allele frequency is sometimes reecorded as a tuple or list
    # Enforce biallelic assumption
    # Extract the first (only) element and cast to float
    if type(record.INFO['AF']) in [type((None,)),type([])]:
      record.INFO['AF'] = float(record.INFO['AF'][0])

    # Add attribute fields to INFO
    record.INFO["ID"]     = record.ID
    record.INFO["CHROM"]  = record.CHROM
    record.INFO["START"]  = int(record.POS)
    record.INFO["REF"]    = record.REF
    record.INFO["ALT"]    = record.ALT[0] # Enforce biallelic assumption
    record.INFO["QUAL"]   = record.QUAL
    record.INFO["FILTER"] = str(record.FILTER)

    # Infer the derived allele from the ancestral
    # If ancestral not specified, assume minor allele is derived
    if record.INFO["ALT"] == record.INFO["AA"]:
      record.INFO["DA"] = record.INFO["REF"]
    else:
      record.INFO["DA"] = record.INFO["ALT"]      
    
    # If CSQ missing for this record or all records, dummy code
    if not csq_headers or not record.INFO["CSQ"]:
      record.CSQ = [{"ID":   record.INFO["ID"],
                    "CHROM": record.INFO["CHROM"],
                    "START": record.INFO["START"],
                    "END":   record.INFO["END"],
                    "Feature":          None,
                    "ENSP":             None,
                    "CANONICAL":        None,
                    "Allele":           None,
                    "Consequence":      None,
                    "cDNA_position":    None,
                    "CDS_position":     None,
                    "Protein_position": None,
                    "Ref_AminoAcid":    None,
                    "Alt_AminoAcid":    None,
                    "Ref_Codon":        None,
                    "Alt_Codon":        None,
                    "PolyPhen":         None,
                    "SIFT":             None,
                    "BIOTYPE":          None,
                    "DOMAINS":          None,
                    "SWISSPROT":        None}]
      record.INFO['EXISTING']  = None
      record.INFO["GENE"]      = None
      record.INFO["HGNC"]      = None
      record.INFO["SWISSPROT"] = None
    else:
      # Extract the consequences
      record.CSQ = self._parse_csq(csq_headers,record.INFO['CSQ'])
      # Add some "consequence" info to the variant info
      record.INFO["EXISTING"] = record.CSQ[0]['Existing_variation']
      record.INFO["GENE"] = record.CSQ[0]['Gene']
      if record.CSQ[0]['SYMBOL_SOURCE'] == 'HGNC':
        record.INFO["HGNC"] = record.CSQ[0]["SYMBOL"]
      else: record.INFO["HGNC"] = None
      # Correct missing variant IDs
      if not record.INFO["ID"]:
        if record.INFO["EXISTING"]:
          record.INFO["ID"] = record.INFO["EXISTING"]
        else:
          record.INFO["ID"] = "%(CHROM)s:%(START)d-%(END)d"%(record.INFO)
      # Add some variant info to the consequences
      for csq in record.CSQ:
        csq["ID"] = record.INFO["ID"]
        csq["CHROM"] = record.INFO["CHROM"]
        csq["START"] = record.INFO["START"]
        csq["END"] = record.INFO["END"]
    return record

  def load_vcf(self,fname,vep=True):
    """ Pipe VCF file through VEP and yield VEP rows """
    print "Initializing VCF generator."
    # Parse the VCF output from VEP, save to cache
    cache  = 'data/cache/%s'%os.path.basename(fname)
    if cache.split('.')[-1] != 'gz':
      cache += '.gz'
    if vep:
      parser = vcf.Reader(self.load_vep(fname,'vcf',cache),prepend_chr=True)
    else:
      parser = vcf.Reader(filename=fname,prepend_chr=True)
      if "CSQ" not in parser.infos or not parser.infos['CSQ']: # CSQ may be pre-computed
        # If no CSQ present, insert dummy structure
        parser.infos['CSQ'] = bed.Consequence() # dummy object
    # Determine Info headers
    info_headers = parser.infos.keys()
    # Determine Consequence headers
    csq_headers  = parser.infos['CSQ'].desc.split(': ')[-1].split('|')
    snpcount = 0
    for record in parser:
      # Do not load non-PASS variants
      if not record.FILTER:
        # If position refers to a haplotype chromosome, ignore
        if record.CHROM in ['chr1','chr2','chr3',
                          'chr4','chr5','chr6','chr7','chr8',
                          'chr9','chr10','chr11','chr12','chr13',
                          'chr14','chr15','chr16','chr17','chr18',
                          'chr19','chr20','chr21','chr22',
                          'chrX','chrY','chrMT']:
          snpcount += 1
          yield self.record_parser(record,info_headers,csq_headers)
    ## We are now allowing Synonymous SNPs to be mapped ##
    print "\nTotal SNPs (syn+nonsyn) in %s: %d"%(fname,snpcount)
    # print "Nonsynonymous SNPs in %s: %d"%(fname,nscount)

  def load_pedmap(self,fname):
    """ Convert PED/MAP to VCF, pipe VCF through VEP and load VEP output """
    print "load_pedmap not implemented"

  def load_bed(self,fname,id_type="id",vep=True,indexing=None):
    """ Load data from BED """
    print "Initializing BED generator."
    # Parse the VCF output from VEP, save to cache
    cache  = 'data/cache/%s.gz'%os.path.basename(fname)
    if cache.split('.')[-1] != 'gz':
      cache += '.gz'
    if vep:
      parser = vcf.Reader(self.load_vep(fname,id_type,cache),prepend_chr=True)
    else:
      parser = bed.Reader(fname,indexing=indexing,prepend_chr=True)
    # Determine Info headers
    info_headers = parser.infos.keys()
    # Determine Consequence headers
    if vep:
      csq_headers  = parser.infos['CSQ'].desc.split(': ')[-1].split('|')
    else: csq_headers = []

    snpcount = 0
    for record in parser:
      snpcount += 1
      # If position refers to a haplotype chromosome, ignore
      if record.CHROM not in ['chr1','chr2','chr3',
                        'chr4','chr5','chr6','chr7','chr8',
                        'chr9','chr10','chr11','chr12','chr13',
                        'chr14','chr15','chr16','chr17','chr18',
                        'chr19','chr20','chr21','chr22',
                        'chrX','chrY','chrMT']:
        continue
      yield self.record_parser(record,info_headers,csq_headers)
    print "Total SNPs (syn+nonsyn) in %s: %d"%(fname,snpcount)

  def load_vcffile(self,fname,io,buffer_size=1):
    """ Creates a supplementary table for original VCF datafile """
    parser = vcf.Reader(filename=fname,prepend_chr=True)
    # Extract header fields
    var_header  = ["chr","start","name","ref","alt","qual","filter"]
    # Types for the standard fields
    header_types  = ["VARCHAR(100)","BIGINT","VARCHAR(100)","VARCHAR(100)"]
    header_types += ["VARCHAR(100)","DOUBLE","VARCHAR(100)"]
    # Correct any overlaps with the default header
    for name in parser.infos.keys():
      if name.lower() in var_header:
        parser.infos["%s2"%name] = parser.infos[name]
        del parser.infos[name]
    # Extract info fields
    info_header = parser.infos.keys()
    # Extract and convert info types
    type_conv = {"Integer":"BIGINT","Float":"DOUBLE","Flag":"TINYINT","String":"TEXT"}
    info_types  = [type_conv[info.type] for info in parser.infos.values()]
    for char in ['.',' ','/','\\','(',')','[',']','-','!','+','=']:
        info_header = [f.replace(char,'_') for f in info_header]
    header = var_header + info_header #+ csq_header
    # Use the first row to infer data types for the INFO fields
    record = parser.next()
    types    = header_types + info_types #+ csq_types
    # Set default values for each type
    defaults = {"BIGINT":0,"DOUBLE":0.0,"TINYINT":0,"TEXT":"''","VARCHAR(100)":"''"}
    notnull  = set(["TINYINT","TEXT","VARCHAR(100)","DOUBLE"])
    # don't worry about type conversion for now, Nones are causing an issue
    formatter = {"BIGINT":"%s","DOUBLE":"%s","TEXT":"%s","TINYINT":"%s","VARCHAR(100)":"%s"}
    queryf = "(%s)"%','.join([formatter[f] for f in types])
    # Generate a create statement for this table
    table_def = ["`%s` %s DEFAULT %s %s"%(header[i].lower(),types[i],defaults[types[i]],
                  "NOT NULL" if types[i] in notnull else "") \
                  for i in range(len(header))]
    # Include as many non-TEXT columns in primary key as allowed (16)
    # Additional INFO (like END) may justify duplicates within the standard VCF fields
    query = "CREATE TABLE IF NOT EXISTS %s.%s (%s, PRIMARY KEY(%s))"
    pk = ','.join([h for i,h in enumerate(header) if types[i]!="TEXT"][:16])
    query = query%(io.dbname,self.dname,', '.join(table_def),pk)
    # Create the table
    io.secure_command(query)
    # Populate the table with contents of VCF file
    query_head  = "INSERT IGNORE INTO %s.%s "%(io.dbname,self.dname)
    query_head += "(%s) VALUES "%','.join(['`%s`'%h for h in header])
    query = query_head
    def record2row(record,infos):#,csq_header=None):
      row  = [record.CHROM,record.POS,record.ID]
      row += [record.REF,record.ALT,record.QUAL,record.FILTER]
      # Only retain the most frequent alternate allele (force biallelic)
      row += [record.INFO[f] if f in record.INFO else None for f in infos.keys()]
      # Replace any empty lists with None
      row =  [r if type(r)!=list or len(r)<1 else r[0] for r in row]
      row =  [r if r else None for r in row]
      return row
    # Add the first record to insert rows
    rows   = record2row(record,parser.infos)
    query += " "+queryf
    for i,record in enumerate(parser):
      # Buffered upload - Flush as specified in config
      if not (i+1) % buffer_size:
        try:
          io.secure_command(query,rows)
        except:
          print "Isolating the problem row..."
          rowlen = queryf.count("%")
          for i in range(0,len(rows),rowlen):
            row = rows[i:i+rowlen]
            query = query_head+" "+queryf
            io.secure_command(query,row)
        rows  = []
        query = query_head+" "+queryf
      else:
        query += ", "+queryf
      rows.extend(record2row(record,parser.infos))
    if rows:
      # If rows in buffer at EOF, flush buffer
      try:
        io.secure_command(query,rows)
      except:
        print "Isolating the problem row..."
        rowlen = queryf.count("%")
        for i in range(0,len(rows),rowlen):
          row = rows[i:i+rowlen]
          query = query_head+" "+queryf
          io.secure_command(query,row)
    # Return the number of rows uploaded
    return i

  def load_bedfile(self,fname,io,delim='\t',indexing=None,vepprep=True):
    """ Creates a supplementary table for original BED datafile """
    print "Uploading original datafile to supplementary database."
    with open(fname,'rb') as fin:
      header = fin.readline().strip().split('\t')
      header[0:4] = ["chr","start","end","name"]
      for char in ['.',' ','/','\\','(',')','[',']','-','!','+','=']:
        header = [f.replace(char,'_') for f in header]
      reader = csv.reader(fin,delimiter='\t')
      row1   = reader.next()
      prepend_chr = True if row1[0][:3] != "chr" else False
    # Remove header from bed file
    types = ["VARCHAR(100)","INT","INT","VARCHAR(100)"]
    for i,col in enumerate(row1):
      if i == 3:
        # allele specification, VEP default
        if re.search('\S+\/\S+',col):
          id_type = 'default'
        # BED format with existing rsIDs
        elif col.startswith('rs'):
          id_type = 'id'
        # BED format with HGVS notation
        else:
          id_type = 'hgvs'
      if i > 3:
        try:
          col = float(col)
          types.insert(i,"DOUBLE")
        except ValueError,TypeError:
          if len(col) > 150:
            types.insert(i,"TEXT")
          else:
            types.insert(i,"VARCHAR(250)") # reasonably large varchar
    table_def = ["`%s` %s"%(header[i],types[i]) for i in range(len(header))]
    query = "CREATE TABLE IF NOT EXISTS %s.%s (%s, PRIMARY KEY(chr,start,end,name))"
    query = query%(io.dbname,self.dname,','.join(table_def))
    io.secure_command(query)
    query  = "LOAD DATA LOCAL INFILE '%s' INTO TABLE %s.%s "%(fname,io.dbname,self.dname)
    if delim != '\t':
      query += r"FIELDS TERMINATED BY '%s'"%delim
    else:
      query += r"FIELDS TERMINATED BY '\t'"
    query += "IGNORE 1 LINES"
    io.secure_command(query)
    if prepend_chr:
      query = "UPDATE %s.%s SET chr=CONCAT('chr',chr)"%(io.dbname,self.dname)
      io.secure_command(query)
    # Make any necessary indexing conversions
    if indexing == 'ucsc':
      query = "UPDATE IGNORE %s SET start=start+1, end=end+1 ORDER BY start,end DESC"%self.dname
      io.secure_command(query)
    elif indexing == 'ensembl':
      query = "UPDATE IGNORE %s SET end=end+1 ORDER BY start,end DESC"%self.dname
      io.secure_command(query)
    if vepprep:
    # Retain only rsID if not VEP-default format
      if id_type != 'default':
        if delim == '\t':
          os.system("sed '1d' %s | cut -f4 > %sVEPPREP.txt"%(fname,fname.replace('.','_')))
        else:
          os.system("sed '1d' %s | cut -d'%s' -f4 > %sVEPPREP.txt"%(fname,delim,fname.replace('.','_')))
        fname = "%sVEPPREP.txt"%fname.replace('.','_')
    return fname,id_type

  def load_vep(self,fname,intype='vcf',outfile=None):
    """ Yield VEP rows """
    print "Initializing VEP generator."
    cmd = [f for f in self.vep_cmd]
    cmd[2] = fname
    if intype == 'default':
      # No format specification required
      cmd = cmd[0:3]+cmd[5:]
    else:
      # Provide the correct format argument (vcf;id;hgvs)
      cmd[4] = intype
    print ' '.join(cmd)
    try:
      # Call VEP and capture stdout in realtime
      p = sp.Popen(cmd,stdout=sp.PIPE,bufsize=1)
      # Filter variants without consequence annotations
      p1 = p
      p = sp.Popen(["grep","^#\|CSQ"],stdin=p1.stdout,stdout=sp.PIPE,bufsize=1)
      if self.dname == '1kg3':
        # Pipe output to vcf_fst for Fst calculations
        p2 = p
        p = sp.Popen(["bash","lib/vcf_fst.sh"],stdin=p2.stdout,stdout=sp.PIPE,bufsize=1)
      if outfile:
        fout = gzip.open(outfile,'wb') # Open cache for writing
      for line in iter(p.stdout.readline,b''):
        if outfile: # Write to cache before yielding
          fout.write(line)
        yield line
      if outfile:
        fout.close() # Close the cache
    except KeyboardInterrupt:
      msg = "ERROR (PDBMapData) Keyboard interrupt. Canceling VEP...\n"
      sys.stderr.write(msg)
      raise
    except:
      msg = "ERROR (PDBMapData) Unknown VEP error. Examine stacktrace...\n"
      sys.stderr.write(msg)
      raise
    finally:
      try: p.kill() # kill any running subprocess
      except: pass # ignore any subprocess exceptions

  def _parse_csq(self,csq_header,csq_values):
    """ Creates a dictionary from VEP CSQ desc and each row of values """
    # grab column headers
    res = []
    for row in csq_values:
      csq = {}
      ## Edits to make while parsing fields (inner loop)
      for i,field in enumerate(row.split('|')):
        csq[csq_header[i]] = field
        # Check the consequence and reformat
        if csq_header[i] == "Consequence":
          cons = field.split('&')
          ## We are now allowing Synonymous SNPs to be mapped ##
          # if not any([con in self._parse_csq.nonsyn for con in cons]):
          #   csq = None
          #   break # Ignore this row. Invalid consequence.
          csq['Consequence'] = ';'.join(cons)
        # Set any empty strings to None and continue
        elif csq[csq_header[i]] == '':
          csq[csq_header[i]] = None
          continue
        # Reformat the reference and alternate amino acids
        elif csq_header[i] == "Amino_acids":
          if not field: # No Value
            ref,alt = None,None
          elif len(field.split('/')) < 2: # Only one amino acid recorded
            ref,alt = field,field
          else: # Both specified
            ref,alt = field.split('/')[0:2]
          csq["Ref_AminoAcid"] = ref
          csq["Alt_AminoAcid"] = alt
        # Reformat the reference and alternate codons
        elif csq_header[i] == "Codons":
          if not field: # No value
            ref,alt = None,None
          elif len(field.split('/')) < 2: # Only one codon recorded
            ref,alt = field,field
          else: # Both specified
            ref,alt = field.split('/')[0:2]
          csq["Ref_Codon"] = ref
          csq["Alt_Codon"] = alt
        # Reformat the existing variation
        elif csq_header[i] == "Existing_variation":
          csq["Existing_variation"] = csq["Existing_variation"].replace('&',',')
        # Reformat the domains
        elif csq_header[i] == "DOMAINS":
          csq['DOMAINS'] = csq['DOMAINS'].replace('&',',')
        # Convert canonical flag to boolean
        elif csq_header[i] == "CANONICAL":
          csq["CANONICAL"] = 1 if csq["CANONICAL"] == "YES" else 0
        ##Transform valid 0-indexed positions to 1-indexed positions
        #UPDATE: Despite initial appearances, positions seem to be correct
        elif csq_header[i] == "Protein_position":
          if csq['Protein_position'] and '?' not in csq['Protein_position']:
            csq["Protein_position"] =  int(csq['Protein_position'].split('-')[0])# + 1
          else: csq["Protein_position"] = None
        elif csq_header[i] == "cDNA_position":
          if csq['cDNA_position'] and '?' not in csq['cDNA_position']:
            csq["cDNA_position"] =  int(csq['cDNA_position'].split('-')[0])# + 1
          else: csq["cDNA_position"] = None
        elif csq_header[i] == "CDS_position":
          if csq['CDS_position'] and '?' not in csq['CDS_position']:
            csq["CDS_position"] =  int(csq['CDS_position'].split('-')[0])# + 1
          else: csq["CDS_position"] = None

      if not csq: continue # Invalid conseqeunce. Skip the entry.

      ## Edits to make after parsing fields completely (outer loop)
      # Transcript is not canonical if no value was given
      if "CANONICAL" not in csq: csq["CANONICAL"] = 0
      # Define ref/alt amino acids if not specified
      if "Amino_acids" not in csq or not csq["Amino_acids"]:
        csq["Ref_AminoAcid"] = None
        csq["Alt_AminoAcid"] = None
      # Define ref/alt codons if not specified
      if "Codons" not in csq or not csq["Codons"]:
        csq["Ref_Codon"] = None
        csq["Alt_Codon"] = None
      # For any field not specified, add with value None
      for header in csq_header:
        if header not in csq:
          csq[header] = None
        elif csq[header] == '':
          csq[header] = None
      res.append(csq)
    return res
    # The Variant Effect Predictor will return all consequences for a
    # variant if any one of those consequences passes the filter. This
    # list is used to filter the erroneous consequences of otherwise
    # relevant variants.
  _parse_csq.nonsyn = ['initiator_codon_variant','inframe_insertion',
                      'inframe_deletion','missense_variant','stop_gained',
                      'stop_lost','frameshift_variant','splice_region_variant',
                      'transcript_ablation']

  def _parse_row(self,line):
    """ Parse one VEP output line into a vep row dictionary """
    line = line.lower()
    row  = line.strip().split('\t')
    vep  = {}

    # Safe type casting
    def sfloat(key,d=None):
      if d:
        if key not in d: return None
        else: key = d[key]
      try: return float(key) 
      except: return None
    def sint(key,d=None):
      if d:
        if key not in d: return None
        else: key = d[key]
      try: return int(key) 
      except: return None
    def sbool(key,d=None):
      if d:
        if key not in d: return None
        else: key = d[key]
      if key == 'yes': return True
      if key == 'no': return False
      try: return bool(key) 
      except: return None

    print "Debug: %s"%' '.join(row)

    vep["name"] = row[0]
    loc = row[1].split(':')  # split chr and loc
    vep["chr"] = loc[0]
    loc = loc[1].split('-') # split start and end
    vep["start"] = sint(loc[0])
    vep["end"] = sint(loc[0])+1 if len(loc) < 2 else sint(loc[1])
    vep["ref_allele"] = row[2]
    vep["gene"] = row[3]
    vep["feature"] = row[4]
    vep["feature_type"] = row[5]
    vep["consequence"] = row[6]
    vep["cdna_pos"] = None if not row[7] else sint(row[7])
    vep["cds_pos"] = None if not row[8] else sint(row[8])
    vep["protein_pos"] = None if not row[9] else sint(row[9])
    aas = row[10].split('/')
    vep["ref_amino_acid"] = '-' if len(aas) < 2 else aas[0]
    vep["var_amino_acid"] = '-' if len(aas) < 2 else aas[1]
    codons = row[11].split('/')
    vep["ref_codon"] = '-' if len(codons) < 2 else codons[0]
    vep["var_codon"] = '-' if len(codons) < 2 else codons[1]
    vep["variation"] = row[12]
    # Process any extra columns
    extra = row[13].split(';') if len(row) > 13 else []
    extra = dict([tuple(x.split('=')) for x in extra])
    vep["aa_maf"]  = sfloat("aa_maf",extra)
    vep["afr_maf"] = sfloat("afr_maf",extra)
    vep["amr_maf"] = sfloat("amr_maf",extra)
    vep["asn_maf"] = sfloat("asn_maf",extra)
    vep["ea_maf"]  = sfloat("ea_maf",extra)
    vep["eur_maf"] = sfloat("eur_maf",extra)
    if 'gmaf' in extra: 
      gmaf = extra['gmaf'].split(',')[0]
      ma,maf = gmaf.split(':')
    else: ma,maf = '-',''
    vep["var_allele"] = ma
    vep["gen_maf"] = sfloat(maf) 
    vep["biotype"] = extra.get("biotype",'-')
    vep["canonical"] = sbool("canonical",extra)
    vep["ccds"] = extra.get("ccds",'-')
    vep["clin_sig"] = extra.get("clin_sig",'-')
    vep["distance"] = sint("distance",extra)
    vep["domains"] = extra.get("domains",'-')
    vep["ensp"]  = extra.get("ensp",'-')
    vep["exon"]  = extra.get("exon",'-')
    vep["hgvsc"] = extra.get("hgvsc",'-')
    vep["hgvsp"] = extra.get("hgvsp",'-')
    vep["intron"] = extra.get("intron",'-')
    vep["pubmed"] = extra.get("pubmed",'-')
    vep["polyphen"] = sfloat("polyphen",extra)
    vep["sift"] = sfloat("sift",extra)
    vep["gene_alias"] = extra.get("symbol",'-')
    return(vep)

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
