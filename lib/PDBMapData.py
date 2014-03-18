#!/usr/bin/env python27
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
import sys,os,csv
import subprocess as sp
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import vcf

class PDBMapData():
  
  def __init__(self,vep="variant_effect_predictor.pl",plink="plink",j=1):
    self.vep = vep
    if not os.path.exists(vep):
      msg = "ERROR: (PDBMapData) VEP location invalid: %s"%vep
      raise Exception(msg)
    # Check for a dbconn file
    registry = "%s/dbconn.conf"%os.path.dirname(vep)
    if not os.path.exists(registry):
      msg = "WARNING: (PDBMapData) Not registry specified. Using Ensembl.\n"
      sys.stderr.write(msg)
      registry = None
    # Construct the VEP command
    self.vep_cmd = [self.vep,'-i','','--database']
    if registry:
      # Use a local Ensembl database
      self.vep_cmd.extend(['--registry',registry])
    # Disable the progress bars
    self.vep_cmd.extend(['--no_progress'])
    # Specify the species for faster db queries
    self.vep_cmd.extend(['--species','homo_sapiens'])
    # Increase buffer size to improve runtime (default 5,000)
    self.vep_cmd.extend(['--buffer_size','300000']) # ~15GB
    # Annotate with functional info/prediction
    self.vep_cmd.extend(['--sift','s','--polyphen','s','--regulatory'])
    # Annotate with variant, gene, protein, and domain identifiers
    self.vep_cmd.extend(['--check_existing','--symbol','--protein','--domains'])
    # Annotate transcript with canonical bool and biotype
    self.vep_cmd.extend(['--canonical','--biotype'])
    # Retain only nonsynonymous coding variants
    self.vep_cmd.extend(['--vcf','--filter','coding_change'])
    # Send to stdout and don't generate a summary stats file
    self.vep_cmd.extend(['--no_stats','-o','stdout'])
    if j > 1:
      # Fork to decrease runtime
      self.vep_cmd.extend(['--fork',str(j)])

  def load_vcf(self,fname):
    """ Pipe VCF file through VEP and yield VEP rows """
    print "Initializing VCF generator."
    # Parse the VCF output from VEP
    parser = vcf.Reader(self.load_vep(fname))
    # Determine Info headers
    info_headers = parser.infos.keys()
    # Determine Consequence headers
    csq_headers  = parser.infos['CSQ'].desc.split(': ')[-1].split('|')
    for record in parser:
      print " # Processing %s #"%record.ID
      # Adjust some specific fields
      if "END" not in record.INFO:
        record.INFO["END"] = int(record.POS) + 1
      for header in info_headers:
        if header not in record.INFO:
          record.INFO[header] = None
      record.INFO['AC'] = record.INFO['AC'][0]
      if record.INFO['SNPSOURCE']:
        record.INFO['SNPSOURCE'] = ';'.join(record.INFO['SNPSOURCE'])
      if not record.INFO['AMR_AF']: record.INFO['AMR_AF'] = 0.0
      if not record.INFO['ASN_AF']: record.INFO['ASN_AF'] = 0.0
      if not record.INFO['AFR_AF']: record.INFO['AFR_AF'] = 0.0
      if not record.INFO['EUR_AF']: record.INFO['EUR_AF'] = 0.0
      if record.INFO['AA']:
        record.INFO['AA'] = record.INFO['AA'].upper()
      # Add attribute fields to INFO
      record.INFO["ID"] = record.ID
      record.INFO["CHROM"] = "chr%s"%record.CHROM
      record.INFO["START"] = record.POS
      record.INFO["REF"] = record.REF
      record.INFO["ALT"] = record.ALT[0] # only the most common
      record.INFO["QUAL"] = record.QUAL
      record.INFO["FILTER"] = str(record.FILTER)
      # Extract the consequences
      record.CSQ = self._parse_csq(csq_headers,record.INFO['CSQ'])
      # Add some "consequence" info to the variant info
      record.INFO["EXISTING"] = record.CSQ[0]['Existing_variation']
      record.INFO["GENE"] = record.CSQ[0]['Gene']
      if record.CSQ[0]['SYMBOL_SOURCE'] == 'HGNC':
        record.INFO["HGNC"] = record.CSQ[0]["SYMBOL"]
      else: record.INFO["HGNC"] = None
      # Add some variant info to the consequences
      for csq in record.CSQ:
        csq["ID"] = record.INFO["ID"]
        csq["CHROM"] = record.INFO["CHROM"]
        csq["START"] = record.INFO["START"]
        csq["END"] = record.INFO["END"]
      yield record

  def load_pedmap(self,fname):
    """ Convert PED/MAP to VCF, pipe VCF through VEP and load VEP output """
    print "load_pedmap not implemented"

  def load_bed(self,fname,vep=False):
    """ Load data from BED """
    print "Initializing BED generator."
    if vep:
      for record in vcf.Reader(self.load_vep(fname)):
        yield record
    else:
      print "load_bed not implemented without VEP"

  def load_vep(self,fname):
    """ Yield VEP rows """
    print "Initializing VEP generator."
    cmd = [f for f in self.vep_cmd]
    cmd[2] = fname
    try:
      # Call VEP and capture stdout
      p = sp.Popen(cmd,stdout=sp.PIPE)
      for line in iter(p.stdout.readline,b''):
        yield line
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
          if not any([con in self._parse_csq.nonsyn for con in cons]):
            csq = None
            break # Ignore this row. Invalid consequence.
          csq['Consequence'] = ';'.join(cons)
        # Set any empty strings to None and continue
        elif csq[csq_header[i]] == '':
          csq[csq_header[i]] = None
          continue
        # Reformat the reference and alternate amino acids
        elif csq_header[i] == "Amino_acids":
          if not field: # No Value
            ref,alt = None,None
          elif len(field) < 2: # Only one value
            ref,alt = field,field
          else: # Both specified
            ref,alt = field.split('/')[0:2]
          csq["Ref_AminoAcid"] = ref
          csq["Alt_AminoAcid"] = alt
        # Reformat the reference and alternate codons
        elif csq_header[i] == "Codons":
          if not field: # No value
            ref,alt = None,None
          elif len(field) < 2: # Only one value
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
