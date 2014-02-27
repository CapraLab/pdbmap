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

class PDBMapData():
  
  def __init__(self,vep="variant_effect_predictor.pl",plink="plink"):
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
    self.vep_cmd = "%s -i %s --database "%(self.vep,"%s")
    if registry:
      self.vep_cmd += "--registry %s "%registry
    self.vep_cmd += "--everything -o stdout"

  def load_vcf(self,fname):
    """ Pipe VCF file through VEP and yield VEP rows """
    for row in self.load_vep(fname,run=True):
      yield row

  def load_pedmap(self,fname):
    """ Convert PED/MAP to VCF, pipe VCF through VEP and load VEP output """
    print "load_pedmap not implemented"

  def load_bed(self,fname,vep=False):
    """ Load data from BED """
    if vep:
      for row in self.load_vep(fname,run=True):
        yield row
    else:
      print "load_bed not implemented without VEP"

  def load_vep(self,fname,run=False):
    """ Yield VEP rows """
    # Define the VEP column fields
    cols = ["uploaded_variation","location","ref_allele","gene","feature","feature_type"]
    cols.extend(["consequence","cdna_pos","cds_pos","protein_pos",])
    cols.extend(["amino_acids","codons","variation"])
    cols.extend(["aa_maf","afr_maf","amr_maf","asn_maf","ea_maf","eur_maf","biotype",])
    cols.extend(["canonical","ccds","clin_sig","distance","domains"])
    cols.extend(["ensp","exon","gen_maf","hgvsc"])
    cols.extend(["hgvsp","intron","pubmed","polyphen","sift","gene_alias"])
    # If another file format needs to passed through VEP
    if run:
      cmd = self.vep_cmd % fname
      # Call VEP and capture stdout
      p = sp.Popen(cmd,shell=True,stdout=sp.PIPE)
      generator = p.stdout.readline
    else:
      # Open the file for reading
      generator = open(fname,'rb')
    # Capture the output line-by-line
    for line in iter(generator):
      if line[0] == "#": continue
      yield _parse_row(line)
    print "load_vep"

  def _parse_row(line):
    """ Parse one VEP output line into a vep row dictionary """
    line = line.lower()
    row  = line.strip().split('\t')
    vep  = {}
    vep["name"] = row[0]
    loc = row[1].split(':')
    vep["chr"] = loc[0]
    vep["start"] = int(loc[1])
    vep["end"] = int(loc[1])+1
    vep["ref_allele"] = row[2]
    vep["gene"] = row[3]
    vep["feature"] = row[4]
    vep["feature_type"] = row[5]
    vep["consequence"] = row[6]
    vep["cdna_pos"] = int(row[7])
    vep["protein_pos"] = int(row[8])
    amino_acids = row[9].split('/')
    vep["ref_amino_acid"] = amino_acids[0]
    vep["var_amino_acid"] = amino_acids[1]
    codons = row[10].split('/')
    vep["ref_codon"] = codons[0]
    vep["var_codon"] = codons[1]
    vep["variation"] = row[11]
    extra = row[12].split(';')
    extra = dict([tuple(x.split('=')) for x in extra])
    vep["aa_maf"] = float(extra["aa_maf"])
    vep["afr_maf"] = float(extra["afr_maf"])
    vep["amr_maf"] = float(extra["amr_maf"])
    vep["asn_maf"] = float(extra["asn_maf"])
    vep["ea_maf"] = float(extra["ea_maf"])
    vep["eur_maf"] = float(extra["eur_maf"])
    vep["gen_maf"] = float(extra["gmaf"])
    vep["biotype"] = extra["biotype"]
    vep["canonical"] = bool(extra["canonical"])
    vep["ccds"] = extra["ccds"]
    vep["clinical_sig"] = extra["clin_sig"]
    vep["dist2trans"] = extra["distance"]
    vep["domains"] = extra["domains"]
    vep["ensp"] = extra["ensp"]
    vep["exon"] = extra["exon"]
    vep["hgvsc"] = extra["hgvsc"]
    vep["hgvsp"] = extra["hgvsp"]
    vep["intron"] = bool(extra["intron"])
    vep["pubmed"] = extra["pubmed"]
    vep["polyphen"] = extra["polyphen"]
    vep["sift"] = extra["sift"]
    vep["gene_alias"] = extra["SYMBOL"]
    return(vep)

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
