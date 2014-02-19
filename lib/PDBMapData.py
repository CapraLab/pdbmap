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
#                : upload of annotation files to the PDBMap database.
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv
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

  def load_vcf(self,dname,fname):
    """ Pipe VCF file through VEP and load VEP output """
    print "load_vcf"

  def load_pedmap(self,dname,fname):
    """ Convert PED/MAP to VCF, pipe VCF through VEP and load VEP output """
    print "load_pedmap"

  def load_bed(self,dname,fname):
    """ Load data from BED """
    print "load_bed"

  def load_vep(self,fname,run=False):
    """ Load data from VEP, run VEP on input if specified """
    if run:
      cmd = self.vep_cmd % fname
      # Make the call to VEP
      # Capture the output line-by-line
      # Identify the column headers
      # Split the final column into individual columns
    # How to return the data to pdbmap for upload without storing?
    print "load_vep"

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
