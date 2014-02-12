#!/usr/bin/env python27
#
# Project        : PDBMap
# Filename       : PDBMapAlignment.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-12
# Description    : Defines the PDBMapTranscript class for calculating and
#                : representing an alignment between a PDBMapStructure and
#                : a PDBMapTranscript.
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv

class PDBMapAlignment():

  def __init__(self):
	 pass

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)