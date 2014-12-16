#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : PDBMapStats.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-27
# Description    : Calculates LD and Fst statistics for variant pairs
#                : intersected with PDBMap. Requires population information.
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv

class PDBMapStats():

  dens_mat = None # Class variable for show_desnity
  
  def __init__(self,io):
    """ Initialize the PDBMapStats object. """
    self.io = io

  

# Copied from biolearn
def multidigit_rand(digits):
  import random
  randlist = [random.randint(1,9) for i in xrange(digits)]
  multidigit_rand = int(''.join([str(x) for x in randlist]))
  return multidigit_rand