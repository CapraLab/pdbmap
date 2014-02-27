#!/usr/bin/env python27
#
# Project        : PDBmap
# Filename       : PDBMapIntersect.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-27
# Description    : Service class designed to intersect structural and 
#                : genomic, and sequence datasets. Does not alter original 
#                : datasets. Calculates intersection with intersectBed and 
#                : stores the crossreference in the PDBMap.IntersectGenome 
#                : database.
#=============================================================================#

# See main check for cmd line parsing
from lib import PDBMapIO

class PDBMapIntersect():
  def __init__(self):
    pass

  @classmethod
  def intersect(slabel,dlabel,dtype='Genomic'):
    # Intersects a structure set with a dataset
    # dtype = ['Genomic','Protein','Structural']
    structures = PDBMapIO.PDBMapIO.get_structures(slabel)
    if dtype == 'Genomic':
      data = PDBMapIO.PDBMapIO.get_genomic_data(dlabel)


# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
