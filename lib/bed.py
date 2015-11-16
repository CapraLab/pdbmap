# Project        : PDBMap
# Filename       : bed.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-04-27
# Description    : BED file parser designed to emulate PyVCF.
#=============================================================================#

import sys,io,csv,gzip

class Reader:
  def __init__(self,fsock=None,filename=None,compressed=False,indexing='ucsc',
                prepend_chr=False,delimiter='\t',header=True):
    self.prepend_chr = prepend_chr
    self.indexing    = indexing
    if fsock and filename:
      raise Exception("Cannot process both fsock and filename")
    elif fsock:
      if compressed:
        raise Exception("Cannot read compressed file stream")
      else:
        fin = io.open(fsock,'rb')
    elif filename:
      if compressed:
        fin = gzip.open(filename,'rb')
      else:
        fin = open(filename,'rb')
    self.reader = csv.reader(fin,delimiter=delimiter)
    self.peek = self.reader.next()
    if header:
        self.header = self.peek
        self.NCOLS  = len(self.peek)
        self.peek    = None # drop the header
    else:
        self.header = ["CHROM","START","END","ID"]
        self.NCOLS  = len(self.peek)
        if self.NCOLS != len(self.header):
          sys.stderr.write("BED file contains additional columns but no header.\n")
          self.header += ["V%d"%i for i in range(self.NCOLS-len(self.header))]
    self.infos = dict((f,[]) for f in self.header)
    self.infos["CSQ"] = Consequence()

  def __iter__(self):
    if self.peek:
      row = self.peek
      self.peek = None
      if self.prepend_chr and not row[0].startswith('chr'):
        row[0] = 'chr'+row[0]
      yield Record(self.header,row,self.indexing)
    for row in self.reader:   
      if self.prepend_chr and not row[0].startswith('chr'):
        row[0] = 'chr'+row[0]
      yield Record(self.header,row,self.indexing)

class Record:
  def __init__(self,info,values,indexing):
    self.CHROM  = values[0]
    self.START  = int(values[1]) + int((indexing == 'ucsc'))
    self.END    = int(values[2]) + int((indexing == 'ucsc'))
    self.POS    = int(values[1]) + int((indexing == 'ucsc'))
    self.ID     = values[3]
    self.REF    = values[info.index("REF")]    if "REF"    in info else None
    self.ALT    = values[info.index("ALT")]    if "ALT"    in info else [None]
    self.QUAL   = values[info.index("QUAL")]   if "QUAL"   in info else None
    self.FILTER = values[info.index("FILTER")] if "FILTER" in info else None
    self.FORMAT = values[info.index("FORMAT")] if "FORMAT" in info else None
    self.INFO   = dict((f,values[i]) for i,f in \
                  enumerate(info) if f not in dir(self))
    if "CSQ" not in self.INFO:
      self.INFO["CSQ"] = ""
    self.samples = []

class Consequence:
  def __init__(self):
    self.desc = ""