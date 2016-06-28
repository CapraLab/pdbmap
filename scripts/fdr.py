#!/usr/bin/env python2.7
# Mike Sivley
# 11-16-2015

import csv,sys,argparse
desc   = "Appends FDR-corrected p-value column to results file"
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("--infile",nargs='?',type=argparse.FileType('rb'),
                    default=sys.stdin,help="Results containing p-value column")
parser.add_argument('--outfile', nargs='?',type=argparse.FileType('wb'),
                    default=sys.stdout)
parser.add_argument("pcol",type=int,
                    help="Column containing p-values")
parser.add_argument("alpha",type=float,default=0.05,
                    help="Alpha threshold for correction")
parser.add_argument("--sep",nargs='?',type=str,default="\t",
                    help="Results delimiter")
parser.add_argument("--header",nargs='?',type=bool,default=False,
                    help="Header flag (default False)")
args = parser.parse_args()
args.sep = args.sep.decode("string-escape") 

from statsmodels.sandbox.stats.multicomp import multipletests as fdr
reader  = csv.reader(args.infile,delimiter=args.sep)
if args.header: 
  header = reader.readline()
results = [row for row in reader]
pvals   = [float(row[args.pcol]) for row in results]

_,padj,_,_ = fdr(pvals,alpha=args.alpha,method="fdr_bh")
writer = csv.writer(args.outfile,delimiter=args.sep)
if args.header:
  header.append("adj")
  writer.writerow(header)
for i,row in enumerate(results):
  row.append(padj[i])
  writer.writerow(row)
      
