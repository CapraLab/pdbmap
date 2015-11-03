def pdb_count(df):
  return len(df.drop_duplicates(["structid","biounit"]))

def annotate_unp(fname="all_pi.txt"):
  import csv
  # Annotate the sphere analysis results with UniProt IDs
  with open("/Volumes/sivleyrm/pdbmap/data/pi_analysis/chain2unp.txt",'rb') as fin:
    fin.readline() # burn the header
    reader = csv.reader(fin,delimiter='\t')
    chain2unp = dict([((sid,bio,chain),unp) for sid,bio,chain,unp in reader])
  # Create a new file appending UniProt IDs to the dPi results
  with open("unp_%s"%fname,'wb') as fout:
    writer = csv.writer(fout,delimiter='\t')
    with open(fname,'rb') as fin:
      reader = csv.reader(fin,delimiter='\t')
      writer.writerow(["unp"]+reader.next()) # update the header
      for row in reader:
        sid,bio,chain = row[0],row[1],row[3] # skip the model
        if (sid,bio,chain) in chain2unp:
          writer.writerow([chain2unp[(sid,bio,chain)]]+row)
  # Load and return the UniProt-annotated results
  return "unp_%s"%fname

def load_results(fname="all_pi.txt",unp=False):
  """ Loads results from file and applies basic filters """
  import pandas as pd,numpy as np
  if unp:
    fname = annotate_unp(fname)
  df = pd.read_csv(fname,header=0,sep='\t')
  if "nbrsnps1" in df.columns:
    unit = "spheres" 
  elif "structid" in df.columns:
    unit = "structures"
  else:
    unit = "genes"
  if unit=="spheres":
    print "\nTotal 2+ SNP spheres (%d PDBs) examined:\t\t%d"%(pdb_count(df),len(df))
  # Skip the neighbor SNP check if loading whole-structure results
  if "nbrsnps1" in df.columns and "nbrsnps2" in df.columns:
    # Are there spheres where neighbor SNPs listed as NaN? If so, make sure their SNP count is 0
    df.ix[df["nbrsnps1"].astype(str)=="nan","snpcnt1"] = 0
    df.ix[df["nbrsnps2"].astype(str)=="nan","snpcnt2"] = 0
  # Identify spheres where deltaPi is NaN. Corresponds to cMAFs of 0.5 and 1.0. Unknown cause.
  df = df[~np.isnan(df["pi1"])]
  df = df[~np.isnan(df["pi2"])]
  print "Total %s after filtering NaN:\t\t\t%d"%(unit,len(df))
  # Filter any missed spheres with fewer than 2 SNPs
  df = df[(df["snpcnt1"]>1) | (df["snpcnt2"]>1)]
  print "Total %s after checking 2+ SNPs:\t\t\t%d"%(unit,len(df))
  # Filter any spheres with a value of 0
  df["absdpi"] = abs(df["dpi"])
  df = df[df["absdpi"]>0]
  print "Total %s after removing non-differentiated:\t%d"%(unit,len(df))
  if unit=="spheres":
    print "%d spheres (%d PDBs) passed all criteria"%(len(df),pdb_count(df))
  elif unit=="structures":
    print "%d structures passed all criteria"%len(df)
  else:
    print "%d genes passed all criteria"%len(df)
  return df

def filter_fet(df,p=0.01,adj=True):
  """ Filters spheres with adjusted p-value exceeding threshold """
  if adj:
    tdf = df[df["fet_padj"]<p]
  else:
    tdf = df[df["fetp"]<p]
  print "%d spheres (%d PDBs) with FET p.adj<%g"%(len(tdf),pdb_count(tdf),p)
  return tdf

def collapse(df):
  """ Collapses spheres with identical SNP sets to the minimum dPi """
  import numpy as np
  tdf = df.sort("absdpi",ascending=True).replace(np.nan,-999).groupby(["structid","biounit","nbrsnps1","nbrsnps2"],as_index=False).nth(0).replace(-999,np.nan)
  print "%d non-redundant spheres (%d PDBs)."%(len(tdf),pdb_count(tdf))
  return tdf

def tails(df,p=5,tail='both'):
  """ Returns spheres in the tails of the dPi distribution """
  import numpy as np
  p  = p if p>=1 else p*100
  if tail == 'both':
    pl  = np.percentile(df['dpi'].values,p/2.)
    pu  = np.percentile(df['dpi'].values,100-p/2.)
    print "  Lower %5.1f percentile: % g"%(p/2.,pl)
    print "  Upper %5.1f percentile: % g"%(100-p/2.,pu)
    tdf = df[(df['dpi']<pl) | (df['dpi']>pu)]
  elif tail == 'lower':
    pl  = np.percentile(df['dpi'].values,p)
    print "  Lower %5.1f percentile: % g"%(p,pl)
    tdf = df[(df['dpi']<pl)]
  elif tail == 'upper':
    pu  = np.percentile(df['dpi'].values,100-p)
    print "  Upper %5.1f percentile: % g"%(100-p,pu)
    tdf = df[(df['dpi']>pu)]
  else:
    raise Exception("Invalid tail type. Valid options: lower,uppper,both")
  print "%d spheres (%d PDBs) in the %d%% tail(s) of the dPi distribution"%(len(tdf),pdb_count(tdf),p)
  return tdf

def save_sids(df,fname="sphere_sids.txt"):
  """ Saves unique set of SIDs from dataframe to file """
  import numpy as np
  tdf = df.drop_duplicates(["structid","biounit"])[["structid","biounit"]]
  np.savetxt(fname,tdf,fmt="%s",delimiter="\t")

