#!/usr/bin/env python2.7

# Load libraries
import glob,gzip,os,sys
import numpy as np
from numpy.random import randn
from itertools import combinations
import random
import pandas as pd
from scipy import stats
from scipy.spatial import KDTree
from sklearn.metrics import confusion_matrix
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
pd.set_option('display.max_columns', 500)

# Ensure that we are in the main pdbmap directory and not a subdirectory or elsewhere
os.chdir('../../pdbmap')

# Load header information
header_fname  = 'scripts/sphere_header.txt'
with open(header_fname,'rb') as hfin:
    header = [h.strip() for h in hfin.readlines()]
# Load the results file into a pandas data frame
fnames = glob.glob('/dors/capra_lab/sivleyrm/pdbmap/results/sliding_sphere_10/split/obs/bystruct/*.txt.gz')
rlist  = [[r.split('-')[0][-4:],r.split('-')[1].split('.')[0]] for r in fnames]

# Read all structures separately
dfs = [(i,pd.read_table(gzip.open(fname,'rb'),delimiter='\t',skiprows=1,names=header)) for i,fname in enumerate(fnames)]

# Add biological assembly identifiers
for i,df in dfs:
    df['structid'],df['biounit'] = [rlist[i][0] for j in range(len(df))],[rlist[i][1] for j in range(len(df))]

# Filter out problem structures (any structure where an entire columns is NULL)
dfs = [df for i,df in dfs if not np.any(df.isnull().sum()==len(df))]

# Merge all dataframes into a single dataframe
df = pd.concat(dfs)
del dfs

# Drop duplicate residues
df.drop_duplicates(['structid','biounit','x','y','z'],inplace=True)

# Observed values
info_cols      = header[:21]       # Residue information
count_cols     = header[20:36]     # SNP counts 
cdaf_cols      = header[36:41]     # Cumulative derived allele frequencies
ddaf_cols      = header[41:69:3]   # Mean directional deltaDAF
ddaf_cols     += header[42:70:3]   # 
absddaf_cols   = header[43:71:3]   # Mean absolute deltaDAF
fst_cols       = header[71:82]     # Fst
w_count_cols   = header[82:98]     # Windowed SNP counts
w_cdaf_cols    = header[98:103]    # Windowed Cumulative derived allele frequencies
w_ddaf_cols    = header[103:131:3] # Windowed mean directional deltaDAF
w_ddaf_cols   += header[104:132:3] # 
w_absddaf_cols = header[105:133:3] # Windowed mean absolute deltaDAF
w_fst_cols     = header[133:144]   # Windowed Fst

# Observed p-values
def pval_cols(cols):
    return ["%s_pval"%col for col in cols]
count_pval_cols     = pval_cols(count_cols)
cdaf_pval_cols      = pval_cols(cdaf_cols)
ddaf_pval_cols      = pval_cols(ddaf_cols)
absddaf_pval_cols   = pval_cols(absddaf_cols)
fst_pval_cols       = pval_cols(fst_cols)
w_count_pval_cols   = pval_cols(w_count_cols)
w_cdaf_pval_cols    = pval_cols(w_cdaf_cols)
w_ddaf_pval_cols    = pval_cols(w_ddaf_cols)
w_absddaf_pval_cols = pval_cols(w_absddaf_cols)
w_fst_pval_cols     = pval_cols(w_fst_cols)

# Nested column sets
value_cols  = [count_cols,cdaf_cols,ddaf_cols,absddaf_cols,fst_cols,w_ddaf_cols,w_absddaf_cols,w_fst_cols]
pvalue_cols = [count_pval_cols,cdaf_pval_cols,ddaf_pval_cols,absddaf_pval_cols,fst_pval_cols,w_ddaf_pval_cols,w_absddaf_pval_cols,w_fst_pval_cols]

# Shift negative Fst values to 0
for col in fst_cols + w_fst_cols:
    df.ix[df[col]<0,col] = 0

# Correct poorly implemented >1 SNP per sphere restriction
# If there are fewer than 2 shared SNPs in the sphere, all measurements
# should be set to 0, and all p-values to 1.0
for i,scount in enumerate(count_cols[-10:]): # Shared SNP counts (spheres)
    df.ix[df[scount]<2,[c[i] for c in value_cols[2:5]]]  = 0.
    df.ix[df[scount]<2,[c[i] for c in pvalue_cols[2:5]]] = 1.

for i,wcount in enumerate(w_count_cols[-10:]): # Shared SNP counts (windows)
    df.ix[df[wcount]<2,[c[i] for c in value_cols[5:]]]  = 0.
    df.ix[df[wcount]<2,[c[i] for c in pvalue_cols[5:]]] = 1.

# Identifies overlapping spheres forming an ROI
# and assigns a biounit-unique RID to the set
def collapse_roi(df,structid,biounit,pmetric,metric):
    """ Significant spheres must have a p-value < 0.01 and an observed value >0.05 """
    dft = df.ix[(df[pmetric]<0.01) & (df[metric]>0.05) & (df.structid==structid) & (df.biounit==biounit),
                ["structid","biounit","model","chain","seqid","icode","x","y","z"]]
    if dft.empty: return dft # biounit contains no significant spheres
    dft = dft.copy().reset_index() # copy and reset index for significant-only residues
    crd = [list(c) for c in np.array(dft[["x","y","z"]],dtype=np.float64)]
    kdt = KDTree(crd)
    cdt = [c for c in crd] # candidate set
    dft["%s_roi"%pmetric] = np.nan    # initialize column for RIDs
    # Define the nearest-neighbor walk
    def nnsearch(c,kdt,crd,frnt):
        nbrs  = [n for n in kdt.query_ball_point(c,15) if crd[n] in frnt]
        for n in nbrs:
            frnt.remove(crd[n]) # remove from candidate set
        return [crd.index(c)]+nbrs+[g for n in nbrs for g in nnsearch(crd[n],kdt,crd,frnt)] 
    # Run the nearest-neighbor walk and assign ROIs
    i = 0
    while len(cdt)>0:
        for n in nnsearch(cdt.pop(0),kdt,crd,cdt):
            dft.ix[n,"%s_roi"%pmetric] = i
        i += 1
    # Return the ROI-annotated copy of df
    return dft

for pcolset in pvalue_cols[3:5]: # abs(dDAF) and Fst
    for p,pcol in enumerate(pcolset):
        outdir = 'results/roi' 
        os.system('mkdir -p %s'%outdir)
        if p!=int(sys.argv[1]): continue
        col = '_'.join(pcol.split('_')[:-1])
        fout = open('%s/%s.txt'%(outdir,col),'wb')
        fout.write("structid\tbiounit\troi\tnumspheres\n")
        dfroi = []
        for i,(sid,biounit) in enumerate(rlist):
            #print "\rROI determination for %20s %3d%% complete (%s.%s)"%(pcol,int(float(i)/len(rlist)*100),sid,biounit),
            dfroi.append(collapse_roi(df,sid,biounit,pcol,col))
        #print "\rROI determination for %20s 100% complete"%(pcol)
        dfroi  = pd.concat(dfroi)
        if not dfroi.empty:
            dfroig = dfroi.groupby(["structid","biounit","%s_roi"%pcol])
            for _,group in sorted(dfroig,key=lambda x: len(x[1]),reverse=True):
                fout.write("%s\t%s\t%3d\t%3d\n"%(group["structid"].iloc[0],group["biounit"].iloc[0],group["%s_roi"%pcol].iloc[0],len(group)))
        fout.close()
