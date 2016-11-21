#!/usr/bin/env python2.7

import numpy as np
from scipy.spatial.distance import pdist, squareform
import pandas as pd
from scipy.stats.mstats import zscore
from scipy.stats import norm,percentileofscore
from math import factorial
from collections import OrderedDict
import random
import argparse
import csv
import time

parser = argparse.ArgumentParser()
parser.add_argument('-coord_file', type=str,
                  help='PDB coordinate file on which simulations will be done.')
parser.add_argument('--radii', type=int, default=20, help='Max radius to be simulated as distance threshold.')
parser.add_argument('--n_points', type=int, default=10, help='Max number of points to sample within each cluster.')
parser.add_argument('--noise', default=False, action='store_true')
parser.add_argument("--ppart",type=int,default=1,
                    help="Number of parallel partitions")
parser.add_argument("--ppidx",type=int,default=0,
                    help="Assigned partition")
parser.add_argument("--outdir",type=str,default=".results/K_Simulation/",
                    help="Alternate output path (detail recommended)")
args = parser.parse_args()

K_parameters = ['structid','chain','rsa','radius','size']
K_stats = ['K','K_p','K_z','minT','maxT','meanT','sigT']
cols = ['structid','biounit','model','chain','seqid','icode',
       'x','y','z','rsa','ss','?','pathogenic']
coords = ['x','y','z']
stored_results = []
radii = np.arange(5, args.radii + 1, 1)
n_points = np.arange(3, args.n_points + 1, 1)

def perm(y,N):
  """ Support function for Kest simulations """
  for i in range(N):
    yield np.random.permutation(y)
def Kest(D,y,T=[],P=999):
  """ Ripley's K-Function Estimator for Spatial Cluster Analysis (w/ Positional Constraints) (w/o Edge Correction)
      D: Distance matrix for all possible point pairs (observed and unobserved)
      y: Weight vector for all possible points (un-observed points must have NaN weight)
      T: Distance thresholds
      P: Number of permutations for simulated confidence envelope
      Caveat: Current implementation only handles positive weights"""
  assert(P!=1)
  y = np.array(y,dtype=np.float64) # for NaN compatibility
  weighted = (y==0).sum()+(y==1).sum() != y[~np.isnan(y)].size
  if not weighted: # convert 0/1 to nan/1
    y[y==0]         = np.nan
    y[~np.isnan(y)] = 1.
  Do  = D[~np.isnan(y),:][:,~np.isnan(y)] # Distance of observed points
  yo  = y[~np.isnan(y)]  # Values of observed points
  R   = y.size           # Number of protein residues
  N   = yo.size          # Number of observed points
  if weighted:
    if P:
      Kest.DT = [Do>=t for t in T] # Precompute distance masks
    Y  = np.outer(yo,yo) # NaN distance diagonal handles i=j pairs
    Y /= Y.sum() # normalize by all-pairs product-sum
    K  =  np.array([np.ma.masked_where(dt,Y).sum() for dt in Kest.DT])
  else:
    K = np.array([(Do<t).sum() for t in T],dtype=np.float64) / (N*(N-1))
  if P:
    if weighted:
      # If weighted, shuffle values for observed points
      K_perm = np.array([Kest(Do,yp,T,P=None) for yp in perm(yo,P)])
    else:
      # If unweighted, shuffle positions for all points
      K_perm = np.array([Kest(D,yp,T,P=None) for yp in perm(y,P)])
    # Add the observed K vector to the permutation matrix
    K_perm = np.concatenate(([K],K_perm))
    # Calculate the simulated z-score 
    K_z = zscore(K_perm)[0]
    if all(np.isnan(K_z)):
      # If all weights are equal, set K_z to 0. rather than NaN
      K_z = np.zeros(K_z.shape[0])
    # Calculate one-sided permutation p-value given K directionality
    p = []
    for i,z in enumerate(K_z):
      if z>0:
        p.append(min(1.,2.*(1.-percentileofscore(K_perm[:,i],K[i],'strict')/100.)))
      else:
        p.append(min(1.,2.*(1.-percentileofscore(-K_perm[:,i],-K[i],'strict')/100.)))
    K_p = np.array(p)
    # Min/Max/Mean T with nominally significant p-values (p<0.01)
    if any(K_p < 0.01):
      minT,maxT,meanT = np.min(T[K_p<0.01]),np.max(T[K_p<0.01]),np.mean(T[K_p<0.01])
    else:
      minT = maxT = meanT = np.nan
    # T where p-value is most significant (min(p))
    sigT = T[np.argmax(np.abs(K_z))]
    return K,K_p,K_z,minT,maxT,meanT,sigT
  else:
    return K
Kest.DT = []

# Marks 10% of residues in the coordinate dataframe to be pathogenic. 
def simulate_noise(df):
  tdf = df.copy()
  size = int(len(df) * 0.01)
  size = size if size != 0 else 1
  sdf = df.sample(n=size)
  sdf['pathogenic'] = True
  tdf.update(sdf)
  return tdf

# Builds a list of appropriate T thresholds for Kest to test.
# Returns a range of T's from the minimum intervariant distance in a cluster
# to the maximum intervariant distance in a cluster if half of the maximum
# intervariant distance is less than the minimum intervariant distance. Otherwise,
# return half of the max intervariant distance. 
def calc_intervariant_distances(cluster_df, noise_df):
  cdm = squareform(pdist(cluster_df.values))
  cdm = cdm[np.nonzero(cdm)]
  ndm = squareform(pdist(noise_df.values))
  ndm = ndm[np.nonzero(ndm)]
  mdf = pd.concat([cluster_df, noise_df])
  mdm = squareform(pdist(mdf.values))
  mdm = mdm[np.nonzero(mdm)]
  c_median_dist = np.median(cdm)
  n_median_dist = np.median(ndm)
  m_median_dist = np.median(mdm)
  minT = cdm.min()
  maxT = 0.5*cdm.max() if 0.5*cdm.max() > minT else cdm.max()
  print 'Median inter-residue distance between clustered residues: %s'%(c_median_dist)
  print 'Median inter-residue distance between noise residues: %s'%(n_median_dist)
  print 'Median inter-residue distance between all pathogenic residues: %s\n'%(m_median_dist) 
  return np.ceil(np.arange(minT, maxT, 1))

# Simulates n_points variants within a radius rad around the point in space specified
# as the index of a row of a coordinates dataframe of residues in a protein structure.
def simulate_variants(df, centroid_idx, n_points, rad):
  tdf = df.copy()
  # Obtain the 3D coordinates of the centroid.
  centroid = tdf.iloc[centroid_idx].values[:3]
  # Compute the distance of each residue in the structure to the centroid.
  dist_to_centroid = [np.linalg.norm(p - centroid) for p in tdf[coords].values]
  noise_df = df[df['pathogenic'] == True]
  noise_idx = list(noise_df.index)
  # Label candidate residues to be marked as pathogenic. To be a candidate, the residue must be within 
  # rad distance of the centroid and must not be already labeled as noise.
  weights = [1 if (dist <= rad and i not in noise_idx) else 0 for i, dist in enumerate(dist_to_centroid)]
  # If at least n_points candidates exist, mark n_points residues as pathogenic. 
  if np.count_nonzero(weights) < n_points:
    msg = 'There are {} < {} residues inside a cluster of radius {}\n'
    raise Exception(msg.format(np.count_nonzero(weights), n_points, rad))
  sdf = tdf.sample(n=n_points, weights=weights)
  sdf['pathogenic'] = True
  # Update dataframe with new, simulated pathogenic residues.
  tdf.update(sdf)
  # Weight vector of all noise and simulated pathogenic residues in the structure.
  K_weights = [1 if marker == True else 0 for marker in tdf['pathogenic'].values]
  T_range = calc_intervariant_distances(sdf[coords], noise_df[coords])
  return tdf, K_weights, T_range

def format_Kest(results):
  res = []
  for element in results:
    try:
      res.append(element.tolist())
    except:
      res.append(element)
  return res

def write_results(results):
  with open(args.outdir + 'simulation_patch_noise.tsv', 'ab') as fout:
    writer = csv.writer(fout, delimiter='\t')
    writer.writerow(parameters + K_results)
# <--------------------------------------------------------------------------------------> # 

structs = []

with open(args.coord_file, 'rb') as structfile:
  structs = [line.rstrip() for line in structfile]

if args.ppart > 1:
  structs = [s for i,s in enumerate(structs) if i % args.ppart == args.ppidx]
  print "Partition %d contains %d structures."%(args.ppidx,len(structs)); print
  time.sleep(args.ppidx%50)

# Parse out relevant fields from coordinate file
for struct in structs:
  with open(struct, 'rb') as infile:
    struct_df = pd.read_csv(infile, sep='\t', header=None, names=cols)
    struct_df = struct_df[['structid','chain','x','y','z','rsa','ss']]
    coord_df = struct_df[['x','y','z']]
    coord_df['pathogenic'] = False

  N = len(coord_df)
  sample_points = np.arange(0, int(0.05 * N), 1)

  if args.noise:
    coord_df = simulate_noise(coord_df)

  for count in sample_points: # Simulate clusters around 5% of the structure's residues
    center_idx = np.random.randint(0, N)
    struct = struct_df.ix[center_idx, ['structid','chain','ss','rsa']]
    msg = "Kest for structure {}-{}, with sample residue {} secondary structure {} "
    msg += "and relative solvent accessibility {}\n"
    print msg.format(struct[0], struct[1], center_idx, struct[2], struct[3])
    for rad in radii: # Simulate len(radii) clusters with radius rad
      for n in n_points: # Simulate i...n (including sampled point) sized clusters within a radius rad around each center.
        print "Simulating radii %s and number of residues per cluster %s\n"%(rad, n)
        parameters = [struct[0], struct[1], struct[3], rad, n]
        print 'Parameters: '
        print parameters; print
        try:
          K_df, K_weights, T_range = simulate_variants(coord_df, center_idx, n, rad)
        except Exception as e:
          print e
          continue
        K_matrix = squareform(pdist(K_df[coords].values))
        K_results = Kest(K_matrix, K_weights, T=T_range)
        K_results = format_Kest(K_results)
        write_results(parameters + K_results)
