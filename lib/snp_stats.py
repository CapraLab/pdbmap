#!/usr/bin/python27

# Command Line Arguments:
# calc_fst.r pop-id-map pedmap krange snp-pdbid-loc

import sys,os,csv,subprocess,math
from itertools import combinations,product

def main():

  ## Check cmd line arguments
  if len(sys.argv) < 2:
    print "usage: \n\tcalc_fst.py pop-id-map pedmap krange snp-pdbid-loc"
    print "\tpop-id-map is a tab delimited file with [population, sampleid] columns"
    print "\tsnp-pdbid-loc is a tab delimited file with [gene, snp_name, pdbid, chain, x, y , z] columns"
    sys.exit(1)

  ## Read cmd line arguments
  pop_id_map = sys.argv[1]
  pedmap     = sys.argv[2]
  kbounds = (0,0)
  if len(sys.argv) > 2:
    kbounds    = sys.argv[3].split(':')
  krange = range(int(kbounds[0]),int(kbounds[1])+1)
  snp_groups = None
  if len(sys.argv) > 3:
    snp_groups = sys.argv[4]

  ## Preprocess the pedmap for hierfstat if necessary
  if not os.path.isfile("%s-tab1234.ped"%pedmap):
    preprocess_pedmap(pedmap)

  ## Read the populations
  pop_id_dict = read_pop(pop_id_map)

  ## Read all SNPs, and their chr,bp
  snp_map = read_map(pedmap)

  ## Remove SNPs with MAF == 0.0
  snp_map = qc_freq(snp_map,pedmap)
  print "\n# %d SNPs found in PEDMAP."%len(snp_map)

  ## Read SNP groups if specified
  if snp_groups:
    group_snp_dict,snp_loc,snp_list = read_groups(snp_groups,snp_map)
    print "# Evaluating %d overlapping SNPs."%len(snp_list)
  # Only group-associated SNPs present in the snp_map are recorded
  # As such, snp_list holds the union of the pedmap and pdbmap SNPs

  ## Compute pairwise distance and LD (R^2) for all SNPs
  print "\n# Calculating LD statistics and 1D distances..."
  snp_ld,snp_bp_dist = calc_ld_dist(pedmap,snp_list)

  ## Compute pairwise structural distance for all SNPs
  print "# Calculating 3D distances..."
  snp_3d_dist = calc_3d_dist(snp_loc)

  ## Compute the global and individual Fst scores, initialize snp_fstat
  print "# Calculating global Fst estimates..."
  snp_fstat = calc_global_fst(snp_list,pop_id_map,pedmap)
  if 1 in krange: krange.remove(1) # k=1 implicitly calculated

  ## Compute Fst for k-wise tuples (k > 1), update snp_fstat
  print "# Calculating k-wise Fst estimates..."
  snp_fstat = calc_kwise_fst(snp_fstat,pop_id_map,pedmap,krange,snp_groups,group_snp_dict)
  # snp_fstat = calc_kwise_fst(snp_fstat,pop_id_map,pedmap,krange,snp_list=snp_list)

  ## Aggregate all the measurements and statistics
  print "# Aggregating results..."
  res = aggregate_stats(snp_bp_dist,snp_3d_dist,snp_ld,snp_fstat)

  ## Save all results
  print "# Saving results..."
  save_results(os.path.basename(pedmap),res)

def save_results(pedmap,res):
  # create results directory
  resdir = 'results/%s'%pedmap
  if not os.path.exists(resdir):
    os.makedirs(resdir)
  # copy the scripts used
  os.system("cp -f lib/snp_stats.py %s/"%resdir)
  os.system("cp -f `which calc_fst.r` %s/"%resdir)
  # record the original command
  with open("%s/command.txt"%resdir,'wb') as fout:
    fout.write(" ".join(sys.argv))
  # save the results
  with open("%s/%s.stats"%(resdir,pedmap),'wb') as fout:
    writer = csv.writer(fout,delimiter='\t')
    for row in res:
      writer.writerow(row)

def aggregate_stats(snp_1d_dist,snp_3d_dist,snp_ld,snp_fst):
  import pdb
  # all snps with 3d structural information and a partner
  res = []
  snps = snp_3d_dist.keys()
  snps.sort() # only sort just before printing results
  print "rsID\tFst(ind)\tn(1d)\tn(3d)\tFst(n(1d))\tFst(n(3d))\tLD(n(1d))\tLD(n(3d))"
  try:
    for snp in snps:
      ind_fst = snp_fst[snp][1]
      # get the snp name of the nearest snp in 1d and 3d
      n1d = min(snp_1d_dist[snp],key=snp_1d_dist[snp].get)
      # for all partner snps B, find shortest distance b/w A and B
      # and select the A->B pair with the shortest distance overall
      n3d = min([(min(snp_3d_dist[snp][snpb][0]),snpb) for snpb in snp_3d_dist[snp]])[1]
      # Lookup Fst for nearest partners
      fstn1d = snp_fst[snp][2][(n1d,)]
      fstn3d = float("NaN") if n3d == "NA" else snp_fst[snp][2][(n3d,)]
      # Lookup LD for nearest partners
      ldn1d  = snp_ld[snp][n1d]
      ldn3d  = float("NaN") if n3d == "NA" else snp_ld[snp][n3d]
      stats  = [snp,ind_fst,n1d,n3d,fstn1d,fstn3d,ldn1d,ldn3d]
      print '\t'.join([str(stat) for stat in stats])
      res.append(stats)
  except Exception as e:
    print "failed on snp: %s"%snp
    print "error:\n%s"%e
    print "entering debug..."
    pdb.set_trace()
    sys.exit(1)
  return res

def read_pop(pop_id_map):
  with open(pop_id_map,'rb') as fin:
    pop_id_dict = {}
    # population : [sample_id,...]
    reader = csv.reader(fin,delimiter='\t')
    for row in reader:
      pop_id_dict.setdefault(row[1],[]).append(row[0])
  populations = pop_id_dict.keys()
  print "\n# Populations: %s"%" ".join([str((i,pop)) for i,pop in enumerate(populations)])
  return pop_id_dict

def read_map(pedmap):
  snp_map = {}
  map_file = "%s-tab1234.map"%pedmap
  with open(map_file,'rb') as fin:
    reader = csv.reader(fin,delimiter='\t')
    for row in reader:
      # SNP name -> chr,bp
      snp_map[row[1]] = (row[0],row[3])
  return snp_map

def qc_freq(snp_map,pedmap):
  # If LD has been calculated before, load from log, else compute
  if not os.path.isfile("plink.frq %s.frq"%pedmap):
    try:
      cmd = "plink --file %s --freq --noweb"%pedmap
      p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      out,err = p.communicate() # ensure process has finished
      os.system("cp -f plink.frq %s.frq"%pedmap)
    except: raise
    finally:
      os.system("rm -f plink.*") # remove plink output files
  else:
    print "# Loading frequency from log..."
  with open("%s.frq"%pedmap,'rb') as fin:
    header = fin.readline()
    for line in fin:
      row = line.split()
      snp,frq = row[1],float(row[4])
      if frq < 0.0000001: # "zero"
        del snp_map[snp]
    return snp_map
  
def read_groups(snp_groups,snp_map):
  snp_loc = {}
  all_snps = set([])
  union_snps  = set([])
  if snp_groups:
    group_snp_dict = {}
    # pdbid : [snp,...]
    with open(snp_groups,'rb') as fin:
      fin.readline() # burn the header
      # gene,snp_name,pdbid,chain,x,y,z
      reader = csv.reader(fin,delimiter='\t')
      for row in reader:
        all_snps.add(row[1])
        # Only read SNPs present in the pedmap
        if row[1] in snp_map:
          # Only read unique SNPs (intentional duplication in file)
          if row[1] not in group_snp_dict.setdefault(row[2],[]):
            group_snp_dict[row[2]].append(row[1]) # store the pdbid-snp association
          if row[1] not in group_snp_dict.setdefault(row[0],[]):
            group_snp_dict[row[0]].append(row[1]) # store the gene-snp association
          # Record the SNP (qualified by pdbid and chain) 3D location
          snp_loc.setdefault(row[1],{})[(row[2],row[3])] = [float(x) for x in row[4:]]
          union_snps.add(row[1])
  print "# %d SNPs found in PDBMap."%len(all_snps)
  # Pull those SNPs 
  snp_list = list(union_snps)
  return group_snp_dict,snp_loc,snp_list

def calc_ld_dist(pedmap,snp_list):
  snp_bp_dist  = {}
  snp_ld = {}

  # If LD has been calculated before, load from log, else compute
  if not os.path.isfile("%s.ld"%pedmap):
    try:
      print "# Calculating LD with plink..."
      cmd = "plink --file %s --r2 --ld-window 250000000 --ld-window-kb 250000 --ld-window-r2 0.0 --noweb"%pedmap
      p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      out,err = p.communicate() # ensure process has finished
      os.system("cp -f plink.ld %s.ld"%pedmap)
    except: raise
    finally:
      os.system("rm -f plink.*") # remove plink output files
  else:
    print "# Loading LD from log..."
  with open("%s.ld"%pedmap,'rb') as fin:
    header = fin.readline()
    # CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
    for line in fin:
      row = line.split()
      # Only include pedmap-pdbmap union SNPs
      # and if the SNPs are on the same chromosome
      if row[2] in snp_list and row[5] in snp_list and row[0]==row[3]:
        row[1],row[4],row[6] = float(row[1]),float(row[4]),float(row[6])
        # Distances
        snp_bp_dist.setdefault(row[2],{})[row[5]] = math.fabs(row[4]-row[1])
        snp_bp_dist.setdefault(row[5],{})[row[2]] = math.fabs(row[4]-row[1])
        # LD (R^2)
        snp_ld.setdefault(row[2],{})[row[5]] = row[6]
        snp_ld.setdefault(row[5],{})[row[2]] = row[6]

  snpAs = snp_bp_dist.keys()
  for snpA in snpAs:
    snpBs = snp_bp_dist[snpA].keys()

  print "# # LD SNPs:",len(snp_ld)
  print "# # 1D SNPs:",len(snp_bp_dist)
  return snp_ld,snp_bp_dist
    

def calc_global_fst(snp_list,pop_id_map,pedmap):
  snp_fstat = {}
  tup = snp_list
  # Write N-tuples to temp file
  tempfile = "%d.temp"%multidigit_rand(5)
  try:
    with open(tempfile,'wb') as fout:
      fout.write(' '.join(tup))
    # Calculate Fst
    cmd = ["calc_fst.r",pop_id_map,pedmap,tempfile]
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    parser = process_parser(p)
    global_fst,ind_fst,pop_freqs = parse_all(parser)
    # Collapse Fst and Frequencies into one dictionary
    snp_popgen = {}
    print "# # Fst scores:",len(ind_fst)
    print "# # Frequencies:",len(pop_freqs)
    keys = pop_freqs.keys()
    for key in keys:
      snp_popgen[key] = (ind_fst[key],pop_freqs[key])

    snps = snp_popgen.keys()
    for snp in snps:
      snp_fstat[snp] = {1:snp_popgen[snp][0]}
    print "# # Fst SNPs:",len(snp_fstat)
    print "# Global Fst: ",global_fst
    return snp_fstat # not necessary, but improves clarity
  except: raise
  finally:
    os.system("rm -f %s"%tempfile) # delete temp file

def calc_kwise_fst(snp_fstat,pop_id_map,pedmap,krange,snp_groups=None,group_snp_dict=None,snp_list=None):
  # Evaluate k-tuples for all k in krange
  for k in krange:

    # If these k-tuples have been computed before, load k-wise Fst from file
    if os.path.isfile('%s_%d-tuple.fst'%(pedmap,k)):
      print "# Loading %d-wise Fst from log..."%k
      with open('%s_%d-tuple.fst'%(pedmap,k),'rb') as fin:
        reader = csv.reader(fin,delimiter='\t')
        for row in reader:
          snp_fstat[row[0]].setdefault(k,{})[tuple(row[2:])] = float(row[1])
      continue # to next k
    with open('%s_%d-tuple.fst'%(pedmap,k),'wb'): pass # create k-wise Fst file

    # Generate k-tuples
    print "# Calculating tuples for k=%d..."%k
    combns = []
    if snp_groups:
      for group in group_snp_dict.iterkeys():
        new_combns = [x for x in combinations(group_snp_dict[group],k)]
        combns.extend(new_combns)
    elif snp_list:
        combns = [x for x in combinations(snp_list,k)]
    else: raise Exception("Must provide snp_list or snp_groups")

    # Some snp-snp duplication possible between groups, reduce to unique
    unique_combns = [tup for tup in combns if tup[0] < tup[1]]
    unique_combns.extend([tup for tup in combns if tup[1] < tup[0] and tup not in unique_combns])
    del combns # immediately free memory

    # Write the k-tuples to temp file
    tempfile = "%d.temp"%multidigit_rand(5)
    try:
      with open(tempfile,'wb') as fout:
          writer = csv.writer(fout,delimiter=' ')
          for combn in unqiue_combns:
            writer.writerow(combn)

      # Calculate Fst
      cmd = ["calc_fst.r",pop_id_map,pedmap,tempfile]
      p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      parser = process_parser(p)
      tuple_fstat = parse_tuples(parser,unique_combns)
      # Store all k-wise data into snp_fstat and log
      with open('%s_%d-tuple.fst'%(pedmap,k),'a') as fout:
        tups = tuple_fstat.keys()
        for tup in tups:
          fst = tuple_fstat[tup][6]
          for snp in tup:
            partners = tuple([x for x in tup if x!=snp])
            if snp in snp_fstat:
              snp_fstat[snp].setdefault(k,{})[partners] = fst
              fout.write("%s\t%f\t%s\n"%(snp,fst,"\t".join(partners)))
    except: raise
    finally:
      os.system("rm -f %s"%tempfile)
  
  return snp_fstat

def calc_3d_dist(snp_loc):
  snp_3d_dist = {}
  snp_pairs = [x for x in combinations(snp_loc.keys(),2)]
  for snp_pair in snp_pairs:
    snpA_locs = snp_loc[snp_pair[0]]
    snpB_locs = snp_loc[snp_pair[1]]
    # only retain pairs if the pdbids match
    loc_pairs = [x for x in product(snpA_locs.iteritems(),snpB_locs.iteritems()) if x[0][0][0]==x[1][0][0]]
    for loc_pair in loc_pairs:
      # (snpA,pdb-chain,(x,y,z)) , (snpB,pdb-chain,(x,y,z))
      xd = loc_pair[0][1][0]-loc_pair[1][1][0]
      yd = loc_pair[0][1][1]-loc_pair[1][1][1]
      zd = loc_pair[0][1][2]-loc_pair[1][1][2]
      d  = math.sqrt(xd*xd + yd*yd + zd*zd)
      # snpA -> snpB -> (distance,loc_pair)
      snp_3d_dist.setdefault(snp_pair[0],{})
      snp_3d_dist[snp_pair[0]].setdefault(snp_pair[1],[])
      # Add distance, keep locations for reference
      snp_3d_dist[snp_pair[0]][snp_pair[1]].append((d,loc_pair))
  lonely_snps = [snp for snp in snp_loc if snp not in snp_3d_dist]
  # Dummy code an NA partner and NaN distance for SNPs with no 3D partners
  for snp in lonely_snps:
    snp_3d_dist[snp] = {"NA":[(float("NaN"),None)]} 
  print "# # 3D SNPs:",len(snp_3d_dist)
  return snp_3d_dist

def parse_all(parser):
  for line in parser:
    line = line.strip()
    if not line:
      continue
    if line[0] == "$":
      pop_freqs = parse_freq(parser,line[1:])
    elif line == "Fst":
      ind_fst = parse_fst(parser)
    elif line.split(' ')[0] == "Ho":
      global_fstats = parse_overall(parser)
  return global_fstats,ind_fst,pop_freqs

def parse_tuples(parser,combns):
  results = {}
  combn_count = 0
  for line in parser:
    line = line.strip()
    if not line or line == "Fst":
      continue
    elif line == "##":
      break
    elif line.split(' ')[0] == "Ho":
      continue
    if len(line.split()) == 1:
      results[combns[combn_count]] = [float("NaN")]*11
    else:
      results[combns[combn_count]] = [max(float(x),0) for x in line.split()]
    combn_count += 1
  # Ho Hs Ht Dst Htp Dstp Fst Fstp Fis Fis Dest
  return results

def parse_overall(parser):
  for line in parser:
    line = line.strip()
    if not line or line == "Fst":
      continue
    elif line == "##":
      break
    elif line.split(' ')[0] == "Ho":
      continue
    results = [float(x) for x in line.split(' ')]
  # Ho Hs Ht Dst Htp Dstp Fst Fstp Fis Fis Dest
  return results

def parse_fst(parser):
  results = {}
  for line in parser:
    line = line.strip()
    if not line:
      continue
    elif line == "##":
      break
    snp_name,fst = line.split()
    snp_name  = snp_name
    results[snp_name] = max(float(fst),0)
  return results

def parse_freq(parser,snp_name):
  results = {}
  pop_freqs = {}
  snp_name = snp_name.replace('.',':')
  for line in parser:
    line = line.strip()
    if not line:
      continue
    elif line == "##":
      if pop_freqs:
        # save current SNP
        results[snp_name] = pop_freqs
      break
    elif line[0] == '$':
      if pop_freqs:
        # save current SNP
        results[snp_name] = pop_freqs
      # begin parsing next SNP
      snp_name  = line[1:].replace('.',':')
      pop_freqs = {}
    elif line[0] == 'x':
      continue
    else:
      line      = line.strip()
      pop_freqs[line[0]] = line[1:].strip().split(' ')
  if pop_freqs:
    # save current SNP
    results[snp_name] = pop_freqs
  return results

# def string_parser(content):
def process_parser(p):
  while True:
    line = p.stdout.readline()
    if not line: break
    yield line

def tup_eq(a,b):
  return all(x in b for x in a) and all(x in a for x in b)

def preprocess_pedmap(pedmap):
  os.system("""plink --noweb --file %s --allele1234 --tab --recode --out %s-tab1234"""%(pedmap,pedmap))
  os.system("""sed -e 's/[ ]//g' %s-tab1234.ped > %s-tab1234TEMP.ped"""%(pedmap,pedmap))
  os.system("""mv -f %s-tab1234TEMP.ped %s-tab1234.ped"""%(pedmap,pedmap))

# Support function for temp files
def multidigit_rand(digits):
  import random
  randlist = [random.randint(1,9) for i in xrange(digits)]
  multidigit_rand = int(''.join([str(x) for x in randlist]))
  return multidigit_rand

if __name__ == "__main__":
  main()
