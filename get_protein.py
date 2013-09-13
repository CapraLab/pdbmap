#!/usr/bin/python27
def load_transmap(transmap_file):
	with open(transmap_file,'rU') as fin:
		reader = csv.reader(fin,delimiter='\t')
		transmap = {}
		for (unp,id_type,trans) in reader:
			unp = unp.strip()
			if trans in transmap:
				transmap[trans].append((unp,id_type))
			else:
				transmap[trans] = [(unp,id_type)]
	return transmap
## -------------------------------------------------- ##
import sys,csv,subprocess

var_file = sys.argv[1]

transmap_file = '/labs/twells/sivleyrm/pdbmap/transmap.dat'
transmap = load_transmap(transmap_file)

perlout = subprocess.check_output(['./get_protein.pl',var_file])
perl_split = perlout.split('\n')
vars_with_trans = [line.split('\t') for line in perl_split if line]

writer = csv.writer(sys.stdout,delimiter='\t')
for row in vars_with_trans:
	chrom = row[0]
	start = row[1]
	end   = row[2]
	name  = row[3]
	trans = row[4]
	if trans in transmap:
		for (unp,id_type) in transmap[trans]:
			writer.writerow([chrom,start,end,name,trans,unp])