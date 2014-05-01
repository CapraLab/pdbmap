#!/usr/bin/python27

# Improvements:
#
# 1. Update so that it loops through all PDBIDs in the scores file
# 2. Dynamically fetch PDB structures in the loop
# 3. Document the required format for the scores file
# 4. Automatically center on high/low scored regions and generate png

import sys,csv,math
if __name__ != '__main__':
	from pymol import cmd
max_dist = -999	# Used by show_density and dist

def overwrite_bfactor(pdbid,chain,resi,value):
	pdbid = pdbid.lower()
	if chain == 'A':
		# Some modbase models do not have chains, and have been dummy coded to A
		selection = "(%s and (chain %s or chain '') and resi %d)"%(pdbid,chain,resi)
	else:
		selection = "(%s and chain %s and resi %d)"%(pdbid,chain,resi)
	exp = "b=%f"%value
	command = "alter %s, %s"%(selection,exp)
	if __name__ != "__main__":
		cmd.alter(selection,exp)
	print command
	return command

def reset_bfactor(pdbid,val=0.0):
	exp = "b=%d"%val
	cmd.alter(pdbid,exp)

def overwrite_bfactors(pdbid,scores,col=4,resis=None,discrete=False,binary=False,displaytype='cartoon',surface=False,var_spheres=False,spec_range=None):
	col -= 1 # correct for 0-indexing
	pdbid = pdbid.lower()
	if __name__ == '__main__':
		cmd.fetch(pdbid)
	if isinstance(scores,str): # if filename, read file
		fin = open(scores,'rU')
		reader = csv.reader(fin,delimiter='\t')
		if resis:
			scores = [row for row in reader if row[0].lower()==pdbid and int(row[2]) in resis]
		else:
			scores = [row for row in reader if row[0].lower()==pdbid]
		fin.close()
	for score in scores:
		if not score[col] or score[col] == 'NULL':
			score[col] = -1
	uniq_scores = [float(score[col]) for score in scores if float(score[col]) > 0]
	if spec_range and uniq_scores:
		min_score = spec_range[0]
		max_score = spec_range[1]
	elif uniq_scores:
		min_score = min(uniq_scores) - 0.0001
		max_score = max(uniq_scores)
	else:
		msg = "ERROR: %s contains no 1000 Genomes variants."%pdbid
		raise Exception(msg)

	if __name__ != "__main__":
		reset_bfactor(pdbid,min_score-10)
	commands = [overwrite_bfactor(row[0],row[1],int(row[2]),float(row[col])) for row in scores]
	include = "br. b > %f"%min_score
	exclude = "br. b < %f"%min_score
	print("Minimum score: %s"%min_score)
	print("Maximum score: %s"%max_score)

	# Generate models
	if surface:
		cmd.create("surface_obj",pdbid)
		cmd.color("grey","surface_obj")
	cmd.hide("everything","all")
	cmd.show(representation=displaytype)
	cmd.set(name="cartoon_discrete_colors",value="on")
	if var_spheres:
			cmd.show(selection="name CA and br. b > -1.1",representation="spheres")
			cmd.set(name="sphere_scale",value=1)
	if surface:
		cmd.show(selection="surface_obj",representation="surface")
		cmd.set(selection="surface_obj",name="transparency",value=0.5)

	# Assign colors
	cmd.bg_color("white")
	cmd.color("grey","all")
	if binary:
		cmd.color("red", selection="br. b=1")
	elif discrete:
		discrete_colors = ["","firebrick","density","forest","tv_yellow","deeppurple","orange","orange","deepteal","white","black","grey"]
		for i in range(1,int(max_score)+1):
			cmd.color(discrete_colors[i],selection="br. b=%d"%i)
	else:
		cmd.spectrum("b","blue_red",selection=include,minimum=min_score,maximum=max_score)
	cmd.color("grey",exclude)
	cmd.orient("all")
	if __name__ == '__main__':
		cmd.save("%s.annotated.pdb"%pdbid)
	return commands

def show_density(pdbid,pdbmap_file,variants_file,radius):
	radius = float(radius)
	if not show_density.dens_mat:
		# columns: chain,seqres,x,y,z
		fin = open(pdbmap_file,'rU')
		fin.readline() # burn the header
		reader = csv.reader(fin,delimiter='\t')
		residues = [row[1:] for row in reader if row[0].lower()==pdbid]
		fin = open(variants_file,'rU')
		fin.readline() # burn the header
		reader = csv.reader(fin,delimiter='\t')
		variants = [row[1:] for row in reader if row[0].lower()==pdbid]

		# Compute the distance "matrix"
		dist_mat = {}
		for res_chain,res_seqres,x,y,z in residues:
			res_pos = (x,y,z)
			for var_chain,var_seqres,i,j,k,chrom,start,end,var_name in variants:
				var_pos = (i,j,k)
				if (res_chain,res_seqres) not in dist_mat:
					dist_mat[(res_chain,res_seqres)] = []
				dist_mat[(res_chain,res_seqres)].append(dist((x,y,z),(i,j,k)))

		# Convert the distance matrix into a density "matrix"
		show_density.dens_mat = []
		for chain,seqres in dist_mat:
			distances = dist_mat[(chain,seqres)]
			for k in range(int(max_dist)):
				num_neighbors = len([x for x in distances if x <= k])
				show_density.dens_mat.append((chain,seqres,k,num_neighbors)) 
					
	# Reduce to requested k
	k_dense_map = [[pdbid,row[0],row[1],row[3]] for row in show_density.dens_mat if row[2]<=radius]
	overwrite_bfactors(pdbid,k_dense_map,displaytype='mesh')

# Once computed, store permanently, extracting different k values
show_density.dens_mat = None

def dist(ipos,jpos):
	global max_dist
	dist_x = abs(float(ipos[0]) - float(jpos[0])) ** 2
	dist_y = abs(float(ipos[1]) - float(jpos[1])) ** 2
	dist_z = abs(float(ipos[2]) - float(jpos[2])) ** 2
	xyz_dist = math.sqrt(dist_x+dist_y+dist_z)
	max_dist = max(max_dist,xyz_dist)
	return xyz_dist

if __name__=="__main__":
	sys.argv = [x.lower() for x in sys.argv]
	commands = overwrite_bfactors(sys.argv[1],sys.argv[2])
	print "Commands:"
	for command in commands:
		print command
else:
	cmd.extend("overwrite_bfactors",overwrite_bfactors)
