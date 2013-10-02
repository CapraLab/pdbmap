#!/bin/sh

# $1 - PDB ID
# $2 - Variant Set
# $3 - Density radius

pdbid=`echo $1 | tr '[:upper:]' '[:lower:]'`

var_res_file="/scratch/sivleyrm/pdbmap/scratch/${pdbid}_variant_residues.txt";
var_pos_file="/scratch/sivleyrm/pdbmap/scratch/${pdbid}_variant_positions.txt";
pdbmap_file="/scratch/sivleyrm/pdbmap/scratch/${pdbid}_pdbmap.txt";

echo "Running MySQL queries...";
mysql -hgwar-dev -umike -pcheezburger -e "USE pdbmap_v7; SELECT DISTINCT pdbid,chain,chain_seq,1 FROM Intersect_Variants_$2 WHERE pdbid='$pdbid' ORDER BY chain,chain_seq;" > $var_res_file
mysql -hgwar-dev -umike -pcheezburger -e "USE pdbmap_v7; SELECT DISTINCT pdbid,chain,chain_seq,x,y,z FROM GenomePDB WHERE pdbid='$pdbid' ORDER BY chain,chain_seq;" > $pdbmap_file
mysql -hgwar-dev -umike -pcheezburger -e "USE pdbmap_v7; SELECT DISTINCT pdbid,chain,chain_seq,x,y,z FROM Intersect_Variants_$2 WHERE pdbid='$pdbid' ORDER BY chain,chain_seq;" > $var_pos_file
echo "MySQL queries complete";

cmd1="run overwrite_pymol.py; "
cmd1+="load /scratch/sivleyrm/pdb/structures/all/pdb/pdb$pdbid.ent,$pdbid; "
cmd1+="overwrite_bfactors('$pdbid','$var_res_file',binary=True); "
cmd1+="orient;mset 1 x360; movie.roll 1,180,1,axis=y; movie.roll 181,360,1,axis=x; "
cmd1+="save ${pdbid}_vars.pdb; save ${pdbid}_vars.pse;"
/usr/local/bin/pymol -c -d "$cmd1"

cmd2="run overwrite_pymol.py; "
cmd2+="load /scratch/sivleyrm/pdb/structures/all/pdb/pdb$pdbid.ent,$pdbid; "
cmd2+="show_density('$pdbid','$pdbmap_file','$var_pos_file',radius=$3); "
cmd2+="orient;mset 1 x360; movie.roll 1,180,1,axis=y; movie.roll 181,360,1,axis=x; "
cmd2+="save $1_$3_density.pdb; save $1_$3_density.pse;"
/usr/local/bin/pymol -c -d "$cmd2"
