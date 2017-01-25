#!/bin/sh
#
# Project        : PDBMap / PathProx
# Filename       : psb_template.sh
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Program in Personalized Structural Biology,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2016-12-15
# Description    : Runs a PathProx analysis for a PSB structure of interest
#=============================================================================#
# Slurm Parameters
#SBATCH -J PSB_PathProx_%(pdbid)s_%(chain)s
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem=5GB
#SBATCH -o ../qsub/psb_pathprox/psb_pathprox.%%A.out
#=============================================================================#

./pathprox.py -c ../config/v12.config %(pdbid)s %(mutation)s --chain=%(chain)s --add_pathogenic --add_exac --radius=D --overwrite --outdir=/dors/capra_lab/projects/psb_collab/%(collab)s/%(project)s/PathProx/%(pdbid)s_%(chain)s/;
./pathprox.py -c ../config/v12.config %(pdbid)s %(mutation)s --chain=%(chain)s --add_somatic    --add_exac --radius=D --overwrite --outdir=/dors/capra_lab/projects/psb_collab/%(collab)s/%(project)s/PathProx/%(pdbid)s_%(chain)s/;
