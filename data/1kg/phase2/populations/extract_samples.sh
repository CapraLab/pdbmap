#!/bin/sh
#
# Project        : PDBMap
# Filename       : extract_samples.sh
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-01-30
# Description    : Extracts the sample IDs for each contiental ancestry
#                : and HapMap population from the 1KG panel file.
#=============================================================================#
# PBS Parameters
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l nodes=vision.mc.vanderbilt.edu
#PBS -l mem=5000mb
#PBS -l walltime=1:00:00:00
#=============================================================================#

# Continental Ancestries
grep AFR phase1_integrated_calls.20101123.ALL.panel | cut -f1 > AFR.samples.list
grep AMR phase1_integrated_calls.20101123.ALL.panel | cut -f1 > AMR.samples.list
grep ASN phase1_integrated_calls.20101123.ALL.panel | cut -f1 > ASN.samples.list
grep EUR phase1_integrated_calls.20101123.ALL.panel | cut -f1 > EUR.samples.list

# HapMap Populations
grep ASW phase1_integrated_calls.20101123.ALL.panel | cut -f1 > ASW.samples.list
grep CEU phase1_integrated_calls.20101123.ALL.panel | cut -f1 > CEU.samples.list
grep CHB phase1_integrated_calls.20101123.ALL.panel | cut -f1 > CHB.samples.list
grep CHS phase1_integrated_calls.20101123.ALL.panel | cut -f1 > CHS.samples.list
grep CLM phase1_integrated_calls.20101123.ALL.panel | cut -f1 > CLM.samples.list
grep FIN phase1_integrated_calls.20101123.ALL.panel | cut -f1 > FIN.samples.list
grep GBR phase1_integrated_calls.20101123.ALL.panel | cut -f1 > GBR.samples.list
grep IBS phase1_integrated_calls.20101123.ALL.panel | cut -f1 > IBS.samples.list
grep JPT phase1_integrated_calls.20101123.ALL.panel | cut -f1 > JPT.samples.list
grep LWK phase1_integrated_calls.20101123.ALL.panel | cut -f1 > LWK.samples.list
grep MXL phase1_integrated_calls.20101123.ALL.panel | cut -f1 > MXL.samples.list
grep PUR phase1_integrated_calls.20101123.ALL.panel | cut -f1 > PUR.samples.list
grep TSI phase1_integrated_calls.20101123.ALL.panel | cut -f1 > TSI.samples.list
grep YRI phase1_integrated_calls.20101123.ALL.panel | cut -f1 > YRI.samples.list
