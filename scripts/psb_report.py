#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : psb_report.py
# Authors        : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2017-01-22
# Description    : Generates structure and spatial reports 
#                : for PSB collaborations.
#=============================================================================#
# from __future__ import print_function
import pandas as pd
pd.set_option("display.max_columns",100)
pd.set_option("display.width",1000)
import numpy as np
from time import sleep
import os,sys,argparse,ConfigParser
from glob import glob
sys.path.append("..")
from lib.amino_acids import longer_names
from lib.PDBMapIO import PDBMapIO
from lib.PDBMapProtein import PDBMapProtein
from slurm import SlurmJob
from weasyprint import HTML,CSS
from jinja2 import Environment, FileSystemLoader


conf_parser = argparse.ArgumentParser(add_help=False)
conf_parser.add_argument("-c","--config",
help="PDBMap configuration profile for database access", metavar="FILE")
args, remaining_argv = conf_parser.parse_known_args()
defaults = {
	"dbhost"       : "chgr2.accre.vanderbilt.edu",
	"dbname"       : "pdbmap_v12",
	"dbuser"       : "script_access",
	"dbpass"       : "capralab",
	"slabel"       : "uniprot-pdb",
	"idmapping"    : "/dors/capra_lab/data/uniprot/idmapping/HUMAN_9606_idmapping_UNP-RefSeq-PDB-Ensembl.tab",}
if args.config:
	config = ConfigParser.SafeConfigParser()
	config.read([args.config])
	defaults.update(dict(config.items("Genome_PDB_Mapper")))

desc  = "Given a gene name or UniProt identifier and mutation, "
desc += "this script will generate a report detailing the structural "
desc += "coverage and various structural and spatial properties of "
desc += "available structural and genetics datasets."
parser = argparse.ArgumentParser(description=desc)
parser.set_defaults(**defaults)
parser.add_argument("collaboration",type=str,help="Collaboration ID")
parser.add_argument("project",type=str,help="Project ID")
parser.add_argument("entity",type=str,help="Gene ID or SwissProt AC")
parser.add_argument("mutation",type=str,help="HGVS mutation string (e.g. S540A)")
parser.add_argument("--overwrite",action="store_true",default=False,help="If true, overwrite previous results.")
args = parser.parse_args()

args.outdir = "/dors/capra_lab/projects/psb_collab/%s/%s"%(args.collaboration,args.project)
if not os.path.exists(args.outdir):
	os.makedirs(args.outdir)
if not os.path.exists("%s/PathProx"%args.outdir):
	os.makedirs("%s/PathProx"%args.outdir)

def load_io(args):
	PDBMapProtein.load_idmapping(args.idmapping)
	io = PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,slabel=args.slabel)
	return io

def detect_entity(io,entity):
	print "\nAttempting to identify entity: %s..."%entity
	etype = io.detect_entity_type(entity)
	if etype == "unp":
		unp  = entity
		hgnc = PDBMapProtein.unp2hgnc(unp)
		print "UniProt: %s => HGNC: %s..."%(unp,hgnc)
	elif etype: # HGNC gene ID
		hgnc = entity
		unp  = etype
		print "HGNC: %s => UniProt: %s..."%(hgnc,unp)
	else:
		msg = "Unrecognized ID. Gene name or UniProt/Swiss-Prot AC required.\n"
		sys.stderr.write(msg); sys.exit(1)
	return unp,hgnc

def get_structures(io,unp,mutation):
	sql = """
	SELECT IF(c.method IS NOT NULL,c.method,d.method) as method,a.structid,a.chain,
					COUNT(distinct b.seqid) as nresidues,
					CONCAT(CAST(min(b.seqid) as CHAR),"-",CAST(max(b.seqid) AS CHAR)) AS coverage,
					c.resolution, d.identity,
					UPPER(d.pdbid) as template
	FROM Chain a
	INNER JOIN Residue b
	ON a.label=b.label AND a.structid=b.structid AND a.chain=b.chain AND a.biounit=b.biounit AND a.model=b.model
	LEFT JOIN Structure c
	ON a.label=c.label AND a.structid=c.pdbid
	LEFT JOIN Model d
	ON a.label=d.label AND a.structid=d.modelid
	WHERE a.biounit=0 AND a.model=0 AND a.label=%s AND a.unp=%s
	GROUP BY a.structid,a.chain
	HAVING %s BETWEEN MIN(b.seqid) AND MAX(b.seqid)
	ORDER BY c.pdbid IS NOT NULL,nresidues DESC,c.resolution ASC,d.identity DESC
	"""
	# res = list(io.secure_query(sql,(io.slabel,unp.upper(),mut_pos),cursorclass='Cursor'))
	try:
		mut_pos = int(mutation[1:-1])
	except:
		# Convert mutation argument to 1-letter codes
		args.mutation = ''.join([longer_names[mutation[:3].upper()],mutation[3:-3],longer_names[mutation[-3:].upper()]])
		mut_pos = int(args.mutation[1:-1])
		
		# FIXME: Incorrectly
	df  = pd.read_sql(sql,io._con,params=(io.slabel,unp.upper(),mut_pos))
	# df  = pd.DataFrame(res,columns=["structid","chain"])
	return df

# Identify the structures available for this gene
print "Identifying structures of %s with coverage of %s..."%(args.entity,args.mutation)
io = load_io(args)
unp,gene = detect_entity(io,args.entity)
df = get_structures(io,unp,args.mutation)
df.columns = ["Method","Structure ID","Chain",
							"Residues","Coverage","Resolution (PDB)",
							"Sequence Identity (ModBase)","PDB Template (ModBase)"]
df = df.sort_values(by=["Method","Residues","Resolution (PDB)","Sequence Identity (ModBase)"],
						ascending=[False,False,True,False])

def generate_initial_report(df):
	""" Generate a report detailing available structural coverage """
	# Create a multiindex for viewing
	tdf = df.set_index(["Method","Structure ID","Chain"])
	env = Environment(loader=FileSystemLoader('.'))
	template = env.get_template("psb_report.html")
	template_vars = {"title" : "PSB Structure Report - %s"%gene,
									 "gene"  : "%s (%s)"%(gene,unp),
									 "psb_structure_results" : tdf.to_html()}
	html_out = template.render(template_vars)
	fname = "%s/%s_structure_initial_report.pdf"%(args.outdir,gene)
	print "\nWriting initial report to %s..."%fname
	HTML(string=html_out).write_pdf(fname,stylesheets=["typography.css"])#,CSS(string="@page { size: A3 landscape; }")])
# Generate an initial report with structure information
generate_initial_report(df)

# Submit and track PathProx analyses
def launch_job(structid,chain,mutation):
	if not os.path.exists("../qsub/psb_pathprox"):
		os.makedirs("../qsub/psb_pathprox")
	params = {"collab":args.collaboration,"project":args.project,"pdbid":structid,"chain":chain,"mutation":mutation}
	return SlurmJob(open("psb_template.slurm",'rb').read()%params)
def processed(sid,chain):
	# Return true if both ClinVar and COSMIC complete flags are present
	return list(glob("%s/PathProx/%s_%s/*.complete"%(args.outdir,sid,chain)))>1	
jobs = [launch_job(sid,ch,args.mutation) for sid,ch in df[["Structure ID","Chain"]].values if args.overwrite or not processed(sid,ch)]
if jobs:
	print "\nLaunched %d PathProx jobs to SLURM queue...\n"%len(jobs)
	ds = [j.get_info() for j in jobs]
	dm = {}
	for k in ["JobId","JobName","JobState","RunTime","NodeList"]:
		dm[k] = [d[k] for d in ds]
	dq = pd.DataFrame.from_dict(dm)
	sys.stdout.write('%s\n\n'%dq.to_string())
	while not all([j.is_finished() for j in jobs]):
		sys.stdout.write("\rWaiting on jobs to finish... (%d/%d complete)"%(sum([j.is_finished() for j in jobs]),len(jobs)))
		sys.stdout.flush()
		sleep(15) # Sleep for 15 seconds
	print ""
	if not all([j.is_complete() for j in jobs]):
		msg = "\nWarning: Not all jobs completed without error.\n"
		sys.stderr.write(msg)

# Compile PathProx results into report (one page per structure+chain+mutation)
def structure_detail(df,structid,chain,mutation):
	""" Generate the HTML for a single page structure report """
	mut_pos = int(mutation[1:-1])
	# Subset the dataframe to the results of this chain
	structure  = {"structid":structid,"chain":chain,"mutation":mutation}
	structure["results"] = df[(df["Structure ID"]==structid) & (df["Chain"]==chain)\
													].set_index(["Method","Structure ID","Chain"]).to_html()
	pp_dir     = "%s/PathProx/%s_%s"%(args.outdir,structid,chain)
	cln_prefix = "%s_%s_clinvar_exac_D"%(structid,chain)
	csm_prefix = "%s_%s_cosmic_exac_D"%(structid,chain)
	# Load the ClinVar-ExAC PathProx results for this mutation
	if os.path.exists("%s/%s_pathprox.txt"%(pp_dir,cln_prefix)):
		pp = pd.read_csv("%s/%s_pathprox.txt"%(pp_dir,cln_prefix),sep='\t',header=0)
		structure["pp_cln_mut"] = pp.ix[pp["unp_pos"]==mut_pos,"pathprox"].values[0]
	else:
		structure["pp_cln_mut"] = np.nan # PathProx analysis did not run for ClinVar
	# Load the COSMIC-ExAC PathProx results for this mutation
	if os.path.exists("%s/%s_pathprox.txt"%(pp_dir,csm_prefix)):
		pp = pd.read_csv("%s/%s_pathprox.txt"%(pp_dir,csm_prefix),sep='\t',header=0)
		structure["pp_csm_mut"] = pp.ix[pp["unp_pos"]==mut_pos,"pathprox"].values[0]
	else:
		structure["pp_csm_mut"] = np.nan # PathProx analysis did not run for COSMIC
	# Load the variant plots
	if os.path.exists("%s/%s_neutral.png"%(pp_dir,cln_prefix)):
		structure["exac_var"] = "%s/%s_neutral.png"%(pp_dir,cln_prefix)
	else:
		structure["exac_var"] = "%s/%s_neutral.png"%(pp_dir,csm_prefix)
	structure["cln_var"]    = "%s/%s_pathogenic.png"%(pp_dir,cln_prefix)
	structure["csm_var"]    = "%s/%s_pathogenic.png"%(pp_dir,csm_prefix)
	# Load the Ripley's K results for ClinVar and COSMIC
	if os.path.exists("%s/%s_neutral_K_plot.png"%(pp_dir,cln_prefix)):
		structure["exac_K"]   = "%s/%s_neutral_K_plot.png"%(pp_dir,cln_prefix)
	else:
		structure["exac_K"]   = "%s/%s_neutral_K_plot.png"%(pp_dir,csm_prefix)
	structure["cln_K"]      = "%s/%s_pathogenic_K_plot.png"%(pp_dir,cln_prefix)
	structure["csm_K"]      = "%s/%s_pathogenic_K_plot.png"%(pp_dir,csm_prefix)
	# Load the Ripley's D results for ClinVar and COSMIC
	structure["cln_D"]      = "%s/%s_D_plot.png"%(pp_dir,cln_prefix)
	structure["csm_D"]      = "%s/%s_D_plot.png"%(pp_dir,csm_prefix)
	# Load the PathProx Mapping and Performance for ClinVar
	structure["cln_pp"]     = "%s/%s_pathprox.png"%(pp_dir,cln_prefix)
	structure["cln_roc"]    = "%s/%s_pathprox_roc.png"%(pp_dir,cln_prefix)
	structure["cln_pr"]     = "%s/%s_pathprox_pr.png"%(pp_dir,cln_prefix)
	# Load the PathProx Mapping and Performance for COSMIC
	structure["csm_pp"]     = "%s/%s_pathprox.png"%(pp_dir,csm_prefix)
	structure["csm_roc"]    = "%s/%s_pathprox_roc.png"%(pp_dir,csm_prefix)
	structure["csm_pr"]     = "%s/%s_pathprox_pr.png"%(pp_dir,csm_prefix)
	return structure

def generate_final_report(df):
	""" Generate the HTML for the final report """
	tdf = df.set_index(["Method","Structure ID","Chain"])
	env = Environment(loader=FileSystemLoader('.'))
	template = env.get_template("psb_report.html")
	# Generate vars for individual structure pages
	structures = [structure_detail(df,sid,ch,args.mutation) for sid,ch in df[["Structure ID","Chain"]].values]
	# Generate vars for overall html template
	template_vars = {"title" : "PSB Structure Report - %s"%gene,
									 "gene"  : "%s (%s)"%(gene,unp),
									 "psb_structure_results" : tdf.to_html(),
									 "Structure_Detail" : structures}
	html_out = template.render(template_vars)
	fname = "%s/%s_structure_final_report.pdf"%(args.outdir,gene)
	print "\nWriting final report to %s..."%fname
	HTML(string=html_out).write_pdf(fname,stylesheets=["typography.css"])#,CSS(string="@page { size: A3 landscape; }")])
generate_final_report(df)
