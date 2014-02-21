CREATE TABLE IF NOT EXISTS GenomicData (
# Standard columns
name VARCHAR(100),
label VARCHAR(100),
chr VARCHAR(10),
start INT,
end INT,
ref_allele VARCHAR(10),
gene VARCHAR(20),
gene_alias TEXT, # ("Extra column") list of non-Ensembl gene IDs
feature VARCHAR(100),
feature_type VARCHAR(100),
consequence VARCHAR(100),
cdna_pos INT,
cds_pos INT,
protein_pos INT,
ref_amino_acid VARCHAR(50),
var_amino_acid VARCHAR(50),
ref_codon VARCHAR(50),
var_codon VARCHAR(50),
variation VARCHAR(100),
# "Extra" columns
aa_maf DOUBLE,
afr_maf DOUBLE,
amr_maf DOUBLE,
asn_maf DOUBLE,
ea_maf DOUBLE,
eur_maf DOUBLE,
gen_maf DOUBLE, # global (general) allele frequency
biotype VARCHAR(100), # of transcript
canonical boolean, # is canonical transcript?
ccds boolean, # is CCDS transcript?
clinical_sig VARCHAR(100), # dbSNP clinical significance
dist2trans INT, # distance to transcript
domains TEXT, # formatted list as string
ensp VARCHAR(100), # Ensembl protein ID
exon VARCHAR(10), # what is this ratio?
intron VARCHAR(10), # what is this ratio?
hgvsc VARCHAR(100), # HGVS coding sequence name
hgvsp VARCHAR(100), # HGVS protein sequence name
pubmed TEXT, # list of PMIDs citing variation
polyphen DOUBLE,
sift DOUBLE,
PRIMARY KEY(label,name,chr,start,end),
KEY(name,chr,start,end),
KEY(chr,start,end),
KEY(consequence),
KEY(gene),
KEY(ensp),
KEY(feature),
KEY(gen_maf),
KEY(aa_maf),
KEY(afr_maf),
KEY(amr_maf),
KEY(asn_maf),
KEY(ea_maf),
KEY(eur_maf)
)