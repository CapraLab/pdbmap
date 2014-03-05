CREATE TABLE IF NOT EXISTS GenomicData (
# Standard columns
label VARCHAR(100), # Dataset label
chr VARCHAR(10),
start INT, # Start site, always specified
end INT,   # End site, specified by END or assumed start + 1
name VARCHAR(100),  # Provided name
variation VARCHAR(100), # Known variation names
vtype VARCHAR(50),  # Variant Type
svtype VARCHAR(50), # Structural variant type
ref_allele VARCHAR(50),
alt_allele VARCHAR(50),
svlen INT, # Difference in length between ref and alt alleles
quality DOUBLE, # Not sure how this is measured
avgpost DOUBLE, # MaCH/Thunder: Average posterior probability
rsq DOUBLE,     # MaCH/Thunder: Genotype imputation quality
erate DOUBLE,   # MaCH/Thunder: Per-marker mutation rate
theta DOUBLE,   # MaCH/Thunder: Per-marker transition rate
ldaf DOUBLE,    # MLE allele frequency accounting for LD
ac INT, # Alternate Allele Count
an INT, # Total Allele Count
aa VARCHAR(50), # Ancestral Allele
maf DOUBLE,     # Allele Frequency: Global (AC/AN)
amr_af DOUBLE,  # Allele Frequency: American
asn_af DOUBLE, # Allele Frequency: Asian
afr_af DOUBLE, # Allele Frequency: African
eur_af DOUBLE, # Allele Frequency: European
ens_gene VARCHAR(50), # Ensembl Gene identifier
hgnc_gene VARCHAR(50), # HGNC Gene identifier
snpsource VARCHAR(50), # Low coverage or Exome?
PRIMARY KEY(label,name,chr,start,end),
KEY(name,chr,start,end),
KEY(chr,start,end),
KEY(maf),
KEY(amr_af),
KEY(asn_af),
KEY(afr_af),
KEY(eur_af)
)