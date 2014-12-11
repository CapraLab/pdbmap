CREATE TABLE IF NOT EXISTS GenomicData (
# Standard columns
label VARCHAR(100), # Dataset label
chr VARCHAR(50),
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
da VARCHAR(50), # Derived Allele
maf DOUBLE,     # Allele Frequency: Global (AC/AN)
amr_af DOUBLE,  # Allele Frequency: American
asn_af DOUBLE, # Allele Frequency: Asian
afr_af DOUBLE, # Allele Frequency: African
eur_af DOUBLE, # Allele Frequency: European
eas_af DOUBLE, # Allele Frequency: East Asian
sas_af DOUBLE, # Allele Frequency: South Asian
ens_gene VARCHAR(50), # Ensembl Gene identifier
hgnc_gene VARCHAR(50), # HGNC Gene identifier
snpsource VARCHAR(50), # Low coverage or Exome?
PRIMARY KEY(label,name,chr,start,end),
KEY(label,chr,start,end),
KEY(label,maf),
KEY(label,amr_af),
KEY(label,asn_af),
KEY(label,afr_af),
KEY(label,eur_af),
KEY(label,eas_af),
KEY(label,sas_af)
)