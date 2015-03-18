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
eas_af DOUBLE, # Allele Frequency: East Asian
sas_af DOUBLE, # Allele Frequency: South Asian
afr_af DOUBLE, # Allele Frequency: African
eur_af DOUBLE, # Allele Frequency: European
ens_gene VARCHAR(50), # Ensembl Gene identifier
hgnc_gene VARCHAR(50), # HGNC Gene identifier
amreas_Nhat DOUBLE, # Fst Numerator: American - East Asian
amrsas_Nhat DOUBLE, # Fst Numerator: American - South Asian
amreur_Nhat DOUBLE, # Fst Numerator: American - European
amrafr_Nhat DOUBLE, # Fst Numerator: American - African
eassas_Nhat DOUBLE, # Fst Numerator: East Asian - South Asian
easeur_Nhat DOUBLE, # Fst Numerator: East Asian - European
easafr_Nhat DOUBLE, # Fst Numerator: East Asian - African
saseur_Nhat DOUBLE, # Fst Numerator: South Asian - European
sasafr_Nhat DOUBLE, # Fst Numerator: South Asian - African
eurafr_Nhat DOUBLE, # Fst Numerator: European - African
allpop_Nhat DOUBLE, # Fst Numerator: All Continental Populations
amreas_Dhat DOUBLE, # Fst Denominator: American - East Asian
amrsas_Dhat DOUBLE, # Fst Denominator: American - South Asian
amreur_Dhat DOUBLE, # Fst Denominator: American - European
amrafr_Dhat DOUBLE, # Fst Denominator: American - African
eassas_Dhat DOUBLE, # Fst Denominator: East Asian - South Asian
easeur_Dhat DOUBLE, # Fst Denominator: East Asian - European
easafr_Dhat DOUBLE, # Fst Denominator: East Asian - African
saseur_Dhat DOUBLE, # Fst Denominator: South Asian - European
sasafr_Dhat DOUBLE, # Fst Denominator: South Asian - African
eurafr_Dhat DOUBLE, # Fst Denominator: European - African
allpop_Dhat DOUBLE, # Fst Denominator: All Continental Populations
amreas_Fst DOUBLE, # Fst: American - East Asian
amrsas_Fst DOUBLE, # Fst: American - South Asian
amreur_Fst DOUBLE, # Fst: American - European
amrafr_Fst DOUBLE, # Fst: American - African
eassas_Fst DOUBLE, # Fst: East Asian - South Asian
easeur_Fst DOUBLE, # Fst: East Asian - European
easafr_Fst DOUBLE, # Fst: East Asian - African
saseur_Fst DOUBLE, # Fst: South Asian - European
sasafr_Fst DOUBLE, # Fst: South Asian - African
eurafr_Fst DOUBLE, # Fst: European - African
allpop_Fst DOUBLE, # Fst: All Continental Populations
snpsource VARCHAR(50), # Low coverage or Exome?
format VARCHAR(50), # Genotype Format
gt TEXT, # Genotypes
gd_id MEDIUMINT NOT NULL AUTO_INCREMENT, # Unique, direct-reference key
PRIMARY KEY(label,name,chr,start,end),
KEY(gd_id),
KEY(label,chr,start,end),
KEY(label,maf),
KEY(label,amr_af),
KEY(label,asn_af),
KEY(label,afr_af),
KEY(label,eur_af),
KEY(label,eas_af),
KEY(label,sas_af),
KEY(label,allpop_Fst)
)