CREATE TABLE IF NOT EXISTS GenomicConsequence (
# Foreign Key to GenomicData
chr VARCHAR(10),
start INT, # Start site, always specified
end INT,   # End site, specified by END or assumed start + 1
name VARCHAR(100), # Provided name
# Consequence data
transcript VARCHAR(100), # Ensembl transcript ID
protein VARCHAR(100),    # Ensembl protein ID
canonical TINYINT, # Canonical transcript?
allele VARCHAR(50),
consequence VARCHAR(100),
cdna_pos INT,
cds_pos INT,
protein_pos INT,
ref_amino_acid VARCHAR(50),
alt_amino_acid VARCHAR(50),
ref_codon VARCHAR(50),
alt_codon VARCHAR(50),
polyphen DOUBLE, # Score only
sift DOUBLE,     # Score only
biotype VARCHAR(100), # of transcript
domains TEXT, # formatted list as string
gc_id MEDIUMINT NOT NULL AUTO_INCREMENT, # Unique, direct-reference key
PRIMARY KEY(transcript,chr,start,end),
KEY(gc_id),
KEY(chr,start,end),
KEY(consequence),
KEY(protein),
KEY(polyphen),
KEY(sift)
)