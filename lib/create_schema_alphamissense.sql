-- AlphaMissense predicts the pathogenicty of missense variants for canonical, and non-canonical 
-- transcript isoforms.  The data arrive in three .tsv.gz files which we directly mysqlimport into
-- corresponding sql tables
-- Source data: https://console.cloud.google.com/storage/browser/dm_alphamissense
--     AlphaMissense_aa_substitutions.tsv
--     AlphaMissense_isoforms_hg38.tsv.gz
--     AlphaMissense_isoforms_aa_substitutions.tsv.gz
--
-- =============================================================
-- After you run the "CREATE TABLE..." commands below, load the sql databases from the raw downloaded alphamissense
-- text files with these respective commands:
-- $ mysqlimport --ignore --ignore-lines=4 --verbose --fields-terminated-by='\t' --local -p pdbmap_v14 AlphaMissense_aa_substitutions.tsv
-- $ mysqlimport --ignore --ignore-lines=4 --verbose --fields-terminated-by='\t' --local -p pdbmap_v14 AlphaMissense_isoforms_aa_substitutions.tsv
--

DROP TABLE IF EXISTS AlphaMissense_aa_substitutions;

CREATE TABLE AlphaMissense_aa_substitutions (
   uniprot_id VARCHAR(50) NOT NULL COMMENT 'UniProtKB accession number of the protein in which the variant induces a single amino-acid substitution',
   protein_variant VARCHAR(10) NOT NULL COMMENT 'AA change in format <RefAA>POS_aa<AltAA> (e.g. V2L)',
   am_pathogenicity FLOAT UNSIGNED COMMENT 'Calibrated AlphaMissense pathogenicity scores in interval [0.0,1.0]',
   am_class VARCHAR(50) COMMENT 'Text classifiation of the protein_variant into likely_benign, likely_pathogenic, or ambiguous',
  PRIMARY KEY (uniprot_id, protein_variant)
) ENGINE=InnoDB COMMENT 'Query by uniprot_id and protein_variant'
CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci;

DROP TABLE IF EXISTS AlphaMissense_isoforms_aa_substitutions;

CREATE TABLE AlphaMissense_isoforms_aa_substitutions (
   transcript_id VARCHAR(50) NOT NULL COMMENT 'ENSMBL transcript ID',
   protein_variant VARCHAR(10) NOT NULL COMMENT 'AA change in format <RefAA>POS_aa<AltAA> (e.g. V2L)',
   am_pathogenicity FLOAT UNSIGNED COMMENT 'Calibrated AlphaMissense pathogenicity scores in interval [0.0,1.0]',
   am_class VARCHAR(50) COMMENT 'Text classifiation of the protein_variant into likely_benign, likely_pathogenic, or ambiguous',
  PRIMARY KEY (transcript_id, protein_variant)
) ENGINE=InnoDB COMMENT 'Query by transcript_id and protein_variant'
CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci;
