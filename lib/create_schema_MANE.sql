-- MANE cross-references link the human genome to transcript databases in a new and convincing manner
-- https://www.ncbi.nlm.nih.gov/refseq/MANE/
--
-- For the pipeline, the MANE cross-references are mined from the text format uniprot_sprot.dat.gz, 
-- or its human-only subset as created by
-- DOWNLOAD_uniprot.bash (see pipeline_download_scripts repository at github)
--
-- Below find the SQL command to create a pipeline table, and populate it from the downloaded uniprot file
--
-- As the VUStruct pipeline evolves, we expect to increasingly integrate these cross-references as preferred
-- over the uniprot assigned cross-references
-- =============================================================
-- After you run the "CREATE TABLE..." command below, continue with extraction and uploaded from
-- your downloaded uniprot files
--
-- $ grep -i 'MANE-Select' data/uniprot/current/uniprot_sprot_human.dat > /tmp/MANE_human.txt
-- The resulting output format is annoying, with a uniprot ID bracket delimited in some rows
-- and a terminating period after the refseq_protein_id
--
-- To convince yourself that the illegible awk below is correct, you can do
-- $ echo 'DR   MANE-Select; ENST00000614654.2; ENSP00000480314.1; NM_199341.4; NP_955373.3. [Q86XI8-2]' | \
-- awk -F' *[;\\[\\]] *'  '{print $2"\t"$3"\t"$4"\t"substr($5,1,length($5)-1)"\t"$6}' 
-- output: ENST00000614654.2	ENSP00000480314.1	NM_199341.4	NP_955373.3	Q86XI8-2
--
-- Proceed to globally convert the uniprot DR row format to 5-column tab delimited format for mariadb-import
-- $ awk -F' *[;\\[\\]] *'  '{print $2"\t"$3"\t"$4"\t"substr($5,1,length($5)-1)"\t"$6}' /tmp/MANE_human.txt > /tmp/MANE_xref
-- $ mariadb-import --ignore --verbose --fields-terminated-by='\t' --local -p  --columns=unp,ID_type,ID pdbmap_v15 /tmp/MANE_xref

CREATE TABLE MANE_xref (
  ensembl_transcript_id VARCHAR(40) NOT NULL COMMENT 'ENST000... Ensembl transcript ID',
  ensembl_protein_id VARCHAR(40) NOT NULL COMMENT 'ENSP000... Ensembl protein ID',
  refseq_transcript_id VARCHAR(40) NOT NULL COMMENT 'NM_nnnnn refseq ID that is cross-referenced',
  refseq_uniprot_protein_id VARCHAR(40) NOT NULL COMMENT 'NP_nnnnn refseq ID that is cross-referenced',
  optional_uniprot_id VARCHAR(40) NOT NULL COMMENT 'the bracketed uniprot ID that appears in some rows of raw data',
  PRIMARY KEY (ensembl_transcript_id)
) ENGINE=InnoDB COMMENT 'Lookup by ENST000 ID verifies MANE cross reference, and returns refseq IDs'

CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci;
