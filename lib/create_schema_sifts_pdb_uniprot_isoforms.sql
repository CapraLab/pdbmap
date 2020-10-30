-- Given either a PDB ID or uniprot ID, this table returns allows reconstruction of the dictionary 
-- as returned by the source rest API documented: https://www.ebi.ac.uk/pdbe/api/doc/sifts.html
-- and as example would be: https://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/6cet
--
-- This mapping supplements the legacy xml, which is detailed by residue but lacks isoform
-- specific uniprot identifiers

CREATE TABLE IF NOT EXISTS sifts_mappings_pdb_uniprot_all_isoforms (
-- {
--  "6cet": {
--    "UniProt": {
--      "O75140": {
  pdbid varchar(20) NOT NULL       COMMENT '4 character rcsb pdb ID',
  uniprot_acc varchar(50) DEFAULT NULL COMMENT 'Ex: O12345-2  6 character uniprot ID and isoform specific id',

--        "identifier": "NPRL3_HUMAN", 
--        "name": "NPRL3_HUMAN", 
  identifier varchar(50) DEFAULT NULL  COMMENT 'id for the uniprot_acc.  Example: DEPD5_HUMAN',
  name varchar(100) DEFAULT NULL       COMMENT 'name associated with the uniprot_acc.  Example: DEPD5_HUMAN',

--        "mappings": [
  mapping_entity_id varchar(50) DEFAULT NULL     COMMENT 'molecule number in mmcif-speak, but rcsb says it need not be an integer',

--            "start": {
--              "author_residue_number": null, 
--              "author_insertion_code": "", 
--              "residue_number": 1
--            }, 
  mapping_start_author_residue_number int(11) DEFAULT NULL COMMENT 'First residue: last residue number in pdb addressing scheme', 
  mapping_start_author_insertion_code varchar(10) DEFAULT NULL COMMENT 'First residue: PDB-style residue insertion code - often NULL', 
  mapping_start_residue_number int(11) DEFAULT NULL COMMENT 'First residue: mmcif style residue index',


--            "end": {
--              "author_residue_number": null, 
--              "author_insertion_code": "", 
--              "residue_number": 380
--            }, 
  mapping_end_author_residue_number int(11) DEFAULT NULL COMMENT 'Last residue: last residue number in pdb addressing scheme', 
  mapping_end_author_insertion_code varchar(10) DEFAULT NULL COMMENT 'Last residue: PDB-style residue insertion code - often NULL', 
  mapping_end_residue_number int(11) DEFAULT NULL COMMENT 'Last residue: mmcif style residue index',

--            "chain_id": "N", 
--            "unp_start": 1, 
--            "unp_end": 380, 
--            "pdb_start": 1, 
--            "pdb_end": 380, 
--            "struct_asym_id": "A", 
--            "identity": 1
  mapping_pdb_chain varchar(20)       NOT NULL   COMMENT 'chain in pdb.  Ex A/B/Z etc',
  mapping_struct_asym_id varchar(20)  NOT NULL   COMMENT 'chain in mmcif.',

  mapping_pdb_start int(11) NOT NULL  COMMENT 'First residue: pdb number',
  mapping_pdb_start_insertion_code varchar(10) DEFAULT NULL  COMMENT 'First residue: pdb insertion code if seen',
  mapping_pdb_end int(11) NOT NULL  COMMENT 'Last residue: pdb residue number',
  mapping_pdb_end_insertion_code varchar(10) DEFAULT NULL  COMMENT 'Last residue: pdb insertion code if seen',

  mapping_unp_start int(11) NOT NULL  COMMENT 'First residue: index of uniprot residue',
  mapping_unp_end int(11) NOT NULL  COMMENT 'Last residue: index of uniprot residue',

  mapping_seq_identity float not null   COMMENT 'sequence identity calculated by ue: PDB-style residue insertion code - often NULL', 

  PRIMARY KEY (pdbid,uniprot_acc,mapping_pdb_chain,mapping_pdb_start),
  KEY (uniprot_acc,pdbid,mapping_unp_start,mapping_pdb_chain)
) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci


-- The best isoforms from sifts come in identical format
-- but it is sanity-saving to keep separate source resources
-- separated in our tables, rather than having an additional "best isoform" column
-- in the above table
CREATE TABLE IF NOT EXISTS sifts_mappings_pdb_uniprot_best_isoforms
  (PRIMARY KEY (pdbid,uniprot_acc,mapping_pdb_chain,mapping_pdb_start),
  KEY (uniprot_acc,pdbid,mapping_unp_start,mapping_pdb_chain))
   AS SELECT * FROM sifts_mappings_pdb_uniprot_all_isoforms where 1=2
