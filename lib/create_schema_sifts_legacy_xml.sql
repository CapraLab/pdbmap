-- Given a PDB ID and chain, and a canonical uniprot_acc, this table provides the SIFTS-supplied
-- sequence alignments between PDB residues (including insertion codes) and transcript
-- position numbers
--
-- PDB files are full of all kinds of troubles, from misisng residues, to inserted residues,
-- to residues in descending order.  The detailed SIFTS xml file mapping is the essential
-- glue that the pdbmap needs to effectively mate pdb 3D structures to the often purer
-- transcript sequences which those structures depict
-- 
-- The above reasons justify prefering curated SIFTS alignments to on-the-fly alignment
-- algorithsm.
-- 
-- Not all PDB residues will be aligned to transcript positions. 
-- Often, transcript positions will not align to PDB residues
--
-- This alignment does not involve isoform specific uniprot identifiers, and care must be taken
-- Both 1) to ensure that the uniprot amino acid sequence that has been aligned mathces
-- your sequence inquestion and 2) that this table is not employed to align non-canonical
-- uniprot identifers to pdb structures.  Those alignments are faciliated by other sifts*
-- tables generated from the Sifts' Rest API

CREATE TABLE IF NOT EXISTS sifts_legacy_xml (
  pdbid varchar(20) NOT NULL       COMMENT '4 character rcsb pdb ID',
  pdb_chain varchar(20)  NOT NULL  COMMENT 'chain.  Ex A/B/Z/1/d etc',
  PDBe_dbResNum int(11) NOT NULL   COMMENT 'The unique PDBe sourced incrementing dbResNum for indexing',
  pdb_resnum int(11) DEFAULT NULL  COMMENT 'pdb residue number, null if missing from .pdb',
  pdb_icode varchar(10) NOT NULL   COMMENT 'pdb insertion code, usually blank - sometimes single Alpha',
  pdb_resname varchar(20) DEFAULT NULL COMMENT 'pdb residue name',
  uniprot_acc varchar(50) DEFAULT NULL COMMENT 'Ex: O12345-2  6 character uniprot ID and isoform specific id',
  uniprot_resnum int(11) DEFAULT NULL  COMMENT 'The 1-starting index into the isoform-specific transcript',
  uniprot_resname varchar(50) DEFAULT NULL COMMENT 'transcript residue.  May vary from pdb',
  ncbi varchar(50) DEFAULT NULL,
  pfam varchar(100) DEFAULT NULL,
  cath varchar(100) DEFAULT NULL,
  scop varchar(100) DEFAULT NULL,
  interpro varchar(200) DEFAULT NULL,
  sscode varchar(50) DEFAULT NULL,
  ssname varchar(50) DEFAULT NULL,
  sft_id BIGINT NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (pdbid,pdb_chain,PDBe_dbResNum),
  KEY (pdbid,pdb_chain,pdb_resnum,pdb_icode,uniprot_acc),
  KEY (sft_id),
  KEY (uniprot_acc,pdb_resnum,pdbid) -- This is NOT a unique key.  One uniprot ID can be represented in many pdbs
)
CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci

