-- Many different 3D chains can be mapped to many different transcripts
-- and vice-versa.  The two tables here capture all these mappings.
-- To retrieve all the residue to residue alignments, LEFT JOIN the
-- second TABLE to the first on column al_id

-- Many different 3D chains can be mapped to many different transcripts
-- and vice-versa.  The two tables here capture all these mappings.
-- To retrieve all the residue to residue alignments, LEFT JOIN the
-- second TABLE to the first on column al_id

CREATE TABLE AlignChainTranscript (
  al_id BIGINT NOT NULL AUTO_INCREMENT COMMENT 'Unique id and master key for AlignChainTranscriptResidues',
  label VARCHAR(20)        NOT NULL COMMENT 'Chain Source: pdb/modbase/swiss/ref90/etc',
  structid VARCHAR(50)     NOT NULL COMMENT 'Unique pdb ID, ENSP..., or swiss ID',
  chain VARCHAR(10)        NOT NULL COMMENT 'Typically one character, but unlimied',
  transcript VARCHAR(50)   NOT NULL COMMENT 'ENST.... Ensembl Transcript ID string',
  -- Columns from the old Alignment table, never populated
  -- as_id BIGINT, # AlignmentScore direct-reference key
  -- res_id BIGINT, # Residue direct-reference key
  -- tr_id BIGINT, # Transcript direct-reference key
  -- ALIGNMENT STATISTICS FOLLOW
  score REAL               COMMENT 'Score from pairwise2.aligned.   PDBMapAlignment.py',
  perc_aligned REAL        COMMENT 'Percentage of residues aligned. PDBMapAlignment.py',
  perc_identity REAL       COMMENT 'Percentage of exact matches.    PDBMapAlignment.py',
  alignment TEXT           COMMENT 'Alignment string. Format: ---ABC--XYZ--- etc',
  PRIMARY KEY (label,structid,chain,transcript),
  KEY (al_id),
  -- KEY(as_id),
  -- KEY(res_id),
  -- KEY(tr_id),
  KEY (transcript) 
) ENGINE=InnoDB COMMENT 'Chains<->Transcripts. JOIN to AlignChainTranscriptResidues';

CREATE TABLE AlignChainTranscriptResidue (
  chain_res_num INT NOT NULL COMMENT 'Residue number in 3D chain',
  chain_res_icode CHAR(1) NOT NULL DEFAULT ' ' COMMENT 'Space or occasional PDB insertion code. Example \'A\'',
  trans_seqid INT NOT NULL COMMENT 'Transcript sequence number',
  al_id BIGINT NOT NULL COMMENT 'Foreign key to join back to AlignChainTranscript',
  PRIMARY KEY(al_id,chain_res_num,chain_res_icode,trans_seqid),
  CONSTRAINT AlignChainTranscriptResidueRequiresAlignChainTranscript FOREIGN KEY(al_id) 
    REFERENCES AlignChainTranscript(al_id)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT UniqueChainResidueNumber UNIQUE INDEX ChainResidueIndex (al_id,chain_res_num,chain_res_icode),
  CONSTRAINT UniqueTranscriptResidue  UNIQUE INDEX TranscriptResidueIndex (al_id,trans_seqid)
) ENGINE=Innodb COMMENT 'Residue Detail Table for AlignChainTranscript';

