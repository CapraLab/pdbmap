-- After Creating the Table, populate it with the pdbmap/scripts/uniparc_parser.py

CREATE TABLE pdbmap_v14.Uniparc (
  uniparc CHAR(13) NOT NULL COMMENT 'Uniparc identifer in format UPI1234567890123',
  md5sum  CHAR(32) NOT NULL COMMENT '128 bit md5sum of the amino acid sequence',
  fasta VARCHAR(60000) NOT NULL COMMENT 'Amino acid sequence',
  PRIMARY KEY (uniparc),
  UNIQUE (md5sum)
  -- UNIQUE (fasta(100)) -- FOR NOW DO NOT INDEX ON the fasta string, because you lookup by md5sum - not AA string
) ENGINE=InnoDB COMMENT 'Map between UniParc identifiers, fasta sequences, and their md5sums'
