-- After Creating the Table, populate it with the pdbmap/scripts/uniparc_parser.py
-- Importantly, the Uniparc identifiers are invariant.  Thus, you should not DROP the table
-- before running the script
--
-- NOTE THAT the fasta sequence can be quite long, and MEDIUMTEXT seems to work on MariaDB 5 and 10 alike
--      VARCHAR(60000) seems to cause all kinds of problems with JOINs by comparison

CREATE TABLE pdbmap_v14.Uniparc (
  uniparc VARCHAR(20) NOT NULL COMMENT 'Uniparc identifer in format UPI1234567890123',
  md5sum  VARCHAR(32) NOT NULL COMMENT '128 bit md5sum of the amino acid sequence',
  fasta MEDIUMTEXT NOT NULL COMMENT 'Amino acid sequence',
  PRIMARY KEY (uniparc),
  UNIQUE (md5sum)  -- The UNIQUE keyword will cause an INDEX to be create on the md5sum
  -- UNIQUE (fasta(100)) -- FOR NOW DO NOT INDEX ON the fasta string, because you lookup by md5sum - not AA string
) ENGINE=InnoDB COMMENT 'Map between UniParc identifiers, fasta sequences, and their md5sums'
CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci;

CREATE TABLE pdbmap_v14.UniparcHuman (
  uniparc VARCHAR(20) NOT NULL COMMENT 'Uniparc identifer in format UPI1234567890123',
  md5sum  VARCHAR(32) NOT NULL COMMENT '128 bit md5sum of the amino acid sequence',
  fasta MEDIUMTEXT NOT NULL COMMENT 'Amino acid sequence',
  PRIMARY KEY (uniparc),
  UNIQUE (md5sum)  -- The UNIQUE keyword will cause an INDEX to be create on the md5sum
  -- UNIQUE (fasta(100)) -- FOR NOW DO NOT INDEX ON the fasta string, because you lookup by md5sum - not AA string
) ENGINE=InnoDB COMMENT 'Map between UniParc identifiers, fasta sequences, and their md5sums'
CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci;
