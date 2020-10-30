-- The idmapping input file is prepared via download and post-processing as
-- described in the data/uniprot/current/ subdirectory.
-- gunzip the file HUMAN_9606_idmapping_sprot.dat.gz to /tmp/Idmapping
--
-- By loading this file into an indexed table, startup time is dramatically reduced
-- ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
-- explains that the three columns are unp, ID_type, and ID
--
-- Tab delimiter allows us to load the idmapping file with the command line, easily
-- for example: mysqlimport mysqlimport --ignore --verbose --fields-terminated-by='\t' --local -p  --columns=unp,ID_type,ID pdbmap_v14 /tmp/Idmapping

CREATE TABLE Idmapping (
  unp VARCHAR(40) NOT NULL COMMENT 'UniProtKB-AC  Typially A12345-nn max, but occasionally longer, especially if not curated',
  ID_type VARCHAR(40) NOT NULL COMMENT 'the type of the right column ID that is being cross-referenced',
  ID VARCHAR(40) NOT NULL COMMENT 'the right column ID that is being cross-referenced',
  PRIMARY KEY (unp,ID_type,ID),
  KEY (ID,ID_type,unp), -- Speeds searches like "Find unps for an ENST sequence"
  KEY (ID_type,unp,ID)  -- Speeds searches like "retrieve all pdb ids"
) ENGINE=InnoDB COMMENT 'Bidirectionally maps uniprot IDs to other IDs'

CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci
