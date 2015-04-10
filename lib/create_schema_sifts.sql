CREATE TABLE `sifts` (
  `pdbid` varchar(100) NOT NULL DEFAULT '',
  `chain` varchar(50) NOT NULL DEFAULT '',
  `sp` varchar(100) DEFAULT NULL,
  `pdb_seqid` int(11) NOT NULL DEFAULT '0',
  `sp_seqid` int(11) DEFAULT NULL,
  PRIMARY KEY (`pdbid`,`chain`,`pdb_seqid`),
  KEY `sp` (`sp`,`sp_seqid`)
);
