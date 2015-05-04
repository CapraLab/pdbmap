CREATE TABLE IF NOT EXISTS `pfam` (
  `pdbid` varchar(100) NOT NULL DEFAULT '',
  `chain` varchar(50) NOT NULL DEFAULT '',
  `seqstart` int(11) NOT NULL DEFAULT '0',
  `seqend` int(11) NOT NULL DEFAULT '0',
  `acc` varchar(50) NOT NULL DEFAULT '',
  `name` varchar(50) DEFAULT NULL,
  `description` text,
  `evalue` double DEFAULT NULL,
  PRIMARY KEY (`pdbid`,`chain`,`seqstart`,`seqend`,`acc`),
  KEY `acc` (`acc`),
  KEY `name` (`name`),
  KEY `evalue` (`evalue`)
  );
