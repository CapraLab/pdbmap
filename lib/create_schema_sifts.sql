CREATE TABLE `sifts` (
  `pdbid` varchar(50) NOT NULL,
  `chain` varchar(50) NOT NULL,
  `resnum` int(11) NOT NULL,
  `icode` varchar(50) NOT NULL,
  `resname` varchar(50) DEFAULT NULL,
  `uniprot_acc` varchar(50) DEFAULT NULL,
  `uniprot_resnum` int(11) DEFAULT NULL,
  `uniprot_resname` varchar(50) DEFAULT NULL,
  `ncbi` varchar(50) DEFAULT NULL,
  `pfam` varchar(100) DEFAULT NULL,
  `cath` varchar(100) DEFAULT NULL,
  `scop` varchar(100) DEFAULT NULL,
  `interpro` varchar(200) DEFAULT NULL,
  `sscode` varchar(50) DEFAULT NULL,
  `ssname` varchar(50) DEFAULT NULL,
  `sft_id MEDIUMINT NOT NULL AUTO_INCREMENT, # Unique, direct-reference key
  PRIMARY KEY (`pdbid`,`chain`,`resnum`,`icode`),
  KEY (`sft_id`)
)

