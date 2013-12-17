CREATE DEFINER=`mike`@`gwar-dev.mc.vanderbilt.edu` PROCEDURE `build_GenomePDB`()
BEGIN

DROP TABLE IF EXISTS GenomePDB;

CREATE TABLE GenomePDB 
(
chr VARCHAR(10),
start BIGINT,
end BIGINT,
name VARCHAR(20),
strand INT,
gene VARCHAR(20), 
transcript VARCHAR(20),
trans_seq INT,
trans_aa1 VARCHAR(1),
pdbid VARCHAR(20),
chain VARCHAR(1),
species VARCHAR(20),
unp VARCHAR(20),
chain_seq INT(11),
chain_aa1 VARCHAR(1),
x DOUBLE,
y DOUBLE,
z DOUBLE,
PRIMARY KEY `transpos` (`pdbid`,`chain`,`chain_seq`,`transcript`,`trans_seq`),
KEY `pdbid` (`pdbid`),
KEY `gene` (`gene`),
KEY `trans` (`transcript`),
KEY `peptide` (`pdbid`,`chain`,`chain_seq`),
KEY `name` (`name`),
KEY `genomic` (`chr`,`start`,`end`));

INSERT INTO GenomePDB
SELECT a.chr,a.start,a.end,a.strand,a.gene,a.transcript,b.trans_seq,a.aa1,d.pdbid,d.chain,d.species,d.unp,b.chain_seq,c.aa1,c.x,c.y,c.z
FROM GenomicCoords AS a
INNER JOIN
Alignment AS b
ON a.transcript=b.transcript AND a.seqres=b.trans_seq
INNER JOIN
PDBCoords AS c
ON b.pdbid=c.pdbid AND b.chain=c.chain AND b.chain_seq=c.seqres
INNER JOIN
PDBInfo as d
ON c.pdbid=d.pdbid AND c.chain=d.chain
ORDER BY pdbid,chain,chain_seq;

END