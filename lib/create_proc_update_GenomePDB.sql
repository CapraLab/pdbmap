CREATE DEFINER=`mike`@`gwar-dev.mc.vanderbilt.edu` PROCEDURE `update_GenomePDB`(new_pdbid VARCHAR(50))
BEGIN
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
WHERE d.pdbid=new_pdbid;
END
