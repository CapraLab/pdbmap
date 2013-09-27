CREATE DEFINER=`mike`@`gwar-dev.mc.vanderbilt.edu` PROCEDURE `update_GenomePDB`(new_pdbid VARCHAR(20))
BEGIN
INSERT INTO GenomePDB
SELECT d.chr,d.start,d.end,d.strand,d.gene,d.transcript,c.trans_seq,d.aa1,a.pdbid,a.chain,a.species,a.unp,c.chain_seq,b.aa1,b.x,b.y,b.z
FROM PDBInfo as a
INNER JOIN
PDBCoords as b
ON a.pdbid=b.pdbid AND a.chain=b.chain
INNER JOIN
Alignment as c
ON b.pdbid=c.pdbid AND b.chain=c.chain AND b.seqres=c.chain_seq
INNER JOIN
GenomicCoords as d
ON c.transcript=d.transcript AND c.trans_seq=d.seqres
GROUP BY a.pdbid,a.chain,b.seqres;
END