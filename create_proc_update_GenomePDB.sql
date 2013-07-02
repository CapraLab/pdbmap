CREATE DEFINER=`mike`@`gwar-dev.mc.vanderbilt.edu` PROCEDURE `update_GenomePDB`(new_pdbid VARCHAR(20))
BEGIN
INSERT INTO GenomePDB
SELECT a.pdbid,a.species,a.chain,a.unp,c.transcript,b.seqres,b.aa1,b.x,b.y,b.z,d.start,d.end,d.chr,d.strand
FROM PDBInfo as a
INNER JOIN PDBCoords as b ON a.pdbid=b.pdbid AND a.chain=b.chain
INNER JOIN PDBTranscript as c ON a.pdbid=c.pdbid AND a.chain=c.chain
INNER JOIN GenomicCoords as d ON c.transcript=d.transcript AND b.seqres=d.seqres
WHERE a.pdbid=new_pdbid
GROUP BY a.pdbid,a.chain,b.seqres
ORDER BY a.pdbid,a.chain,b.seqres;
END