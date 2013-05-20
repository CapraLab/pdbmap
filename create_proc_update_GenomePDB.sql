CREATE DEFINER=`%(dbuser)s`@`%(dbhost)s` PROCEDURE `update_GenomePDB`(new_pdbid VARCHAR(20))
BEGIN
CREATE TABLE IF NOT EXISTS GenomePDB (pdbid VARCHAR(20),species VARCHAR(20),chain VARCHAR(1),unp VARCHAR(20),transcript VARCHAR(20),seqres INT,aa1 VARCHAR(1),x DOUBLE,y DOUBLE,z DOUBLE,start BIGINT,end BIGINT,chr INT,strand INT, PRIMARY KEY(pdbid,chain,seqres));
DROP TABLE IF EXISTS PDB;

CREATE TABLE PDB AS
(SELECT PDBInfo.pdbid,PDBInfo.species,PDBInfo.chain,PDBInfo.unp,PDBCoords.seqres,PDBCoords.aa1,PDBCoords.x,PDBCoords.y,PDBCoords.z,PDBTranscript.transcript
FROM PDBInfo,PDBCoords,PDBTranscript
WHERE PDBInfo.pdbid=new_pdbid AND PDBCoords.pdbid=new_pdbid AND PDBTranscript.pdbid=new_pdbid 
AND PDBInfo.pdbid=PDBCoords.pdbid AND PDBInfo.chain=PDBCoords.chain 
AND PDBInfo.pdbid=PDBTranscript.pdbid AND PDBInfo.chain=PDBTranscript.chain 
AND PDBCoords.pdbid=PDBTranscript.pdbid AND PDBCoords.chain=PDBTranscript.chain);

INSERT INTO GenomePDB
SELECT PDB.pdbid,PDB.species,PDB.chain,PDB.unp,PDB.transcript,PDB.seqres,PDB.aa1,PDB.x,PDB.y,PDB.z,GenomicCoords.start,GenomicCoords.end,GenomicCoords.chr,GenomicCoords.strand
FROM PDB,GenomicCoords
WHERE PDB.pdbid=new_pdbid AND PDB.transcript=GenomicCoords.transcript AND PDB.seqres=GenomicCoords.seqres
GROUP BY PDB.pdbid,PDB.chain,PDB.seqres
ORDER BY PDB.pdbid,PDB.chain,PDB.seqres;

DROP TABLE PDB;
END