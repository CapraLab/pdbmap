CREATE DEFINER=`%(dbuser)s`@`%(dbhost)s` PROCEDURE `build_GenomePDB`()
BEGIN

DROP TABLE IF EXISTS GenomePDB;
DROP TABLE IF EXISTS PDB;

CREATE TABLE PDB AS
(SELECT PDBInfo.pdbid,PDBInfo.species,PDBInfo.chain,PDBInfo.unp,PDBCoords.seqres,PDBCoords.aa1,PDBCoords.x,PDBCoords.y,PDBCoords.z,PDBTranscript.transcript
FROM PDBInfo,PDBCoords,PDBTranscript
WHERE PDBInfo.pdbid=PDBCoords.pdbid AND PDBInfo.chain=PDBCoords.chain AND PDBInfo.pdbid=PDBTranscript.pdbid AND PDBInfo.chain=PDBTranscript.chain AND PDBCoords.pdbid=PDBTranscript.pdbid AND PDBCoords.chain=PDBTranscript.chain);

CREATE TABLE GenomePDB AS
SELECT PDB.pdbid,PDB.species,PDB.chain,PDB.unp,PDB.transcript,PDB.seqres,PDB.aa1,PDB.x,PDB.y,PDB.z,GenomicCoords.start,GenomicCoords.end,GenomicCoords.chr,GenomicCoords.strand
FROM PDB,GenomicCoords
WHERE PDB.transcript=GenomicCoords.transcript AND PDB.seqres=GenomicCoords.seqres
GROUP BY PDB.pdbid,PDB.chain,PDB.seqres
ORDER BY PDB.pdbid,PDB.chain,PDB.seqres;

DROP TABLE PDB;
END