DROP TABLE PDBSide;
DROP TABLE EnsemblSide;
DROP TABLE EnsemblPDB;

CREATE TABLE PDBSide AS
(SELECT PDB.id,PDB.chain,UPT.transcript,PDB.aa3,PDB.seqres,PDB.x,PDB.y,PDB.z
FROM PDB as PDB,PDB_UniProt as PUP,UniProt_Transcript as UPT
WHERE PDB.id=PUP.pdb AND PUP.unp=UPT.unp);

CREATE TABLE EnsemblSide AS
(SELECT TG.transcript,TG.seqres,TG.aa1,TG.start,TG.end,TG.strand,T.chr
FROM UniProt_Transcript as UPT,Transcript_Genome as TG,Transcript as T
WHERE UPT.transcript=TG.transcript AND TG.transcript=T.id);

CREATE TABLE EnsemblPDB AS
SELECT PDBSide.id as pdbid,PDBSide.chain,PDBSide.seqres,PDBSide.aa3,EnsemblSide.aa1,PDBSide.x,PDBSide.y,PDBSide.z,EnsemblSide.transcript,EnsemblSide.start,EnsemblSide.end,EnsemblSide.strand,EnsemblSide.chr
FROM PDBSide INNER JOIN EnsemblSide ON (PDBSide.transcript=EnsemblSide.transcript AND PDBSide.seqres=EnsemblSide.seqres)
GROUP BY PDBSide.chain,PDBSide.seqres
ORDER BY PDBSide.chain,PDBSide.seqres;

SELECT * FROM PDBSide;
SELECT * FROM EnsemblSide WHERE seqres > 33;
SELECT * FROM EnsemblPDB;