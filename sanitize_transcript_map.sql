DROP TABLE PDBSide;

CREATE TABLE PDBSide AS
(SELECT PDB.id,PDB.chain,UPT.transcript,PDB.aa3,PDB.aa1,PDB.seqres,PDB.x,PDB.y,PDB.z
FROM PDB as PDB,PDB_UniProt as PUP,UniProt_Transcript as UPT
WHERE PDB.id=PUP.pdb AND PUP.unp=UPT.unp);

DELETE FROM UniProt_Transcript
WHERE
# Delete if the amino acids don't match at the minimum seqres
UniProt_Transcript.transcript NOT IN 
	(SELECT transcript FROM Transcript_Genome as TG WHERE
		# PDB seqres must match the TG seqres (no out of range)
		TG.seqres=(SELECT seqres FROM PDBSide ORDER BY seqres LIMIT 1)
		AND
		# PDB amino acid must match the TG amino acid (no mismatched peptides)
		TG.aa1=(SELECT aa1 FROM PDBSide ORDER BY seqres LIMIT 1)
	)
OR
# Delete if the amino acids don't match at the maximum seqres
UniProt_Transcript.transcript NOT IN 
	(SELECT transcript FROM Transcript_Genome as TG WHERE
		# PDB seqres must match the TG seqres (no out of range)
		TG.seqres=(SELECT seqres FROM PDBSide ORDER BY seqres DESC LIMIT 1)
		AND
		# PDB amino acid must match the TG amino acid (no mismatched peptides)
		TG.aa1=(SELECT aa1 FROM PDBSide ORDER BY seqres DESC LIMIT 1)
	)