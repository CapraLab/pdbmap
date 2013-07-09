CREATE DEFINER=`%(dbuser)s`@`%(dbhost)s` PROCEDURE `sanitize_transcripts`(new_pdbid VARCHAR(20))
BEGIN

DROP TABLE IF EXISTS minTrans;
DROP TABLE IF EXISTS maxTrans;

# Genomic transcripts that match the minimum seqres for their PDB Sequence
CREATE TABLE minTrans AS (SELECT PDBTranscript.transcript FROM PDBTranscript,GenomicCoords as GC WHERE
		PDBTranscript.pdbid=new_pdbid AND
		# PDB seqres must match the GC seqres (no out of range)
		GC.seqres=(SELECT seqres FROM PDBCoords WHERE PDBCoords.pdbid=new_pdbid AND PDBCoords.pdbid=PDBTranscript.pdbid AND PDBCoords.chain=PDBTranscript.chain ORDER BY seqres LIMIT 1)
		AND
		# PDB amino acid must match the GC amino acid (no mismatched peptides)
		GC.aa1=(SELECT aa1 FROM PDBCoords WHERE PDBCoords.pdbid=new_pdbid AND PDBCoords.pdbid=PDBTranscript.pdbid AND PDBCoords.chain=PDBTranscript.chain ORDER BY seqres LIMIT 1)
	);

# Genomic transcripts that match the maximum seqres for their PDB Sequence
CREATE TABLE maxTrans AS (SELECT PDBTranscript.transcript FROM PDBTranscript,GenomicCoords as GC WHERE
		PDBTranscript.pdbid=new_pdbid AND
		# PDB seqres must match the GC seqres (no out of range)
		GC.seqres=(SELECT seqres FROM PDBCoords WHERE PDBCoords.pdbid=new_pdbid AND PDBCoords.pdbid=PDBTranscript.pdbid AND PDBCoords.chain=PDBTranscript.chain ORDER BY seqres DESC LIMIT 1)
		AND
		# PDB amino acid must match the GC amino acid (no mismatched peptides)
		GC.aa1=(SELECT aa1 FROM PDBCoords WHERE PDBCoords.pdbid=new_pdbid AND PDBCoords.pdbid=PDBTranscript.pdbid AND PDBCoords.chain=PDBTranscript.chain ORDER BY seqres DESC LIMIT 1)
	);

DELETE FROM PDBTranscript
WHERE
pdbid=new_pdbid
AND
(# Delete if the amino acids don't match at the minimum seqres
transcript NOT IN (SELECT transcript FROM minTrans)
OR
# Delete if the amino acids don't match at the maximum seqres
transcript NOT IN (SELECT transcript FROM maxTrans)
);
DROP TABLE minTrans;
DROP TABLE maxTrans;
END
