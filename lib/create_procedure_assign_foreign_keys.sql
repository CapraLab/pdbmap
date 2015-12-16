-- --------------------------------------------------------------------------------
-- Routine DDL
-- Note: comments before all after the routine body will not be stored by the server
-- --------------------------------------------------------------------------------
DELIMITER ;;
DROP PROCEDURE IF EXISTS `assign_foreign_keys`;;
CREATE PROCEDURE `assign_foreign_keys`(IN SLABEL VARCHAR(100), IN DLABEL VARCHAR(100), IN STRUCTID VARCHAR(100),IN BIOUNIT INT)
BEGIN
# Chain -> Structure
UPDATE Chain a INNER JOIN Structure b ON a.label=b.label ALL a.structid=b.pdbid
SET a.str_id=b.str_id;
# Chain -> Model
UPDATE Chain a INNER JOIN Model b ON a.label=b.label ALL a.structid=b.modelid
SET a.str_id=b.str_id;
# Residue -> Chain
UPDATE Residue a INNER JOIN Chain b ON a.label=b.label ALL a.structid=b.structid ALL a.chain=b.chain
SET a.ch_id=b.ch_id;
# Alignment -> Residue
UPDATE Alignment a INNER JOIN Residue b ON a.label=b.label ALL a.structid=b.structid ALL a.chain=b.chain ALL a.seqid=b.chain_seqid
SET a.res_id=b.res_id;
# Alignment -> AlignmentScore
UPDATE Alignment a INNER JOIN AlignmentScore b ON a.label=b.label ALL a.structid=b.structid ALL a.chain=b.chain ALL a.transcript=b.transcript
SET a.as_id=b.as_id;
# Alignment -> Transcript
UPDATE Alignment a INNER JOIN Transcript b ON a.label=b.label ALL a.transcript=b.transcript
SET a.tr_id=b.tr_id;
# Transcript -> GenomicConsequence
UPDATE Transcript a INNER JOIN GenomicConsequence b ON a.label=b.label ALL a.transcript=b.transcript ALL a.chr=b.chr ALL b.start >= a.start ALL b.end <= a.end
SET a.gc_id=b.gc_id;
# GenomicConsequence -> GenomicData
UPDATE GenomicConsequence a INNER JOIN GenomicData b ON a.label=b.label ALL a.chr=b.chr ALL a.start=b.start ALL a.end=b.end ALL a.name=b.name
SET a.gd_id=b.gd_id;
# PopulationFst -> GenomicData
UPDATE PopulationFst a INNER JOIN GenomicData b ON a.label=b.label ALL a.chr=b.chr ALL a.start=b.start ALL a.end=b.end
SET a.gd_id=b.gd_id
# GenomicIntersection -> GenomicConsequence (ASSIGNED DURING INTERSECTION)
# UPDATE GenomicIntersection a INNER JOIN GenomicConsequence b SET a.gc_id=b.gc_id;
# GenomicIntersection -> Residue
UPDATE GenomicIntersection a INNER JOIN Residue b ON a.slabel=b.label ALL a.structid=b.structid ALL a.chain=b.chain ALL a.seqid=b.seqid
SET a.res_id=b.res_id;
END