-- --------------------------------------------------------------------------------
-- Routine DDL
-- Note: comments before and after the routine body will not be stored by the server
-- --------------------------------------------------------------------------------
DELIMITER ;;
DROP PROCEDURE IF EXISTS `assign_foreign_keys`;;
CREATE PROCEDURE `assign_foreign_keys`(IN SLABEL VARCHAR(100), IN DLABEL VARCHAR(100), IN STRUCTID VARCHAR(100),IN BIOUNIT INT)
BEGIN
# Chain -> Structure
UPDATE Chain a INNER JOIN Structure b SET a.str_id=b.str_id;
# Chain -> Model
UPDATE Chain a INNER JOIN Model b SET a.str_id=b.str_id;
# Residue -> Chain
UPDATE Residue a INNER JOIN Chain b SET a.ch_id=b.ch_id;
# Alignment -> Residue
UPDATE Alignment a INNER JOIN Residue b SET a.res_id=b.res_id;
# Alignment -> AlignmentScore
UPDATE Alignment a INNER JOIN AlignmentScore b SET a.as_id=b.as_id;
# Alignment -> Transcript
UPDATE Alignment a INNER JOIN Transcript b SET a.tr_id=b.tr_id;
# Transcript -> GenomicConsequence
UPDATE Transcript a INNER JOIN GenomicConsequence b SET a.gc_id=b.gc_id;
# GenomicConsequence -> GenomicData
UPDATE GenomicConsequence a INNER JOIN GenomicData b SET a.gd_id=b.gd_id;
# GenomicIntersection -> GenomicConsequence
UPDATE GenomicIntersection a INNER JOIN GenomicConsequence b SET a.gc_id=b.gc_id;
# GenomicIntersection -> Residue
UPDATE GenomicIntersection a INNER JOIN GenomicConsequence b SET a.res_id=b.res_id;
END