-- --------------------------------------------------------------------------------
-- Routine DDL
-- Note: comments before and after the routine body will not be stored by the server
-- --------------------------------------------------------------------------------
DELIMITER ;;
DROP PROCEDURE IF EXISTS `get_full_structure`;;
CREATE PROCEDURE `get_full_structure`(IN SLABEL VARCHAR(100), IN DLABEL VARCHAR(100), IN STRUCTID VARCHAR(100),IN BIOUNIT INT)
BEGIN
SELECT
/*LABELS*/g.label as dlabel,a.label as slabel,
/*STRUCTURE*/c.*,
/*MODEL*/d.*,
/*CHAIN*/b.*,
/*RESIDUE*/d.*,
/*VARIANT*/g.*,
/*CONSEQUENCE*/f.*,
/*ALIGN*/h.*,j.*,
/*TRANS*/i.*
FROM Residue as a
LEFT JOIN Chain as b 
ON a.label=b.label AND a.structid=b.structid AND a.biounit=b.biounit AND a.model=b.model AND a.chain=b.chain
LEFT JOIN Structure as c
ON b.label=c.label AND b.structid=c.pdbid
LEFT JOIN Model as d
ON b.label=d.label AND b.structid=d.modelid
LEFT JOIN GenomicIntersection as e
ON a.label=e.slabel AND a.structid=e.structid AND a.chain=e.chain AND a.seqid=e.seqid AND e.dlabel='1kg3'
LEFT JOIN GenomicConsequence as f
ON e.dlabel=f.label AND e.gc_id=f.gc_id AND f.canonical=1
LEFT JOIN GenomicData as g
ON f.label=g.label AND f.chr=g.chr AND f.start=g.start AND f.end=g.end AND f.name=g.name
LEFT JOIN Alignment as h USE INDEX(PRIMARY)
ON a.label=h.label AND a.structid=h.structid AND a.chain=h.chain AND a.seqid=h.chain_seqid AND f.transcript=h.transcript
LEFT JOIN Transcript as i USE INDEX(PRIMARY)
ON h.label=i.label AND h.transcript=i.transcript AND h.trans_seqid=i.seqid
LEFT JOIN AlignmentScore as j
ON h.label=j.label AND h.structid=j.structid AND h.chain=j.chain AND h.transcript=j.transcript
LEFT JOIN pfam as p
ON a.structid=p.pdbid AND a.chain=p.chain AND a.seqid BETWEEN p.seqstart AND p.seqend
AND (f.transcript IS NULL OR f.transcript=h.transcript)
where a.label=SLABEL and g.label=DLABEL
and a.structid=STRUCTID AND a.biounit=BIOUNIT
ORDER BY a.structid,a.biounit,a.model,a.chain,a.seqid;
END