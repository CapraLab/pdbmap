-- --------------------------------------------------------------------------------
-- This procedure queries PDBMap for the best representative biological assembly
-- for each UniProt protein in PDBMap. Biological assemblies must be experimentally
-- determined, have a resolution no higher than 4A, have an alignment identity
-- >90% and contain the longest single chain of that protein.
-- --------------------------------------------------------------------------------
DELIMITER ;;

CREATE PROCEDURE `repr_subset` (IN SLABEL VARCHAR(100), IN MAXRES DOUBLE, IN MINALIGN DOUBLE)
BEGIN
# Report the UniProt ID and its representative structure
select y.label,y.unp,y.structid,y.biounit from
# For each UniProt ID, report the highest count observed in any chain
(select label,unp,max(numres) as numres from
(select b.label,b.structid,b.biounit,b.unp,LENGTH(b.sequence) as numres
from Chain b
inner join Structure c
on b.label=SLABEL
and b.label=c.label and b.structid=c.pdbid
and c.resolution<MAXRES and b.biounit>0
inner join AlignmentScore d
on b.label=d.label and b.structid=d.structid and b.chain=d.chain
and d.perc_identity>MINALIGN
group by b.label,b.structid,b.biounit,b.model,b.chain) z
group by label,unp) x
inner join
# Identify which biological assembly containing the chain with the highest count
(select b.label,b.structid,b.biounit,b.unp,LENGTH(b.sequence) as numres
from Chain b
inner join Structure c
on b.label=SLABEL
and b.label=c.label and b.structid=c.pdbid
and c.resolution<MAXRES and b.biounit>0
inner join AlignmentScore d
on b.label=d.label and b.structid=d.structid and b.chain=d.chain
and d.perc_identity>MINALIGN
group by b.label,b.structid,b.biounit,b.model,b.chain) y
on x.label=y.label and x.unp=y.unp and x.numres=y.numres
group by y.label,y.unp;
END;;