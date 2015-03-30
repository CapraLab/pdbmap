-- --------------------------------------------------------------------------------
-- This procedure queries PDBMap for the best representative biological assembly
-- for each UniProt protein in PDBMap. Biological assemblies must be experimentally
-- determined, have a resolution no higher than 4A, have an alignment identity
-- >90% and contain the longest single chain of that protein.
-- --------------------------------------------------------------------------------
DELIMITER ;;

CREATE PROCEDURE `repr_subset` (IN MAXRES DOUBLE, IN MINALIGN DOUBLE)
BEGIN
# Report the UniProt ID and its representative structure
select y.label,y.unp,y.structid,y.biounit from
# For each UniProt ID, report the highest count observed in any chain
(select label,unp,max(numres) as numres from
(select a.label,b.unp,count(distinct a.seqid) as numres
from Residue a
inner join Chain b
on a.label=b.label and a.structid=b.structid 
and a.biounit=b.biounit and a.model=b.model 
and a.chain=b.chain and a.biounit>0
inner join Structure c
on b.label=c.label and b.structid=c.pdbid
and c.resolution<MAXRES
inner join AlignmentScore d
on b.label=d.label and b.structid=d.structid and b.chain=d.chain
and d.perc_identity>MINALIGN
group by a.label,a.structid,a.biounit,a.model,a.chain) z
group by label,unp) x
inner join
# Identify which biological assembly containing the chain with the highest count
(select a.label,a.structid,a.biounit,b.unp,count(distinct a.seqid) as numres
from Residue a
inner join Chain b
on a.label=b.label and a.structid=b.structid 
and a.biounit=b.biounit and a.model=b.model 
and a.chain=b.chain and a.biounit>0
inner join Structure c
on b.label=c.label and b.structid=c.pdbid
and c.resolution<MAXRES
inner join AlignmentScore d
on b.label=d.label and b.structid=d.structid and b.chain=d.chain
and d.perc_identity>MINALIGN
group by a.label,a.structid,a.biounit,a.model,a.chain) y
on x.label=y.label and x.unp=y.unp and x.numres=y.numres
group by y.label,y.unp;
END;;