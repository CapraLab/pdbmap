-- --------------------------------------------------------------------------------
-- This procedure returns all biological assemblies containing a given
-- UniProt protein, along with the relevant chain.
-- --------------------------------------------------------------------------------
DELIMITER ;;
DROP PROCEDURE IF EXISTS `get_protein`;;
CREATE PROCEDURE `get_protein` (IN SLABEL VARCHAR(100), IN UNP VARCHAR(100), IN MAXRES DOUBLE, IN MINALN DOUBLE)
BEGIN
select a.label,a.structid,a.biounit,b.resolution as res,c.perc_identity as aln,GROUP_CONCAT('#',a.model,'.',a.chain,';') as chains
from Chain a
inner join Structure b
on a.label=slabel and a.unp=unp
and a.label=b.label and a.structid=b.pdbid 
and b.resolution<maxres
inner join AlignmentScore c
on a.label=c.label and a.structid=c.structid and a.chain=c.chain
where c.perc_identity>minaln
group by a.label,a.structid,a.biounit
order by b.resolution asc,c.perc_identity desc;

select a.label,a.structid,a.biounit,b.mpqs,c.perc_identity as aln,GROUP_CONCAT('#',a.model,'.',a.chain,';') as chains
from Chain a
inner join Model b
on a.label=slabel and a.unp=unp
and a.label=b.label and a.structid=b.modelid 
and b.mpqs>(2-(log(maxres))/(2*log(2))) # arbitrary mapping of Resolution to MPQS (4A->~1, 8A->~0.5)
inner join AlignmentScore c
on a.label=c.label and a.structid=c.structid and a.chain=c.chain
where c.perc_identity>minaln
group by a.label,a.structid,a.biounit
order by b.mpqs desc,c.perc_identity desc;
END