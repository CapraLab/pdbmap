-- --------------------------------------------------------------------------------
-- This procedure returns all biological assemblies containing a given
-- UniProt protein, along with the relevant chain.
-- --------------------------------------------------------------------------------
DELIMITER ;;

CREATE PROCEDURE `get_protein` (IN SLABEL VARCHAR(100), IN UNP VARCHAR(100), IN MAXRES DOUBLE, IN MINALN DOUBLE)
BEGIN
select a.label,a.structid,a.biounit,GROUP_CONCAT('#',a.model,'.',a.chain,';') as chains
from Chain a
inner join Structure b
on a.label='uniprot-pdb' and a.unp='Q8IXJ6'
and a.label=b.label and a.structid=b.pdbid 
and b.resolution<4
inner join AlignmentScore c
on a.label=c.label and a.structid=c.structid and a.chain=c.chain
where c.perc_identity>0.99
group by a.label,a.structid,a.biounit;
END