select a.structid,a.biounit,a.model,a.chain,a.seqid,a.icode,GROUP_CONCAT(y.acc) as pfam_acc,GROUP_CONCAT(y.name) as pfam_domain,a.x,a.y,a.z,d.chr,d.start,d.name,c.consequence as csq,d.aa,d.ref_allele as ref,d.%s_af as maf,gt
from Residue a
inner join pdbmap_supp.unp_repr_pdb e
on a.label=e.label and a.structid=e.pdb and a.biounit=e.biounit
inner join Chain x
on a.label=x.label and a.structid=x.structid and a.biounit=x.biounit and a.model=x.model and a.chain=x.chain
left join pfam y
on x.structid=y.pdbid and x.unp=y.unp and a.seqid between y.seqstart and y.seqend
left join GenomicIntersection b
on a.label=b.slabel and a.structid=b.structid and a.chain=b.chain and a.seqid=b.seqid and b.dlabel='1kg3'
left join GenomicConsequence c
on b.dlabel=c.label and b.gc_id=c.gc_id and c.consequence like '%%missense_variant%%'
left join GenomicData d
on c.label=d.label and c.chr=d.chr and c.start=d.start and c.end=d.end and c.name=d.name and c.end-c.start=1
where a.label='uniprot-pdb'
group by a.structid,a.biounit,a.model,a.chain,a.seqid
order by a.structid,a.biounit,a.model,a.chain,a.seqid;
