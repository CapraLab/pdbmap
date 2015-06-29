-- --------------------------------------------------------------------------------
-- Routine DDL
-- Note: comments before and after the routine body will not be stored by the server
-- --------------------------------------------------------------------------------
DELIMITER $$

CREATE PROCEDURE `pdbmap_stats` (IN SLABEL VARCHAR(200), IN DLABEL VARCHAR(200))
BEGIN
# Number of PDB structures
select count(distinct pdbid) as NumPDB from Structure where label=@SLABEL;

# Number of PDB Biological Assemblies
select count(distinct structid,biounit) as NumBio from Chain
where label=@SLABEL and biounit>0;

# Number of ModBase models
select count(distinct modelid) as NumModel from Model where label=@SLABEL;

# Number of proteins covered by biological assemblies/models?
select count(distinct unp) as NumProt from Chain
where label=@SLABEL;

# Number of proteins covered by the PDB
select count(distinct a.unp) as NumProt_PDB from Chain a
left join Structure b
on a.label=b.label and a.structid=b.pdbid
where a.label=@SLABEL and b.pdbid is not null;

# Number of proteins covered by Modbase
select count(distinct a.unp) as NumProt_ModBase from Chain a
left join Model b
on a.label=b.label and a.structid=b.modelid
where a.label=@SLABEL and b.modelid is not null;

# Number of nsSNPs
select count(distinct name) as NumSNPs from GenomicData 
where label=@DLABEL and abs(end-start)=1;

# Number of nsSNPs covered by biological assemblies/models
select count(distinct a.name) as NumSNPs_Mapped from GenomicData a
inner join GenomicConsequence b
on a.label=b.label and a.chr=b.chr and a.start=b.start and a.end=b.end and a.name=b.name
inner join GenomicIntersection c
on b.label=c.dlabel and b.gc_id=c.gc_id
where c.slabel=@SLABEL;

# Number of nsSNPs covered by the PDB
select count(distinct a.name) as NumSNPs_PDB from GenomicData a
inner join GenomicConsequence b
on a.label=b.label and a.chr=b.chr and a.start=b.start and a.end=b.end and a.name=b.name
inner join GenomicIntersection c
on b.label=c.dlabel and b.gc_id=c.gc_id
where c.slabel=@SLABEL and c.structid not like "ENSP%";

# Number of nsSNPs covered by ModBase
select count(distinct a.name) as NumSNPs_ModBase from GenomicData a
inner join GenomicConsequence b
on a.label=b.label and a.chr=b.chr and a.start=b.start and a.end=b.end and a.name=b.name
inner join GenomicIntersection c
on b.label=c.dlabel and b.gc_id=c.gc_id
where c.slabel=@SLABEL and c.structid like "ENSP%";
END