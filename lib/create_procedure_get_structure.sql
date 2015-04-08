-- --------------------------------------------------------------------------------
-- Query a single asymmetric unit, NMR, model, or biological assembly from PDBMap
-- --------------------------------------------------------------------------------
DELIMITER ;;
DROP PROCEDURE IF EXISTS`get_structure`;;
CREATE PROCEDURE `get_structure` (IN SLABEL VARCHAR(100), IN DLABEL VARCHAR(100), IN STRUCTID VARCHAR(100),IN BIOUNIT INT)
BEGIN
SELECT
/*SLABEL*/a.label as slabel,
/*ID*/a.structid,
/*METHOD*/IF(c.method IS NULL,d.method,c.method) as method,
/*STRUCTURE*/c.resolution,
/*MODEL*/d.mpqs,
/*CHAIN*/b.biounit,b.model,b.chain,b.unp,b.hybrid,
/*RESIDUE*/a.seqid,a.icode,a.rescode,a.ss,a.phi,a.psi,a.tco,a.kappa,a.alpha,a.x,a.y,a.z,
/*PFAM*/p.acc as pfamid,p.name as pfam_domain,p.description as pfam_desc,p.evalue as pfam_evalue,
/*ALIGN*/h.chain_seqid as aln_chain_seqid,h.trans_seqid as aln_trans_seqid,j.perc_identity,
/*TRANS*/i.transcript as ens_trans,i.gene as ens_gene,i.protein as ens_prot,i.seqid as ens_prot_seqid,i.rescode as ens_prot_aa,
/*VARIANT*/IF(f.consequence LIKE '%missense_variant%', 1, 0) as issnp,
/*DLABEL*/g.label as dlabel,
/*VARIANT*/g.name as snpid,g.chr,g.start,g.end,g.hgnc_gene,g.ens_gene,g.aa as anc_allele,g.ref_allele,g.alt_allele,g.maf,g.amr_af,g.eas_af,g.sas_af,g.eur_af,g.afr_af,g.allpop_Fst,
/*CONSEQUENCE*/f.gc_id,f.transcript as vep_trans,f.protein as vep_prot,f.protein_pos as vep_prot_pos,f.ref_codon,f.alt_codon,f.ref_amino_acid as vep_ref_aa,f.alt_amino_acid as vep_alt_aa,
/*CONSEQUENCE*/f.consequence,f.polyphen,f.sift,f.biotype
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