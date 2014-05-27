CREATE TABLE IF NOT EXISTS GenomicIntersection (
label VARCHAR(100), # Dataset label
pdbid VARCHAR(50), # Structure ID
chain VARCHAR(10), # Structure chain
seqid INT, # Position in chain sequence
gc_id MEDIUMINT, # GenomicConsequence direct reference key
PRIMARY KEY(label,pdbid,chain,seqid,gc_id),
KEY(gc_id),
KEY(label)
)
