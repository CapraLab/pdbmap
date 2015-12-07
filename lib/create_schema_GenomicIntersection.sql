CREATE TABLE IF NOT EXISTS GenomicIntersection (
dlabel VARCHAR(100), # Dataset label
slabel VARCHAR(100), # Structure label
structid VARCHAR(50), # Structure ID
chain VARCHAR(10), # Structure chain
seqid INT, # Position in chain sequence
gi_id BIGINT NOT NULL AUTO_INCREMENT, # Unique, direct-reference key
gc_id BIGINT, # GenomicConsequence direct reference key
res_id BIGINT, # Residue direct reference key
PRIMARY KEY(slabel,dlabel,structid,chain,seqid,gc_id),
KEY(slabel,structid,chain,seqid,gc_id),
KEY(gi_id),
KEY(res_id),
KEY(gc_id)
)