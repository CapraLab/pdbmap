CREATE TABLE IF NOT EXISTS GenomicIntersection (
dlabel VARCHAR(100), # Dataset label
slabel VARCHAR(100), # Structure label
structid VARCHAR(50), # Structure ID
chain VARCHAR(10) BINARY, # Structure chain
seqid INT, # Position in chain sequence
gc_id BIGINT,  # GenomicConsequence direct reference key
gd_id BIGINT,  # GenomicData direct reference key
res_id BIGINT, # Residue direct reference key
gi_id BIGINT NOT NULL AUTO_INCREMENT, # Unique, direct-reference key
PRIMARY KEY(slabel,dlabel,structid,chain,seqid,gc_id),
KEY(slabel,structid,chain,seqid,gc_id),
KEY(gi_id),
KEY(res_id),
KEY(gc_id),
KEY(gd_id)
)