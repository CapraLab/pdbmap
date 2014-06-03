CREATE TABLE IF NOT EXISTS GenomicIntersection (
dlabel VARCHAR(100), # Dataset label
slabel VARCHAR(100), # Structure label
structid VARCHAR(50), # Structure ID
chain VARCHAR(10), # Structure chain
seqid INT, # Position in chain sequence
gc_id MEDIUMINT, # GenomicConsequence direct reference key
PRIMARY KEY(label,structid,chain,seqid,gc_id),
KEY(gc_id),
KEY(dlabel),
KEY(slabel)
)
