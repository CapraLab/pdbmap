CREATE TABLE IF NOT EXISTS GenomicIntersection (
pdbid VARCHAR(50), # Structure ID
chain VARCHAR(10), # Structure chain
seqid INT, # Position in chain sequence
gc_id MEDIUMINT, # GenomicConsequence direct reference key
