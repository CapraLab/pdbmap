CREATE TABLE IF NOT EXISTS AlignmentScore (
label VARCHAR(100), # Dataset label
structid VARCHAR(50),
chain VARCHAR(10) BINARY,
transcript VARCHAR(50),
score REAL,
perc_aligned REAL,
perc_identity REAL,
alignment TEXT,
as_id BIGINT NOT NULL AUTO_INCREMENT, # Unique, direct-reference key
PRIMARY KEY(label,structid,chain,transcript),
KEY(as_id),
KEY(label,structid,chain),
KEY(label,transcript),
KEY(label,score),
KEY(label,perc_aligned),
KEY(label,perc_identity)
)
