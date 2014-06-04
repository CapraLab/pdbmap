CREATE TABLE IF NOT EXISTS AlignmentScore (
label VARCHAR(100), # Dataset label
structid VARCHAR(50),
chain VARCHAR(10),
transcript VARCHAR(50),
score REAL,
perc_aligned REAL,
perc_identity REAL,
alignment TEXT,
PRIMARY KEY(label,structid,chain,transcript),
KEY(label,structid,chain),
KEY(label,transcript),
KEY(label,score),
KEY(label,perc_aligned),
KEY(label,perc_identity)
)
