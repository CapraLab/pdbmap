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
KEY(structid,chain),
KEY(transcript),
KEY(score),
KEY(perc_aligned),
KEY(perc_identity)
)
