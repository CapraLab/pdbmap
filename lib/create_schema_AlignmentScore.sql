CREATE TABLE IF NOT EXISTS AlignmentScores (
pdbid VARCHAR(50),
chain VARCHAR(10),
transcript VARCHAR(50),
score REAL,
perc_aligned REAL,
perc_identity REAL,
PRIMARY KEY(pdbid,chain,transcript),
KEY(pdbid,chain),
KEY(transcript),
KEY(score),
KEY(perc_aligned),
KEY(perc_identity)
)
