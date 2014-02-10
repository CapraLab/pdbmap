CREATE TABLE IF NOT EXISTS AlignmentScores (
pdbid VARCHAR(50),
chain VARCHAR(10),
transcript VARCHAR(50),
score REAL,
perc_nongap REAL,
perc_match REAL,
alignment TEXT,
PRIMARY KEY(pdbid,chain,transcript),
KEY(pdbid,chain),
KEY(transcript),
KEY(score),
KEY(perc_nongap),
KEY(perc_match)
)