CREATE TABLE IF NOT EXISTS Chain (
pdbid VARCHAR(50),
chain VARCHAR(10),
unp VARCHAR(50),
sequence TEXT,
PRIMARY KEY(pdbid,chain),
KEY(unp)
)