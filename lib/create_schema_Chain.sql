CREATE TABLE IF NOT EXISTS Chain (
label VARCHAR(100), # Dataset label
pdbid VARCHAR(50),
chain VARCHAR(10),
unp VARCHAR(50),
offset INT,
sequence TEXT,
PRIMARY KEY(pdbid,chain),
KEY(unp)
)