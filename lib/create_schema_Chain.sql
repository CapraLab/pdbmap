CREATE TABLE IF NOT EXISTS Chain (
label VARCHAR(100), # Dataset label
structid VARCHAR(50),
chain VARCHAR(10),
unp VARCHAR(50),
offset INT,
hybrid TINYINT,
sequence TEXT,
PRIMARY KEY(structid,chain),
KEY(unp)
)