CREATE TABLE IF NOT EXISTS Residue (
pdbid VARCHAR(50),
chain VARCHAR(10),
residue VARCHAR(10),
seqid INT,
ins INT,
rescode VARCHAR(1),
x DOUBLE,
y DOUBLE,
z DOUBLE,
PRIMARY KEY(pdbid,chain,seqid,ins),
KEY(pdbid,x,y,z)
)