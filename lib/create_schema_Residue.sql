CREATE TABLE IF NOT EXISTS Residue (
pdbid VARCHAR(50),
chain VARCHAR(10),
resname VARCHAR(10),
rescode VARCHAR(1),
seqid INT,
ins INT,
x DOUBLE,
y DOUBLE,
z DOUBLE,
PRIMARY KEY(pdbid,chain,seqid,ins),
KEY(pdbid,x,y,z)
)