CREATE TABLE IF NOT EXISTS Residue (
label VARCHAR(100), # Dataset label
structid VARCHAR(50),
chain VARCHAR(10),
resname VARCHAR(10),
rescode VARCHAR(1),
seqid INT,
icode VARCHAR(10),
x DOUBLE,
y DOUBLE,
z DOUBLE,
PRIMARY KEY(structid,chain,seqid,icode),
KEY(structid,x,y,z)
)
