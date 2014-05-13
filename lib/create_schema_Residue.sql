CREATE TABLE IF NOT EXISTS Residue (
label VARCHAR(100), # Dataset label
structid VARCHAR(50),
biounit INT,
model INT,
chain VARCHAR(10),
resname VARCHAR(10),
rescode VARCHAR(1),
seqid INT,
icode VARCHAR(10),
x DOUBLE,
y DOUBLE,
z DOUBLE,
PRIMARY KEY(label,structid,biounit,model,chain,seqid,icode),
KEY(structid,biounit,x,y,z)
)
