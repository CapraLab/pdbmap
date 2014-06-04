CREATE TABLE IF NOT EXISTS Model (
label VARCHAR(100),
modelid VARCHAR(100),
unp VARCHAR(100),
method VARCHAR(100),
no35 DOUBLE,
rmsd DOUBLE,
mpqs DOUBLE,
evalue DOUBLE,
ga341 DOUBLE,
zdope DOUBLE,
pdbid VARCHAR(50),
chain VARCHAR(10),
PRIMARY KEY(label,modelid),
KEY(label,pdbid),
KEY(label,unp),
KEY(label,mpqs),
KEY(label,zdope)
)