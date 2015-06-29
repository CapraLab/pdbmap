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
identity DOUBLE,
str_id MEDIUMINT NOT NULL AUTO_INCREMENT, # Unique, direct-reference key
PRIMARY KEY(label,modelid),
KEY(str_id),
KEY(label,pdbid),
KEY(label,unp),
KEY(label,mpqs),
KEY(label,zdope),
KEY(label,identity)
)