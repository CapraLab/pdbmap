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
ss VARCHAR(1),
rsa DOUBLE,
phi DOUBLE,
psi DOUBLE,
tco DOUBLE,
kappa DOUBLE,
alpha DOUBLE,
res_id BIGINT NOT NULL AUTO_INCREMENT, # Unique, direct-reference key
ch_id BIGINT, # Chain direct-reference key
PRIMARY KEY(label,structid,biounit,model,chain,seqid,icode),
KEY(res_id),
KEY(ch_id),
KEY(label,structid,biounit,x,y,z),
KEY(label,x,y,z),
KEY(label,ss),
KEY(label,rsa),
KEY(label,phi,psi)
)
