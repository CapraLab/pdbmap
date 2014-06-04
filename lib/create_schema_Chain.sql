CREATE TABLE IF NOT EXISTS Chain (
label VARCHAR(100), # Dataset label
structid VARCHAR(50),
biounit INT,
model INT,
chain VARCHAR(10),
unp VARCHAR(50),
offset INT,
hybrid TINYINT,
sequence TEXT,
PRIMARY KEY(label,structid,biounit,model,chain),
KEY(label,unp)
)