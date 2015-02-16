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
ch_id MEDIUMINT NOT NULL AUTO_INCREMENT, # Unique, direct-reference key
str_id MEDIUMINT, # Structure/Model direct-reference key
PRIMARY KEY(label,structid,biounit,model,chain),
KEY(ch_id),
KEY(str_id),
KEY(label,unp)
)