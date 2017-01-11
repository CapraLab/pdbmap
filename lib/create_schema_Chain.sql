CREATE TABLE IF NOT EXISTS Chain (
label VARCHAR(100), # Dataset label
structid VARCHAR(50),
biounit INT,
model INT,
chain VARCHAR(10) BINARY,
unp VARCHAR(50),
gene VARCHAR(50),
offset INT,
hybrid TINYINT,
sequence TEXT,
ch_id BIGINT NOT NULL AUTO_INCREMENT, # Unique, direct-reference key
str_id BIGINT, # Structure/Model direct-reference key
PRIMARY KEY(label,structid,biounit,model,chain),
KEY(ch_id),
KEY(str_id),
KEY(label,unp)
)