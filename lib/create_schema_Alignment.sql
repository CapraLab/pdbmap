CREATE TABLE IF NOT EXISTS Alignment (
label VARCHAR(100), # Dataset label
structid VARCHAR(50),
chain VARCHAR(10),
chain_seqid INT,
transcript VARCHAR(50),
trans_seqid INT,
PRIMARY KEY(structid,chain,chain_seqid,transcript,trans_seqid),
KEY(structid,chain),
KEY(transcript)
)