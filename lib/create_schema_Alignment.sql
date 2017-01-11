CREATE TABLE IF NOT EXISTS Alignment (
label VARCHAR(100), # Dataset label
structid VARCHAR(50),
chain VARCHAR(10) BINARY,
chain_seqid INT,
transcript VARCHAR(50),
trans_seqid INT,
al_id BIGINT NOT NULL AUTO_INCREMENT, # Unique, direct-reference key
as_id BIGINT, # AlignmentScore direct-reference key
res_id BIGINT, # Residue direct-reference key
tr_id BIGINT, # Transcript direct-reference key
PRIMARY KEY(label,structid,chain,chain_seqid,transcript,trans_seqid),
KEY(al_id),
KEY(as_id),
KEY(res_id),
KEY(tr_id),
KEY(label,structid,chain,chain_seqid),
KEY(label,transcript,trans_seqid)
)