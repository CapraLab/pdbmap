CREATE TABLE IF NOT EXISTS Alignment (
label VARCHAR(100), # Dataset label
structid VARCHAR(50),
chain VARCHAR(10),
chain_seqid INT,
transcript VARCHAR(50),
trans_seqid INT,
al_id MEDIUMINT NOT NULL AUTO_INCREMENT, # Unique, direct-reference key
as_id MEDIUMINT, # AlignmentScore direct-reference key
res_id MEDIUMINT, # Residue direct-reference key
tr_id MEDIUMINT, # Transcript direct-reference key
PRIMARY KEY(label,structid,chain,chain_seqid,transcript,trans_seqid),
KEY(al_id),
KEY(as_id),
KEY(res_id),
KEY(tr_id),
KEY(label,structid,chain,chain_seqid),
KEY(label,transcript,trans_seqid)
)