CREATE TABLE IF NOT EXISTS Transcript (
label VARCHAR(100), # Dataset label
transcript VARCHAR(50),
protein VARCHAR(50),
gene VARCHAR(50),
seqid INT,
rescode VARCHAR(1),
chr VARCHAR(50),
start BIGINT,
end BIGINT,
strand INT,
PRIMARY KEY(label,transcript,seqid),
KEY(chr,start,end),
KEY(gene),
KEY(protein)
)
