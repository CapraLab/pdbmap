CREATE TABLE IF NOT EXISTS Transcript (
transcript VARCHAR(50),
gene VARCHAR(50),
seqid INT,
rescode VARCHAR(1),
chr VARCHAR(20),
start BIGINT,
end BIGINT,
strand INT,
PRIMARY KEY(transcript,seqid),
KEY(chr,start,end),
KEY(gene)
)
