CREATE TABLE IF NOT EXISTS Structure (
pdbid VARCHAR(50),
tier INT,
method VARCHAR(100)
quality DOUBLE,
resolution DOUBLE,
release DATE,
keywords TEXT,
reference TEXT,
PRIMARY KEY(pdbid),
KEY(tier,quality),
KEY(method),
KEY(quality)
KEY(resolution)
)