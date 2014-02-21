CREATE TABLE IF NOT EXISTS Structure (
pdbid VARCHAR(50),
label VARCHAR(100),
tier INT,
method VARCHAR(100),
quality DOUBLE,
resolution DOUBLE,
`name` TEXT,
author TEXT,
deposition DATE,
`release` DATE,
compound TEXT,
keywords TEXT,
reference TEXT,
structure_reference TEXT,
PRIMARY KEY(label,pdbid),
KEY(pdbid),
KEY(tier,quality),
KEY(method),
KEY(quality),
KEY(resolution),
KEY(`release`)
)