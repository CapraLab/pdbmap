CREATE TABLE IF NOT EXISTS Structure (
label VARCHAR(100),
pdbid VARCHAR(50),
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
str_id MEDIUMINT NOT NULL AUTO_INCREMENT, # Unique, direct-reference key
PRIMARY KEY(label,pdbid),
KEY(str_id),
KEY(label,quality),
KEY(label,method),
KEY(label,resolution),
KEY(label,`release`),
KEY(label,deposition)
)