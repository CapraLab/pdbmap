# SQL script to store the SwissModel INDEX file in mySQL
# The corresponding .pdb.gz files are in
# "rootdir/" + unp[0:2] + "/" + unp[2:4] + "/" + unp[4] + "/swissmodel/" +
# "%d_%d"%(start,end) + template_pdbid + coordinate_id + ".pdb.gz"
#
# Additional model quality metrics must be found in the header of the 
# referenced .pdb file.  Chain is loated there as well
#

DROP TABLE IF EXISTS Swiss;

CREATE TABLE IF NOT EXISTS Swiss (
label VARCHAR(100),     # Always swiss
modelid VARCHAR(100),   # Catenation of "%s_%d_%d_%s"%(unp,start,end,template)
unp VARCHAR(100),       # Uniprot
start INT,              # Always an integer in SwissModel INDEX file
end INT,                # Always an integer in SwissModel INDEX file
qmean DOUBLE,           # Swiss-model quality information
qmean_norm DOUBLE,      # Swiss-model quality info
template VARCHAR(50),  # Typically format pdb_id.[1-9].[A-Z]  * Do I need to parse 
coordinate_id VARCHAR(32),   # 24 character guid-like string needed to find pdb filename     
# Next few fields come from REMARK 3 of the SWISS model .pdb file.  They are not in the meta file
pdbid VARCHAR(10),      # REMARK 3: The rcsb-deposited experimental structure used to template the model
chain VARCHAR(10),      # REMARK 3: The original chain letter(s) from the rcsb-deposited template
ostat VARCHAR(50),      # REMARK 3: Text that suggests we have a dimer or higher order to deal with
mthd VARCHAR(50),       # REMARK 3: Method, as in X-ray or NMR and resolution
identity double,        # REMARK 3: Sequence identity transcript to model
url VARCHAR(200),       # Location to save in case direct downlad from Swiss-Model becoes desirable
str_id BIGINT NOT NULL AUTO_INCREMENT, # Unique, direct-reference key  
PRIMARY KEY(label,modelid),
KEY(str_id),
KEY(label,unp)
);

