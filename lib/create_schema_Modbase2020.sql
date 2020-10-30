# SQL script to replicate the modbase Homo_sapiens_2020.summary.txt file in mySQL
# The columns are exactly column headers in
# (/dors/capra_lab/data/modbase/current/)
# ....  Homo_sapiens_2020.summary.txt   exactly as downloaded from modbase
#
# Citation for the columns and models
# "ModBase, a database of annotated comparative protein structure models and associated resources"
# Pieper et al....  Andrej Sali    Nucleic Acids Research, 2013, 111
# doi:10.1093/nar/gkt1144
#
# Additional model quality metrics must be found in the header of the 
# referenced .pdb file.  Chain is located there as well
# 
#############################################################################################
# After you have run this script, you need to populate the table with a few shell commands
# mysqlimport only deals with filename=tablename.  So create a symlink
/*
cd /dors/capra_lab/data/modbase/current
ln -s Homo_sapiens_2020.summary.txt Modbase2020.txt
*/
# Then Run mysqlimport to laod the Modbase2020 table
/*
mysqlimport -h vgi01 -d pdbmap_v14 -u mothcw \
--ignore --ignore-lines=1 --verbose --fields-terminated-by='\t' \
--local -p  /dors/capra_lab/data/modbase/current/Modbase2020.txt
*/
#
#
DROP TABLE IF EXISTS Modbase2020;

CREATE TABLE IF NOT EXISTS Modbase2020 (
run_name VARCHAR(20)        COMMENT 'Not sure - but all were MW-Human_D11',
seq_id   VARCHAR(100)       COMMENT 'Non-Unique precise sequence id (first 4 AAs of seq + md5sum of  all AAs + last 4 AA letters)',
model_id VARCHAR(100)       COMMENT 'Unique hex identifier for every model',
database_id VARCHAR(100)    COMMENT 'Unique ENSP*/NP/etc appended identifiers to find models in filesystem',
target_beg INT              COMMENT 'The starting residue # of the ENST.... sequence that is modelled/covered',
target_end INT              COMMENT 'The final residue # of the ENST... sequence that is modelled/covered',
`sequence identity` DOUBLE  COMMENT 'Percent of resdiues modelled that exactly match the template position',
evalue    DOUBLE            COMMENT 'Threhold for target-template alignment.  See refs',
ga341     DOUBLE            COMMENT 'Quality score in Protein Sci., 20, 2412-2426',
mpqs      DOUBLE            COMMENT 'modpipe quality Score in Nucleic Acids Res., 39, 465-474',
zdope     DOUBLE            COMMENT '"Statistical Potential" in Protein Sci., 15, 2507-2524',
pdb_code  CHAR(4)           COMMENT '4 character pdbe/rcsb template lookup',
pdb_chain VARCHAR(4)        COMMENT 'modeled chain ID from pdb entry',
pdb_beg   VARCHAR(20)       COMMENT 'First PDB residue included in model.  May have insertion code letter suffix',
pdb_end   VARCHAR(20)       COMMENT 'Last PDB residue included in model.  May have insertion code letter suffix',
`hit history` CHAR(4)       COMMENT '4 digit number - not sure',
tsvmod_method CHAR(5)       COMMENT '354K MSALL  22K MSRED   97K MTALL',
tsvmod_no35   DOUBLE        COMMENT 'Predicted Native Overlap 3.5A  See Pritein Sci., 17(11), 1881-1893',
tsvmod_rmsd   DOUBLE        COMMENT 'Predicted Calpha RMSD to native structure See Pritein Sci., 17(11), 1881-1893',

PRIMARY KEY(model_id),
KEY(database_id,target_beg),
KEY(seq_id)
) ENGINE=InnoDB COMMENT 'Modbase2020 Summary file'
CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci;
