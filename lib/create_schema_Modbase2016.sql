# SQL script to create the Modbase2016Model INDEX file in mySQL
# These columns are exactly column headers in
# (/dors/capra_lab/data/modbase/H_sapiens_2016/)
# ....  Homo_sapiens_2016.summary.txt   exactly as downloaded from modbase
#
# Note - The Summary file contains entries for non ENSP... models and sequences.  We ignore those entries
#
# Citation for the columns and models
# "ModBase, a database of annotated comparative protein structure models and associated resources"
# Pieper et al....  Andrej Sali    Nucleic Acids Research, 2013, 111
# doi:10.1093/nar/gkt1144
#
# Additional model quality metrics must be found in the header of the 
# referenced .pdb file.  Chain is located there as well
#
# To load the rows of this file first extract the ENSP* records (filter on column 4)
# Copy over records of interest.  DO NOT Copy over header, as that gets added to SQL unless..
#     grep -P "^.*\t.*\t.*\tENSP.*" \
         /dors/capra_lab/data/modbase/H_sapiens_2016/Homo_sapiens_2016.summary.txt \
         >> /tmp/Modbase2016.tsv

# Run mysqlimport to laod the tabl3e
mysqlimport -h vgi01 -d pdbmap_v14 -u mothcw --ignore --ignore-lines=1 --verbose --fields-terminated-by=
'\t' --local -p  /tmp/Modbase2016.tsv


#   grep -E "
#
#

DROP TABLE IF EXISTS Modbase2016;

CREATE TABLE IF NOT EXISTS Modbase2016 (
run_name VARCHAR(20)        COMMENT 'Not sure - but all were MW-Human_D11',
seq_id   VARCHAR(100)       COMMENT 'Non-Unique id of every modelled sequence (there can be 2+ models per seq_id)',
model_id VARCHAR(100)       COMMENT 'Unique hex identifier for every model',
database_id VARCHAR(100)    COMMENT 'Unique ENSP* strings that are the filenames We ONLY load ENSP* rows',
target_beg INT              COMMENT 'The starting residue # of the ENST.... sequence that is modelled/covered',
target_end INT              COMMENT 'The final residue # of the ENST... sequence that is modelled/covered',
`sequence identity` DOUBLE  COMMENT 'Percent of resdiues modelled that exactly match the template position',
evalue    DOUBLE            COMMENT 'Threhold for target-template alignment.  See refs',
ga341     DOUBLE            COMMENT 'Quality score in Protein Sci., 16, 2412-2426',
mpqs      DOUBLE            COMMENT 'modpipe quality Score in Nucleic Acids Res., 39, 465-474',
zdope     DOUBLE            COMMENT '"Statistical Potential" in Protein Sci., 15, 2507-2524',
pdb_code    CHAR(4)         COMMENT '4 character pdbe/rcsb template lookup',
pdb_chain   VARCHAR(4)      COMMENT 'modeled chain ID from pdb entry',
pdb_beg     INT             COMMENT 'First PDB residue included in model.  May have insertion code letter suffix',
pdb_end     INT             COMMENT 'Last PDB residue included in model.  May have insertion code letter suffix',
`hit history`     CHAR(4)   COMMENT '4 digit number - not sure',
tsvmod_method CHAR(5)       COMMENT '276K MSALL  15K MSRED   102K MTALL',
tsvmod_no35   DOUBLE        COMMENT 'Predicted Native Overlap 3.5A  See Pritein Sci., 17(11), 1881-1893',
tsvmod_rmsd   DOUBLE        COMMENT 'Predicted Calpha RMSD to native structure See Pritein Sci., 17(11), 1881-1893',

PRIMARY KEY(model_id),
KEY(database_id,target_beg)
) ENGINE=InnoDB COMMENT 'Modbase2016 Summary file'
CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci;
