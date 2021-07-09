# PDBMap

PDBMap is a python library which manages the intersection of 
- Genetic Variant Database in .vcf (or .bed) format (Gnomad/Clinvar/COSMIC/etc)
- The ENST* Transcript ID set fand the PERL API of the ENSEMBL human genome database
- Uniprot IDs and attributes, from the curated Uniprot IDMapping file
- SIFTs alignment XML files for canonical transcript IDs and REST API results for non-canonical uniprot transcript ID
- The ENSEMBL Variant Effect Predictor

The datasets which are integrated are a combination of SQL tables created by pdbmap.py, and a variety of raw data files downloaded to the filesystem.

Typically PSB Pipeline users will download the various resources discussed at https://github.com/CapraLab/pipeline_download_scripts and point their .config
files to the SQL Server at Vanderbilt.

Therefore, the below documentation primarily targets PDBMap software maintainers, and those who wish to populate their own local SQL Server.

The PDBMap library (lib/PDBMap*.py) is a critical component of the Personal Structural Biology (PSB) Pipeline.

Loading of variant datasets into a SQL database is accomplished by the pdbmap.py command line.

Legacy features of the library included visualation and analysis of genomic features in context of protein structures.

However, as the library has evolved, Protein structures (pdb/modbase/swiss) are no longer stored in the SQL database.
Rather the lib/PDBMap{Protein,Modbase,Swiss}.py modules load these structures directly from the filesystem.


## PDBMap External Dependencies

PDBMap is a portal between the fields of genetics and structural biology, and as such it relies on several (free) software packages and databases. 
Unfortunately, these cannot be distributed with PDBMap and must be installed separately.

Typically, these dependencies are satisfied with the sequence of steps required to instantiate the Personal Structural Biology pipeline.

A complete list of the packages are provided below:

* [Python 3.x] (We recommend [Anaconda](https://www.anaconda.com/products/individual))
* Biopython (pip install biopython after Python 3.x)
* A SQL database server, MySQL or MariaDB: [MySQL](https://dev.mysql.com/downloads/os-linux.html)
* [Ensembl Core Database](http://www.ensembl.org/info/docs/webcode/mirror/install/ensembl-data.html) (use of the Ensembl public MySQL server is supported, but may result in slow runtimes)
* [Ensembl Perl API](http://www.ensembl.org/info/docs/api/api_git.html)
* [Ensembl Variant Effect Predictor (and cache)](https://github.com/Ensembl/ensembl-tools/tree/release/87/scripts) (a new beta version of VEP is now available on [github](https://github.com/Ensembl/ensembl-vep))
* [UCSF Chimera (Headless)](https://www.cgl.ucsf.edu/chimera/cgi-bin/secure/chimera-get.py?file=alpha/chimera-alpha-linux_x86_64_osmesa.bin) (for visualization)
* [DSSP](http://swift.cmbi.ru.nl/gv/dssp/) (for secondary structure and solvent accessibility)

All of these resources must be installed prior to using PDBMap. Note that all Ensembl resources (PERL API, Variant Effect Predictor (VEP), and VEP cache files)
must use the same genome build and all versions should match. All genomic data loaded into the database must match the Ensembl genome build. All existing resources have been built and maintained using genome build **GRCh38**


It is helpful to install PDBMap on a SLURM cluster. Many PDBMap tasks, such as loading varaints from 24 separate .chrNN.vcf files, lend themselves well to parallelization. SLURM scripts for many common tasks are provided for convenience. Before launching many jobs to the cluster, check that your MySQL server is configured to manage the corresponding number of connections.

## Instantiating the PDBMap MAriaDB/MySQL Database

To instantiate the PDBMap database, a database of tables must be created in SQL.

For routine usage, it is desirable to have a user id which lacks change priviledges.  But you will need extensive credentials 
to CREATE a database, and then CREATE tables and insert data rows.

If you are new to SQL, it could be helpful to install the ENSEMBL sql tables first, since these processes, while complex, 
are more structured than the table installs for the PDBMap library.

We strongly recommend using a full-feature interactive sql client, such as JetBrains DataGrip, or the MySqlWorkbench.  Because
the mysql command line tool is universal, we give the example using mysql.  (It is also true that some operations are faster when run
through mysql on the same compute node/PC that the SQL Server is instantiated on)

Assuming you have extensive CREATE/DROP/INSERT credentials on the sql server, launch the mysql linux client, example:

   $ mysql -h vgi01.accre.vanderbilt.edu -p

You will enter your sql password and your screen will then like something like:

    Enter password: 
    Welcome to the MySQL monitor.  Commands end with ; or \g.
    Your MySQL connection id is 277870
    Server version: 5.5.68-MariaDB MariaDB Server

    Copyright (c) 2000, 2016, Oracle and/or its affiliates. All rights reserved.

    Oracle is a registered trademark of Oracle Corporation and/or its
    affiliates. Other names may be trademarks of their respective
    owners.

    Type 'help;' or '\h' for help. Type '\c' to clear the current input statement.

    mysql> 
    
Continue with 
   
    mysql> CREATE DATABASE pdbmap_v15;
    
where pdbmap_v15 is the database name you are placing in your .config file(s).

## Creating PDBMap Tables

Most table creation is performed by the pdbmap/lib/create_schemea_*.sql scripts.  These critical scripts not only
create initialized tables with well-commented SQL table schemas.  The text of the scripts often includes 
comments about how to populate the scripts using other available utilities or scripts.

A guiding principle of the PDBMap SQL database is that, whenever possible, the Table Rows should directly copy, untranslated, 
source datafile rows as downloaded.  The mysqlimport utility is an excellent tool for this.  When mysqlimport is insufficient, 
short, transparently written python scripts, bridge the gap.

### Uniprot IDMapping File
PDBMap and the pipeline make extensive use of the 3 column IDMapping file, The Download .bash script for data/uniprot

https://github.com/CapraLab/pipeline_download_scripts/uniprot/DOWNLOAD_uniprot.bash

includes a python module which creates HUMAN_9606_idmapping_sprot.dat.gz.  Importantly, this file is a reduced list of Swissprot curated
uniprot IDs, and excludes general Trembl entries.

Run lib/create_schema_Idmapping.sql, either as standard input to mysql, or by pasting into a GUI tool.  The note the instructions 
in the file to load table rows

$ cp data/HUMAN_9606_idmapping_sprot.dat.gz /tmp/Idmapping.gz
$ gunzip /tmp/Idmapping.gz
$ mysqlimport --ignore --verbose --fields-terminated-by='\t' --local -p  --columns=unp,ID_type,ID pdbmap_v15 /tmp/Idmapping

### Uniparc Invariant sequence Idnetifiers
Every uniprot protein identifier is cross-referenced to a uniparc ID which uniquely specifies an amino acid sequence.
From Uniprot release to release, the cross-referecned uniparc ID might well change.  However, the amino acid sequence
mapped to  a Uniprot ID (ex. UPI1234567890123) will NEVER change.

Thus, it is never necessary to "DROP" the Uniparc table before loading.  The scripts/uniparc_parser.py program will always 
add rows with the "ignore flag".

For a true first time table create, the script scripts/create_schema_Uniprot.sql has the needed instructions.

An example commandline to update (or initially populate) the table following a fresh uniprot download is:

$ ./uniparc_parser.py -c /mylocation/config/global.config ...../somewhere/data/uniprot/current/uniparc_active.fasta.gz 100000

See the available --help option for more information on how to run uniparc_parser.py

### Sifts .xml files and 





## Loading Genomic Information into PDBMap
Any genomic dataset can be loaded into PDBMap. By default, scripts are provided to download local copies of variant data from
* The Single Nucleotide Polymorphism Database (dbSNP)
* The 1000 Genomes Project
* The Exome Sequencing Project (ESP)
* The Exome Aggregation Consortium (ExAC)
* The Catalogue of Somatic Mutations in Cancer (COSMIC)
* The GWAS Catalogue
* ClinVar
Each of these datasets can be downloaded by running the `get_<dataset>.sh` script within the corresponding directories.

To load a genomic dataset into PDBMap, use
```
./pdbmap.py -c config/<USER>.config load_data <data_file> <data_name>
OR
./pdbmap.py -c config/<USER>.config --dlabel=<data_name> load_data <data_file> [<data_file> ...]
```
Genetic datasets are often distributed by-chromosome and are thus easily parallelizable. SLURM scripts for some of the default datasets are provided and may be used a templates for designing SLURM scripts for other datasets. To load data in parallel, use
```
sbatch --array=1-24 slurm/load_exac.slurm 
```

## Intersecting Structural and Genomic Information
Once the structural and genomic datasets have each been loaded into PDBMap, they must be intersected to construct the direct mapping from nucleotide to amino acid. This can be a lengthy process, but it must only be performed once for each dataset. Once intersected, queries are very efficient, enabling large-scale, high-throughput analysis of genetic information within its protein structural context. To intersect two datasets, use
```
./pdbmap.py -c config/<USER>.config --slabel=pdb --dlabel=exac intersect
```
This command download the structural and genomic data to flat files indexed by chromosomal position, perform an intersection using `intersectBed`, and upload the results back to the database. If you are working with smaller datasets, you may consider adding the `quick` flag after `intersect`. This will perform the intersection using a MySQL join instead of `intersectBed`, which may decrease runtime. This is highly discouraged for larger datasets.



The join order for the genomic tables is:
```
GenomicData -> GenomicConsequence -> GenomicIntersection
            -> Supplementary Tables
```
Most genetic datasets are uploaded in their original form to supplemental tables prior to being processed and uploaded into PDBMap. Joined with GenomicData on their chromosomal position, these tables allow users to incorporate additional information not supported by the default PDBMap schemas.

A detailed layout of the PDBMap database schema is provided **here**.

## Previous architecture: Historical footnote

Prior to an overhaul by Chris Moth which was started in 2019, the library loaded structural information in SQL tables (chain/residue/etc) 
and pre-aligned to ENSEMBL transcripts.

This ambitious architecture allowed for some advanced global analyses to be more easily performed.  However, it made updating the various structural sources, 
as well as the ENSEMBL Human Genome, more involved.  

These datasets are then intersected using BEDTools or MySQL to create a direct mapping between each nucleotide and each associated amino acid in all associated protein structures. A schematic overview of the [PDBMap Pipeline](./docs/PDBMapPipeline.png) is provided.

Once the PDBMap database has been loaded, the PDBMap library can be imported from other projects and used to interface with and analyze genetic data within solved structures and computational models of human proteins.

### pdbmap.py used to be able to download all external databases.
This is an ambitious feature, as external database file formats and download instructions tend to shift.
For historical reference the old instructions to accomplish this were:


Once complete, PDBMap will prompt users to install all other requried databases and resources. If you would like to install these resources later, or bring the local copies up-to-date, use:
```
./pdbmap.py -c config/<USER>.config --refresh
```
This command will download or refresh a local mirror of and/or necessary files from the following databases:
* RCSB Protein Data Bank (solved protein structures)
* ModBase (computationally predicted homology models)
* UniProt (Swiss-Prot, ID mapping, secondary-to-primary AC mapping)
* SIFTS (residue-level functional annotation and pdb-to-reference sequence alignment)
* PFAM (functional domain annotations)

The location of any of these resources may be changed. Users should update DEFAULT.config with location of all necessary resources. Note that existing copies of these datasets may be used, but the functionallity of `--refresh` may be affected. Parsing of the PDB and ModBase directory structures is also sensitive to change, so consider downloading the datasets with PDBMap and then moving the directories into a shared location; update the configuration file with the new location.


### Loading Structural Information into PDBMap
To load **only** protein structures from the Protein Data Bank into PDBMap, use
```
./pdbmap.py -c config/<USER>.config load_pdb all
```
To load **only** protein structural models from ModBase into PDBMap, use
```
./pdbmap.py -c config/<USER>.config load_model all
```
To load **all** PDB structures and ModBase models for all Swiss-Prot human proteins (recommended), use
```
./pdbmap.py -c config/<USER>.config load_unp all
```
In the database, all PDB structures receive the label `pdb` and all ModBase models receive the label `modbase` unless otherwise specified.

To load the entire PDBMap structural database in parallel using `N` SLURM jobs, update `slurm/load_pdbmap.slurm` with your SLURM account information, then use,
```
sbatch --array=0-N slurm/load_pdbmap.slurm N
```

### Visualizing Genomic Information in Structure
The visualization capabilities of PDBMap are built around Chimera. Any property of a genomic dataset can be visualized by specifying the dataset and property name along with the specified gene, protein, or structure name. For example, to visualize the location of all ExAC missense variants in the first biological assembly of 2SHP, use
```
./pdbmap.py -c config/<USER>.config visualize 2SHP exac . 1
```
A new directory will be added to the `results/` directory containing several files, including a Chimera scene and an automatically generated png image of the variant-to-structure mapping. If you would instead like to color each missene variant by its minor allele frequency, 
```
./pdbmap.py -c config/<USER>.config visualize 2SHP exac maf 1
```
You can also compare the distribution to the synonymous distribution of minor allele frequencies,
```
./pdbmap.py -c config/<USER>.config visualize 2SHP exac maf.synonymous 1
```

### Navigating the PDBMap MySQL Database
The MySQL database is composed of several tables, generally organized into structural tables, genomic tables, the intersection table, and supplementary tables. 
The join order for the structure tables is:
```
Structure
          -> Chain -> Residue -> Alignment -> Transcript
Model              -> AlignmentScore
                 
```
Each table is joined on its common columns. For example, Residue and Chain are joined on matching `slabel`, `structid`, and `chain` columns.
