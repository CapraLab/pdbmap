# PDBMap

PDBMap is a command line tool and database interface designed to facilitate the analysis of genetic features in protein structures. This software includes methods for parsing and uploading structural information from the Protein Data Bank, ModBase, and custom protein structural models. It also includes methods for parsing and uploading genetic datasets in VCF, BED, and other file formats. These datasets are then intersected using BEDTools or MySQL to create a direct mapping between each nucleotide and each associated amino acid in all associated protein structures.

PDBMap is a portal between the fields of genetics and structural biology, and as such it relies on several other software packages. A complete list of necessary packages are provided below:

* [UCSF Chimera (Headless)](https://www.cgl.ucsf.edu/chimera/cgi-bin/secure/chimera-get.py?file=alpha/chimera-alpha-linux_x86_64_osmesa.bin)
* [MySQL](https://dev.mysql.com/downloads/os-linux.html)
* [Ensembl Core and Variation Database](http://www.ensembl.org/info/docs/webcode/mirror/install/ensembl-data.html)
* [Ensembl Perl API](http://www.ensembl.org/info/docs/api/api_git.html)
* [Ensembl Variant Effect Predictor](https://github.com/Ensembl/ensembl-tools/tree/release/87/scripts)
 * A new beta version of VEP is now available on [github](https://github.com/Ensembl/ensembl-vep)
* [DSSP](ftp://ftp.cmbi.ru.nl/pub/software/dssp/)

All of these resources must be installed prior to using PDBMap. Note that all Ensembl resources should use the same genome build and all versions should match. All genomic data loaded into the database must match the Ensembl genome build. All existing resources have been built and maintained using genome build GRCh37/hg19.

All other databases and resources required by PDBMap can be downloaded and installed using the following command:
```
./pdbmap.py -c config/DEFAULT.config --refresh
```
This command will download or refresh a local mirror of and/or necessary files from the following databases:
* RCSB Protein Data Bank (solved protein structures)
* ModBase (computationally predicted homology models)
* UniProt (Swiss-Prot, ID mapping, secondary-to-primary AC mapping)
* SIFTS (residue-level functional annotation and pdb-to-reference sequence alignment)
* PFAM (functional domain annotations)

The location of any of these resources can be changed. Users should create a copy of DEFAULT.config with their MySQL login information and the location of all necessary resources. Note that existing copies of these datasets may be used, but the functionallity of `--refresh` may be affected. Parsing of the PDB and ModBase directory structures is also sensitive to change, so consider allowing PDBMap to maintain these datasets.
