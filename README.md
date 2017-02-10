# PDBMap

PDBMap is a command line tool and database interface designed to facilitate the analysis of genetic features in protein structures. This software includes methods for parsing and uploading structural information from the Protein Data Bank, ModBase, and custom protein structural models. It also includes methods for parsing and uploading genetic datasets in VCF, BED, and other file formats. These datasets are then intersected using BEDTools or MySQL to create a direct mapping between each nucleotide and each associated amino acid in all associated protein structures.

PDBMap is a portal between the fields of genetics and structural biology, and as such it relies on several other software packages. A complete list of necessary installations (that cannot be shipped with the PDBMap distribution) are provided below:

* [UCSF Chimera (Headless)](https://www.cgl.ucsf.edu/chimera/cgi-bin/secure/chimera-get.py?file=alpha/chimera-alpha-linux_x86_64_osmesa.bin)
* [MySQL](https://dev.mysql.com/downloads/os-linux.html)
* [Ensembl Core and Variation Database](http://www.ensembl.org/info/docs/webcode/mirror/install/ensembl-data.html)
* [Ensembl Perl API](http://www.ensembl.org/info/docs/api/api_git.html)
* [Ensembl Variant Effect Predictor](https://github.com/Ensembl/ensembl-tools/tree/release/87/scripts)
 * A new beta version of VEP is now available on [github](https://github.com/Ensembl/ensembl-vep)
* [DSSP](ftp://ftp.cmbi.ru.nl/pub/software/dssp/)

Note that all Ensembl resources should use the same genome build and all versions should match. All genomic data loaded into the database must match the Ensembl genome build. All existing resources have been built and maintained using genome build GRCh37/hg19.
