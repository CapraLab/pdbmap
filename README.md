# PDBMap

PDBMap is a command line tool and database interface designed to facilitate the analysis of genetic features in protein structures. This software includes methods for parsing and uploading structural information from the Protein Data Bank, ModBase, and custom protein structural models. It also includes methods for parsing and uploading genetic datasets in VCF, BED, and other file formats. These datasets are then intersected using BEDTools or MySQL to create a direct mapping between each nucleotide and each associated amino acid in all associated protein structures.

PDBMap is a portal between the fields of genetics and structural biology, and as such it relies on several other software packages. A complete list of necessary installations (that cannot be shipped with the PDBMap distribution) are provided below:

Chimera
Ensembl Database
Ensembl API
Ensembl Variant Effect Predictor
DSSP
...
