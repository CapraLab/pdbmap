#!/bin/sh

# Pull the UniProt primary AC -> dbref idmapping
wget -q -N --reject -nH -nd --timeout=100000 ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz
# Decompress the idmapping
gunzip -c HUMAN_9606_idmapping_selected.tab.gz > HUMAN_9606_idmapping_selected.tab
# Grab only the UniProt, RefSeq, PDB, Ensembl Transcript, and Ensembl Protein Columns
cut -f 1,4,6,20,21 HUMAN_9606_idmapping_selected.tab > HUMAN_9606_idmapping_UNP-RefSeq-PDB-Ensembl.tab

# Pull the UniProt secondary AC -> primary AC mapping

wget -q -N --reject -nH -nd --timeout=100000 ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/docs/sec_ac.txt
# Remove the header and convert whitespace to tab
sed '1,30d' sec_ac.txt | sed 's/ \+ /\t/g' > uniprot_sec2prim_ac.txt

# Cleanup intermediate files
rm -f HUMAN_9606_idmapping_selected.tab
