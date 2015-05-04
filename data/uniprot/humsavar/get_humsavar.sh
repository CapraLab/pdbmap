#!/bin/sh

# Download all UniProtKB protein-modifying polymorphisms and disease mutations
# There are 30 descriptive lines before data
wget -q -N --reject -nH -nd --timeout=100000 http://www.uniprot.org/docs/humsavar.txt
