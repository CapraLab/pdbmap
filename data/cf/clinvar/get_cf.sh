#!/bin/sh

# Downloads variants in the CFTR gene 
# and labeled with an association to 
# Cystic fibrosis from ClinVar
#
# Locations are given in buildGRCh38

# Query ClinVar for CFTR variants
# This step is not automated.
echo "http://www.ncbi.nlm.nih.gov/clinvar/advanced"
echo "Specify the [Gene Name] as CFTR"
echo "Download search results a tabular file"
echo "Save in this directory (do not edit filename)"

echo "Converting results to bed-like format..."

# Convert to bed-like format
dos2unix clinvar_result.txt
head -n1 clinvar_result.txt | awk -F"\t" -v OFS="\t" '{print "chr","start","end","name","gene",$4,$5,$6,$10,"association",$7}' > cf.txt
tail -n +2 clinvar_result.txt | awk -F"\t" -v OFS="\t" '{n=split($9,loc,"-"); if (n>1) {print sprintf("chr%d",$8),loc[1],loc[2],sprintf("CFTR:%s",$3),$2,$4,$5,$6,$10,$1,$7} else print sprintf("chr%d",$8),loc[1],loc[1]+1,sprintf("CFTR:%s",$3),$2,$4,$5,$6,$10,$1,$7}' >> cf.txt
# Remove rows without locations
grep -v "^chr0" cf.txt > cfFIX.txt
mv cfFIX.txt cf.txt
# Remove rows where the phenotype is not Cystic fibrosis
awk -F"\t" -v OFS="\t" '{if ($7=="Cystic fibrosis" || $7=="Phenotype") print}' cf.txt > cfFIX.txt
mv cfFIX.txt cf.bed