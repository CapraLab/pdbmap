#!/bin/sh
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l mem=5000mb
#PBS -l walltime=3:00:00:00

# Download all of dbSNP
# Lift from hg38 to hg19
cd /projects/Bush_eQTL/sivleyrm/projects/pdbmap/data/dbsnp;

# Download the BED files into bed/
mkdir hg38; mkdir hg19; cd hg38;
wget -r --no-parent -N --reject -nH -nd --timeout=100000 ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/All.vcf.gz
wget -r --no-parent -N --reject -nH -nd --timeout=100000 ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/All.vcf.gz.tbi
cd ..

# Lift from hg38 to hg19
for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT;
do
  liftOverVCF.pl hg38/bed_chr_$c.bed /projects/Bush_eQTL/sivleyrm/bin/liftover/hg38ToHg19.over.chain.gz hg19/bed_chr_$c_hg19.bed /dev/null;
  gzip hg19/bed_chr_$c_hg19.bed;
done;
