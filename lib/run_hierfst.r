#!/usr/bin/RScript

# Command Line Arguments:
# ./run_hierfst.r pop-id-map ped map

require(hierfstat)
args <- commandArgs(trailingOnly=T)

# Load populations
df.pop <- read.table(args[1],sep='\t',header=F)

# Convert populations to numeric IDs
pops <- unique(unlist(df.pop[,2]))
numpop <- apply(df.pop,1,function(x){match(x[2],pops)})
df.pop[,2] <- numpop

# Load genotypes
df.gen <- read.table(args[2],sep='\t',header=F,na.strings='-9')

# Load the map
df.map <- read.table(args[3],header=F)
names(df.gen)[1] <- "pop"
names(df.gen)[-1] <- df.map[-1]

# Merge data
df.popgen <- merge(df.pop,df.gen,by="V1")

# Reduce to pop column and genotype columns
df.popgen <- df.popgen[,-c(1,3,4,5,6,7)]

# Run HierFst
fstat <- basic.stats(df.popgen)

