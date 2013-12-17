#!/usr/bin/Rscript

require(fpc)

args <- commandArgs(TRUE)

# cat("enter pdbid: ")
stdin <- file("stdin")
pdbid <- args[1]
# pdbid <- readLines(stdin,1)
fin  <- sprintf('/scratch/sivleyrm/pdbmap/scratch/%s_variant_positions.txt',pdbid)
df <- read.table(fin,header=TRUE,sep='\t')

# Define the hclust output
hclust.prot.tiff <- sprintf('/scratch/sivleyrm/pdbmap/scratch/%s_prot_hclust.tiff',pdbid)
hclust.prot.fout <- sprintf('/scratch/sivleyrm/pdbmap/scratch/%s_prot_hclust.txt',pdbid)
hclust.gene.tiff <- sprintf('/scratch/sivleyrm/pdbmap/scratch/%s_gene_hclust.tiff',pdbid)
hclust.gene.fout <- sprintf('/scratch/sivleyrm/pdbmap/scratch/%s_gene_hclust.txt',pdbid)

# Run hclust in protein space
hclust.res <- hclust(dist(df[,c("x","y","z")],method="euclidean"))
tiff(hclust.prot.tiff)
plot(hclust.res,labels=df$var_name)
dev.off()
cat("# clusters: ")
clus.num <- as.integer(readLines(stdin,1))
clusterids <- cutree(hclust.res,k=clus.num)
clusterids
prot.hclust <- data.frame(pdbid=pdbid,chain=df$chain,seqres=df$chain_seq,clusterid=clusterids,name=df$var_name)
write.table(prot.hclust,hclust.prot.fout,sep='\t',row.names=FALSE,quote=FALSE)

# Run hclust in gene space
hclust.res <- hclust(dist(df[,c("chr","start")],method="euclidean"))
tiff(hclust.gene.tiff)
plot(hclust.res,labels=df$var_name)
dev.off()
cat("# clusters: ")
clus.num <- as.integer(readLines(stdin,1))
clusterids <- cutree(hclust.res,k=clus.num)
clusterids
gene.hclust <- data.frame(pdbid=pdbid,chain=df$chain,seqres=df$chain_seq,clusterid=clusterids,name=df$var_name)
write.table(gene.hclust,hclust.gene.fout,sep='\t',row.names=FALSE,quote=FALSE)

# cat("enter protein eps: ")
# prot.eps <- readLines(stdin,1)
# cat("enter gene eps: ")
# gene.eps <- readLines(stdin,1)
# close(stdin)

# # Define the dbscan output
# prot.fout <- sprintf('%s_dbscan_eps%s.prot.clus',pdbid,prot.eps)
# gene.fout <- sprintf('%s_dbscan_eps%s.gene.clus',pdbid,gene.eps)
# dbscan.prot.pdf <- sprintf('%s_dbscan.prot.pdf',pdbid)
# dbscan.gene.pdf <- sprintf('%s_dbscan.gene.pdf',pdbid)

# # Run dbscan over protein space
# pdf(dbscan.prot.pdf)
# dist.mat <- dist(df[,c("x","y","z")],method="euclidean")
# dist.mat <- as.matrix(dist.mat,nrow=nrow(df),ncol=nrow(df))
# clus.res <- dbscan(dist.mat,eps=prot.eps,MinPts=2,method="dist",showplot=1)
# dev.off()
# prot.clus <- data.frame(pdbid=pdbid,chain=df$chain,seqres=df$chain_seq,clusterid=clus.res$cluster,name=df$var_name)
# write.table(prot.clus,prot.fout,sep='\t',row.names=FALSE,quote=FALSE)

# # Run dbscan over genomic position
# pdf(dbscan.gene.pdf)
# dist.mat <- dist(df[,c("chr","start")],method="euclidean")
# dist.mat <- as.matrix(dist.mat,nrow=nrow(df),ncol=nrow(df))
# clus.res <- dbscan(dist.mat,method="dist",eps=gene.eps,MinPts=2,showplot=1)
# dev.off()
# gene.clus <- data.frame(pdbid=pdbid,chain=df$chain,seqres=df$chain_seq,clusterid=clus.res$cluster,name=df$var_name)
# write.table(gene.clus,gene.fout,sep='\t',row.names=FALSE,quote=FALSE)