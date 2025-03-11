
#######################################不去除0
rm(list = ls())

gene_dat1 <- read.csv("ch.csv",row.names = 1)#cold hardiness
#gene_dat1 <- log(gene_dat1+1)

gene_dat2 <- read.csv("cs.csv",row.names = 1)#cold stress

gene_dat1a <- gene_dat1[seq(1,nrow(gene_dat1),2),]
gene_dat1b <- gene_dat1[seq(2,nrow(gene_dat1),2),]

gene_dat2a <- gene_dat2[seq(1,nrow(gene_dat2),2),]
gene_dat2b <- gene_dat2[seq(2,nrow(gene_dat2),2),]









