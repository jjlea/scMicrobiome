#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])
z <- read.table("zscore_merged.txt",header=T, stringsAsFactors=F)
ID <- read.table(paste0(args[2],"/magma/locations/NCBI37.3.gene.loc"), header=F, stringsAsFactors=F)
k <- ID[match(z$ID, ID$V1),"V6"]
z$ID <- k
z2 <- z[-nrow(z),]
write.table(z2, file="zscore_merged2.txt", sep="\t", row.names=F, quote=F)

