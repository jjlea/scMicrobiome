#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
d <- data.table::fread(args[1], header=T, stringsAsFactors=F)
d2 <- d[-which(d$chr==6 & d$pos > 28477797 & d$pos < 33448354 ),] # hg19/GRCh37
write.table(d2, file=paste(gsub(".tsv","",args[1]),"_removeMHC.tsv", sep=""), sep="\t", quote=F, row.names=F)