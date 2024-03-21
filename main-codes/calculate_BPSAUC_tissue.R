#!/usr/bin/env Rscript

library(DelayedArray)
library("ggthemes")


args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]
myproject_file <- args[3]



setwd(outfile)
source(paste0(myproject_file,“/main-codes/AUCell_code.R"))
source(paste0(myproject_file,“/main-codes/myfun.R")



save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

### load the annotation file for scRNA-seq data
load(paste0(myproject_file,"/data/anno.RData")


############ Or tissue level

num <- dcast(anno, organ_tissue~"num", fun.aggregate = function(x){length(x)})
rownames(num) <- as.character(num$organ_tissue)
tissue_ids <- num$organ_tissue

### read the BPS scores
INFILE <-  paste0(infile,"/1_df_norm_score.tsv")
score <- data.table::fread(file=INFILE, header = T,  stringsAsFactors = FALSE)
score <- data.frame(score)
rownames(score) <- score$cell_id

### rank the cells
score$cell_id <- NULL
rankscore <- buildRankings(score)

#### calculate the AUC scores
# you may change the parameter of threshold for top cells, e.g., choose the top 1% cells by changing parameter '0.05' to '0.1' in the following codes. 

# tissue level
Bac_cells_AUC_tissue <- calc.AUC2(tissue_ids, anno[,c("cell_id","organ_tissue")], rankscore, 0.05)

save(ls(), file=paste0(outfile,"tissue_calculate_BPSAUC.RData"))

