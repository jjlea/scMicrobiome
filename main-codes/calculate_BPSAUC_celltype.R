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

############ For cell-type-level analysis

num <- dcast(anno, free_annotation2~"num", fun.aggregate = function(x){length(x)})
rownames(num) <- as.character(num$free_annotation2)
#filter the cell type with cell number < 200
celltype_ids <- num$free_annotation2[num$num>200] # about 254 cell-types


### read the BPS scores
INFILE <-  paste0(infile,"/1_df_norm_score.tsv")
score <- data.table::fread(file=INFILE, header = T,  stringsAsFactors = FALSE)
score <- data.frame(score)
rownames(score) <- score$cell_id

### rank the cells
score$cell_id <- NULL
rankscore <- buildRankings(score)

#### calculate the AUC scores
# you may change the parameter of threshold for top cells, e.g., choose the top 1%          cells by changing parameter '0.05' to '0.1' in the following codes. 

# cell type level
Bac_cells_AUC_cell_type <- calc.AUC2(celltype_ids, anno[,c("cell_id","free_annotation2")], rankscore, 0.05)

save(ls(), file=paste0(outfile,"celltype_calculate_BPSAUC.RData"))

