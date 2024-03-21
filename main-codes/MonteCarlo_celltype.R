#!/usr/bin/env Rscript

library(parallel)

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]

setwd(outfile)

######### cell type level
load(paste0(infile, "/celltype_calculate_BPSAUC.RData"))

### permutation

cl <- makeCluster(20)
clusterExport(cl, c('celltype_ids','num','anno','rankscore',
                  'calc.AUC_ram2','perm_cal2','getset','AUC.cellSet2','auc','calc.AUC2')) 
#### run
clusterEvalQ(cl,{library(AUCell)
  library(reshape2)
  library(data.table)
  library(DelayedMatrixStats)
  }) 
background <- parLapply(cl, celltype_ids, function(x){
    perm_cal2(x,num, 1000, anno, rankscore, 0.05)})

stopCluster(cl) #关闭集群


###  Monte carlo p-values

library(parallel)
cl <- makeCluster(2)
clusterExport(cl, c('Bac_cells_AUC_cell_type','mc_pvalue')) 

clusterEvalQ(cl,{library(ismev)
  library(scanstatistics)
  library(evd)
  }) 


mc.pvalue_0.05 <- parLapply(cl, rownames(Bac_cells_AUC_cell_type), function(x){
    file <- paste("./ramd/",x ,"_ramd_auc.RData", sep="")
  #  file <- gsub("\"", "",file)
    load(file)
    res <- res[c(1:207),]
    r <- lapply(colnames(Bac_cells_AUC_cell_type)[c(1:207)], function(y){
          p <- mc_pvalue(Bac_cells_AUC_cell_type[x,y], res[y,])
         return(p)
    })
    return(r)
})
stopCluster(cl) 


df.mc.p <- data.frame(matrix(unlist(mc.pvalue_0.05), byrow = T, ncol = length(mc.pvalue_0.05[[1]]) ))
colnames(df.mc.p) <- colnames(Bac_cells_AUC_cell_type)[c(1:207)]
rownames(df.mc.p) <- rownames(Bac_cells_AUC_cell_type)


### FDR correction

df.mc.fdr <- as.data.frame(t( data.frame(apply(df.mc.p, 1,function(x){p.adjust(x, method = "BY")}))))   


### merge tables

df <- merge_tb(Bac_cells_AUC_cell_type[,c(1:207)], df.mc.fdr)
df$mcp.adjust <- NULL
df$`-mcp.adjust.log` <- NULL
df$auc <- as.numeric(df$auc)
colnames(df)[3]<- "mc.pvalue.adjusted"

write.table(df, file="AUCell_df_cell.type.tsv", sep="\t, row.names = T, col.names = T)




