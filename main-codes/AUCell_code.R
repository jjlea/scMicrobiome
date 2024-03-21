library(reshape2)


# x

###############################################################
## 1. 产生matrix的rank

buildRankings <- function(exprMat){
  rowNames <- rownames(exprMat)
  colNames <- colnames(exprMat)
  exprMat <- -exprMat # ro rank decreasingly
  exprMat <- DelayedMatrixStats::colRanks(DelayedArray(exprMat), ties.method="random", preserveShape=TRUE) 
  exprMat <- t(exprMat)
  rownames(exprMat) <- rowNames
  colnames(exprMat) <- colNames
  return(exprMat)
}

############################################################
# 2. get cell ids for each cell-type name
# x: cell-type name; 
# mat: matrix of cell-type * cellname, column1=cell_id, column2 = cell_type

getset <- function(x, mat){
  cells <- mat[which(mat[,2]==x),1]
  return(cells)
}


##################################################################
#3. calc AUC


AUC.cellSet <- function(cellSet, rankings, aucMaxRank)
{ 
  
  cellSet <- unique(cellSet)
  ncells <- length(cellSet)
 
  cellSet <- cellSet[which(cellSet %in% rownames(rankings))]
  gSetRanks <- rankings[which(rownames(rankings) %in% cellSet),,drop=FALSE]
  rm(rankings)
  
  aucThreshold <- round(aucMaxRank)
  maxAUC <- aucThreshold * nrow(gSetRanks)  # database.cell_count 
 
  
  # Apply by columns (i.e. to each ranking)
  auc <- apply(gSetRanks, 2, auc, aucThreshold, maxAUC)
  return(c(auc,  ncells=ncells))
}


AUC.cellSet2 <- function(cellSet, rankings, aucMaxRank){
  cellSet <- unique(cellSet)
  ncells <- length(cellSet)
  cellSet <- cellSet[which(cellSet %in% rownames(rankings))]
  missing <- ncells-length(cellSet)
  
  gSetRanks <- rankings[which(rownames(rankings) %in% cellSet),,drop=FALSE]
  rm(rankings)
  
  aucThreshold <- round(aucMaxRank)
  ########### NEW version:  #######################
  x_th <- 1:nrow(gSetRanks)
  x_th <- sort(x_th[x_th<aucThreshold])
  y_th <- seq_along(x_th)
  maxAUC <- sum(diff(c(x_th, aucThreshold)) * y_th) 
  ############################################
  
  # Apply by columns (i.e. to each ranking)
  auc <- apply(gSetRanks, 2, auc, aucThreshold, maxAUC)
  
  return(c(auc, missing=missing, ncells=ncells))
}




auc <- function(oneRanking, aucThreshold, maxAUC)
{
  x <- unlist(oneRanking)
  x <- sort(x[x<aucThreshold])
  y <- seq_along(x)
  sum(diff(c(x, aucThreshold)) * y)/maxAUC
}



############################################################
### 4. calculate AUC scores
# ids: uniq list of cell_types
# anno: annotation of cell_types, two columns, column1=cell_id, column2=cell_type
# percent: the percentage of cells to enrich, e.g. 0.05=50%, 0.01=1%
calc.AUC <- function(ids, anno,  cells_rankings, percent){
  res <- data.frame(lapply(ids, function(x){
    sub <- getset(x, anno)
    res <- AUC.cellSet(cellSet=sub, rankings=cells_rankings, aucMaxRank=ceiling(percent*nrow(cells_rankings))) # 0.05 表示取top5%的cell
    return(res)
  }))
  colnames(res) <- ids
  res <- data.frame(t(res))
  return(res)
}



calc.AUC2 <- function(ids, anno,  cells_rankings, percent){
  res <- data.frame(lapply(ids, function(x){
    sub <- getset(x, anno)
    res <- AUC.cellSet2(cellSet=sub, rankings=cells_rankings, aucMaxRank=ceiling(percent*nrow(cells_rankings))) # 0.05 表示取top5%的cell
    return(res)
  }))
  colnames(res) <- ids
  res <- data.frame(t(res))
  return(res)
}


calc.AUC3 <- function(ids, anno,  cells_rankings, N){
  res <- data.frame(lapply(ids, function(x){
    sub <- getset(x, anno)
    res <- AUC.cellSet2(cellSet=sub, rankings=cells_rankings, aucMaxRank=N) # N 表示取N个cell
    return(res)
  }))
  colnames(res) <- ids
  res <- data.frame(t(res))
  return(res)
}

# id: cell_type id
# num: dataframe for number of cells of each cell_type id
# n : number of permutations
# anno: annotation of cell_type and cell_id
# rankscore: rankscore
# percent: 0.01 or 0.05

perm_cal <- function(id, num, n, anno, rankscore, percent){
  res <- list()
  for (k in 1:n){
    rannum <- round(runif (num[id,"num"], min = 1, max = 483152))
    ranid <- anno$cell_id[rannum]
    res[[k]] <- round(calc.AUC_ram(ranid, rankscore, percent),5)
    
  }
  res <- data.frame(res)
  colnames(res) <- paste("ramd",c(1:n), sep = "")
  save(res, file=paste(id,"_ramd_auc.RData", sep=""))
  return(id)
}

perm_cal2 <- function(id, num, n, anno, rankscore, percent){
  res <- list()
  for (k in 1:n){
    rannum <- round(runif (num[id,"num"], min = 1, max = 483152))
    ranid <- anno$cell_id[rannum]
    res[[k]] <- round(calc.AUC_ram2(ranid, rankscore, percent),5)
    
  }
  res <- data.frame(res)
  colnames(res) <- paste("ramd",c(1:n), sep = "")
  save(res, file=paste(id,"_ramd_auc.RData", sep=""))
  return(id)
}


calc.AUC_ram <- function(ids,  cells_rankings, percent){
  res <- AUC.cellSet(cellSet=ids, rankings=cells_rankings, aucMaxRank=ceiling(percent*nrow(cells_rankings))) # 0.05 表示取top5%的cell
  return(res)
}


calc.AUC_ram2 <- function(ids,  cells_rankings, percent){
  res <- AUC.cellSet2(cellSet=ids, rankings=cells_rankings, aucMaxRank=ceiling(percent*nrow(cells_rankings))) # 0.05 表示取top5%的cell
  return(res)
}

mc_pvalue <- function(observed, replicates) {
  replicates <- as.vector(as.matrix(replicates))
  if (length(replicates) == 0) {
    return(NULL)
  } else {
    f <- Vectorize(
      function(y) {
        (1 + sum(replicates > y)) / (1 + length(replicates))
      }
    )
    
    return(f(observed))
  }
}



################ merge table 
# x=Bac_cells_AUC_0.05
# y=df.mc.p
merge_tb <- function(x,y){

	x$ncells <- NULL
	x$id <- rownames(x)
	auc <- reshape2::melt(x ,id.vars = "id")
	
	### p value table 
	y <- data.frame(y)
	y$ncells <- NULL
	y$id <- rownames(y)
	qval3 <- reshape2::melt(y, id.vars = "id")
	qval3$mcp.adjust <- p.adjust(qval3$value, method = "fdr")
	qval3$mcp.adjust.log <- -log10(qval3$mcp.adjust)
	
	### merged table
	colnames(qval3) <- c("id","variable","mc.pvalue","mcp.adjust","-mcp.adjust.log")
	df <- merge(qval3, auc, all = T)
	colnames(df)[6] <- "auc"
	df$ID <- unlist(apply(df,1, function(x){paste(x[1],x[2], sep=" <> ")}))
	df$ID <- lapply(df$ID, function(x){gsub(".genes.score.gz","",x)})
	return(df)

}

