library(reshape2)
#' @name Gambler_R
#' @param test_r Rows should be cells and columns should be genes.
#' @param prob_r  trained scibet model  Rows should be genes and columns should be cell types.The genes must be matched with test_r
#' @param ret_tab   return list or matrix
#' @return  'cellType'  or matrix

Gambler_R <- function(test_r, prob_r,gene_list = NULL){
  remain_genes <- intersect(gene_list,colnames(test_r))
  test_r <- log1p(as.matrix(test_r[,remain_genes])) / log(2)
  prob_r <- prob_r[remain_genes,]
  total <- as.matrix(test_r) %*% as.matrix(prob_r)
  cellType <- c()
  for(i in 1:nrow(total)){
    index <- which.max(total[i,])
    cellType <- c(cellType,colnames(total)[index])
  }
  return(cellType)
}

#' @name celltype_predict
#' @param matrix Rows should be cells and columns should be genes.
#' @return  'cellType'

celltype_predict <- function(matrix){
  l <- data("panC_cell_ref",package="SPECTRUM")
  panC_cell_ref <- eval(parse(text = l))
  prob_r <- t(panC_cell_ref)
  out <- Gambler_R(matrix, prob_r, gene_list = rownames(prob_r))
  return(out)
}

#' @name patient_group_predict
#' @param proportion Rows should be samples and columns should be cell types.
#' @return  'group'

patient_group_predict <- function(proportion){
  l <- data("panC_samp_ref",package="SPECTRUM")
  panC_samp_ref <- eval(parse(text = l))
  prob_r <- t(panC_samp_ref)
  proportion <- proportion * 1000000
  out <- Gambler_R(proportion, prob_r, gene_list = rownames(prob_r))
  return(out)
}

#' @name celltype_predict
#' @param matrix Rows should be cells and columns should be genes.
#' @param sampleID should be a vector of samplID matched with each cell in matrix
#' @return  'group'
merge_predict <- function(matrix, sampleID = NULL){
  cell.type <- celltype_predict(matrix)
  df <- as.data.frame(cbind(cell.type,sampleID))
  colnames(df) <- c("cell.type","sampleID")
  proportion <- table(df$sampleID, df$cell.type)
  proportion <- as.data.frame(proportion / rowSums(proportion))
  colnames(proportion) <- c("sampleID","cell.type","Freq")
  proportion <- reshape2::dcast(proportion, sampleID ~ proportion$cell.type, value.var = "Freq")
  rownames(proportion) <- proportion$sampleID
  proportion <- proportion[,-1]
  group.predict <- patient_group_predict(proportion)
  return(group.predict)
}


