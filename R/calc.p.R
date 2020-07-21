#' Calculate P Values for Interactions Based on Permutations
#' 
#' @param data.obj A \link{\code{Cape}} data object
#' @param pval.correction One of "holm", "fdr", "lfdr" or "none", indicating 
#' whether the p value correction method used should be the Holm step-down procedure, 
#' false discovery rate, local false discovery, or no correction rate respectively.
#' 
#' @return The data object is returned with a new table called var_to_var_p_val.
#' This table is the same as var_to_var_influences, but with p value and adjusted
#' p value columns appended. 
#' 
calc.p <- function(data.obj, pval.correction = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")) {
  
  
  pval.correction = pval.correction[1]
  data.obj$pval_correction <- pval.correction	
  
  influences.org <- data.obj$var_to_var_influences
  influences.perm <- data.obj$var_to_var_influences_perm
  
  if(is.null(influences.org)){
    stop("error.prop() with perm = FALSE must be run before running calc.p()")
  }
  
  if(is.null(influences.perm)){
    warning("error.prop() with perm = TRUE must be run before running calc.p()")
  }
  
  marker.mat <- influences.org[,1:2] #a matrix listing the names of all marker combinations
  colnames(marker.mat) <- c("marker1", "marker2")
  
  
  #### Combinine across permutations#####
  #get the t statistics for all permutations
  m12.null <- as.numeric(influences.perm[,3]) / as.numeric(influences.perm[,4])
  m21.null <- as.numeric(influences.perm[,5]) / as.numeric(influences.perm[,6])
  m.null <- c(m12.null, m21.null)
  
  m12 <- as.numeric(influences.org[,3]) / as.numeric(influences.org[,4])
  m21 <- as.numeric(influences.org[,5]) / as.numeric(influences.org[,6])
  m <- c(m12, m21)
  
  #changed calculation of p value to account for the asymmetric m12/m21 distribution
  #I now calculate the p value based on above and below 0 separately
  low.null <- m.null[which(m.null < 0)]; low.fun <- ecdf(abs(low.null))
  high.null <- m.null[which(m.null > 0)]; high.fun <- ecdf(high.null)
  
  all.emp.p <- matrix(NA, ncol = 1, nrow = length(m))
  high.m.locale <- which(m > 0)
  low.m.locale <- which(m < 0)
  
  low.p <- 1-low.fun(abs(m[low.m.locale]))
  high.p <- 1-high.fun(abs(m[high.m.locale]))
  
  all.emp.p[low.m.locale,1] <- low.p
  all.emp.p[high.m.locale,1] <- high.p
  # hist(all.emp.p)
  
  m12 <- matrix(c(marker.mat[,2],marker.mat[,1],as.numeric(as.matrix(influences.org[,3])),as.numeric(as.matrix(influences.org[,4])),(abs(as.numeric(influences.org[,3])) / as.numeric(influences.org[,4]))), ncol = 5)	
  colnames(m12) <- c("Source","Target","Effect","SE","|Effect|/SE")	
  
  m21 <- matrix(c(marker.mat[,1],marker.mat[,2],as.numeric(as.matrix(influences.org[,5])),as.numeric(as.matrix(influences.org[,6])),(abs(as.numeric(influences.org[,5])) / as.numeric(influences.org[,6]))), ncol = 5)
  colnames(m21) <- c("Source","Target","Effect","SE","|Effect|/SE")
  
  p.adjusted <- p.adjust(all.emp.p, method = pval.correction)
  final.table <- rbind(m12, m21)
  final.table <- cbind(final.table, all.emp.p, p.adjusted)
  colnames(final.table)[c(6,7)] <- c("P_empirical", "p.adjusted")
  final.table <- final.table[order(as.numeric(final.table[,"|Effect|/SE"]), decreasing = TRUE),]
  
  data.obj$var_to_var_p_val <- final.table
  
  return(data.obj)
}
