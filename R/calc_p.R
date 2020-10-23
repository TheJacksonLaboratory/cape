#' Calculate P Values for Interactions Based on Permutations
#' 
#' @param data_obj A \code{\link{Cape}} data object
#' @param pval_correction One of "holm", "fdr", "lfdr" or "none", indicating 
#' whether the p value correction method used should be the Holm step-down procedure, 
#' false discovery rate, local false discovery, or no correction rate respectively.
#' 
#' @return The data object is returned with a new table called var_to_var_p_val.
#' This table is the same as var_to_var_influences, but with p value and adjusted
#' p value columns appended.
#' 
#' @importFrom stats p.adjust
#' 
#' @export
calc_p <- function(data_obj, pval_correction = 
                     c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")) {
  
  pval_correction = pval_correction[1]
  data_obj$pval_correction <- pval_correction	
  
  influences_org <- data_obj$var_to_var_influences
  influences_perm <- data_obj$var_to_var_influences_perm
  
  if(is.null(influences_org)){
    stop("error_prop() with perm = FALSE must be run before running calc_p()")
  }
  
  if(is.null(influences_perm)){
    warning("error_prop() with perm = TRUE must be run before running calc_p()")
  }
  
  marker_mat <- influences_org[,1:2] #a matrix listing the names of all marker combinations
  colnames(marker_mat) <- c("marker1", "marker2")
  
  
  #### Combinine across permutations#####
  #get the t statistics for all permutations
  m12_null <- as.numeric(influences_perm[,3]) / as.numeric(influences_perm[,4])
  m21_null <- as.numeric(influences_perm[,5]) / as.numeric(influences_perm[,6])
  m_null <- c(m12_null, m21_null)
  
  m12 <- as.numeric(influences_org[,3]) / as.numeric(influences_org[,4])
  m21 <- as.numeric(influences_org[,5]) / as.numeric(influences_org[,6])
  m <- c(m12, m21)
  
  #changed calculation of p value to account for the asymmetric m12/m21 distribution
  #I now calculate the p value based on above and below 0 separately
  low_null <- m_null[which(m_null < 0)]; low_fun <- ecdf(abs(low_null))
  high_null <- m_null[which(m_null > 0)]; high_fun <- ecdf(high_null)
  
  all_emp_p <- matrix(NA, ncol = 1, nrow = length(m))
  high_m_locale <- which(m > 0)
  low_m_locale <- which(m < 0)
  
  low_p <- 1-low_fun(abs(m[low_m_locale]))
  high_p <- 1-high_fun(abs(m[high_m_locale]))
  
  all_emp_p[low_m_locale,1] <- low_p
  all_emp_p[high_m_locale,1] <- high_p
  # hist(all_emp_p)
  
  m12 <- matrix(c(marker_mat[,2],marker_mat[,1],as.numeric(as.matrix(influences_org[,3])),as.numeric(as.matrix(influences_org[,4])),(abs(as.numeric(influences_org[,3])) / as.numeric(influences_org[,4]))), ncol = 5)	
  colnames(m12) <- c("Source","Target","Effect","SE","|Effect|/SE")	
  
  m21 <- matrix(c(marker_mat[,1],marker_mat[,2],as.numeric(as.matrix(influences_org[,5])),as.numeric(as.matrix(influences_org[,6])),(abs(as.numeric(influences_org[,5])) / as.numeric(influences_org[,6]))), ncol = 5)
  colnames(m21) <- c("Source","Target","Effect","SE","|Effect|/SE")
  
  p_adjusted <- p.adjust(all_emp_p, method = pval_correction)
  final_table <- rbind(m12, m21)
  final_table <- cbind(final_table, all_emp_p, p_adjusted)
  colnames(final_table)[c(6,7)] <- c("P_empirical", "p_adjusted")
  final_table <- final_table[order(as.numeric(final_table[,"|Effect|/SE"]), decreasing = TRUE),]
  
  data_obj$var_to_var_p_val <- final_table
  
  return(data_obj)
}
