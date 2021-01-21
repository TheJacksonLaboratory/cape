#' Checks for unused markers
#' 
#' This function checks for markers that aren't used
#' in cape. For example, markers on sex chromosomes,
#' or mitochondrial markers. It also removes any 
#' markers that are invariant across all individuals.
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object
#' @keywords internal
#'


check_bad_markers <- function(data_obj, geno_obj){
  
  remove_flag = 0
  gene <- get_geno(data_obj, geno_obj)
  
  #we do not scan markers on the sex chromosomes
  #take these out here.
  x_locale <- grep("X", data_obj$chromosome, ignore.case = TRUE)
  y_locale <- grep("Y", data_obj$chromosome, ignore.case = TRUE)
  
  
  #take out markers with missing data
  # num_allele <- apply(gene, locus_dim, function(x) colSums(x, na.rm = TRUE))
  num_allele <- apply(gene, 3, function(x) sum(x, na.rm = TRUE))
  mono_allele <- which(num_allele == 0)

  if(length(x_locale)){
    remove_flag = 1
  }
  
  if(length(y_locale) > 0){
    remove_flag = 1
  }
  if(length(mono_allele) > 0){
    remove_flag = 1
  }

  if(remove_flag == 1){
    cross_obj <- remove_unused_markers(data_obj, geno_obj)
  }else{
    cross_obj <- list(data_obj, geno_obj)
  }

  return(cross_obj)
}