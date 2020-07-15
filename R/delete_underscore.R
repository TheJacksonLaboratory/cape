#' Delete underscores from marker names
#' 
#' Because cape uses underscores to separate marker names
#' from alleles, we need to delete any underscores that 
#' were present in the original marker names. This function
#' does that. 
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object
#' 
#' @return lists containing data.obj and geno.obj
#' with underscores removed from data.obj$geno_names
#' and dimnames(geno.obj). These two elements must be
#' entered separately into \code{link{run.cape}}.
#'
#' @export

delete_underscore <- function(data.obj, geno.obj = NULL){
  
  geno <- get.geno(data.obj, geno.obj)
  
  marker.names <- data.obj$geno_names[[3]]
  under.locale <- grep("_", marker.names)
  
  if(length(under.locale) > 0){
    bad.names <- marker.names[under.locale]
    new.names <- unlist(lapply(strsplit(bad.names, "_"), function(x) paste(x[1:length(x)], collapse = "")))
    
    data.obj$geno_names[[3]][under.locale] <- new.names
    dimnames(geno)[[3]][under.locale] <- new.names
    cat("Removing underscores from marker names\n")
  }	
  
  results <- list(data.obj, geno)
  names(results) <- c("data.obj", "geno.obj")
  
  return(results)
  
}

