#' Take out markers not used in cape
#' 
#' This function removes any markers that are not
#' used in cape. This includes markers on the sex 
#' chromosomes, mitochondrial markers, and any 
#' invariant markers.
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object
#' @param verbose A logical value indicating whether to print progress 
#' to the screen. Default is FALSE
#' 
#' @return an updated \code{\link{Cape}} object (data_obj)
#' 
#' @examples 
#' \dontrun{
#' data_obj <- remove_unused_markers(data_obj, geno_obj)
#' }
#'
#' @export

remove_unused_markers <- function(data_obj, geno_obj, verbose = FALSE){
  
  #we do not scan markers on the sex chromosomes
  #take these out here.
  x_locale <- grep("X", data_obj$chromosome, ignore.case = TRUE)
  if(length(x_locale) > 0){
    message("Removing markers on the X chromosome\n")
    data_obj$chromosome <- data_obj$chromosome[-x_locale]
    data_obj$marker_location <- data_obj$marker_location[-x_locale]
    data_obj$geno_names[[3]] <- data_obj$geno_names[[3]][-x_locale]
    data_obj$marker_num <- data_obj$marker_num[-x_locale]		
  }
  
  y_locale <- grep("Y", data_obj$chromosome, ignore.case = TRUE)
  if(length(y_locale) > 0){
    message("Removing markers on the Y chromosome\n")
    data_obj$chromosome <- data_obj$chromosome[-y_locale]
    data_obj$marker_location <- data_obj$marker_location[-y_locale]
    data_obj$marker_num <- data_obj$marker_num[-y_locale]	
    data_obj$geno_names[[3]] <- data_obj$geno_names[[3]][-y_locale]
  }
  
  m_locale <- grep("M", data_obj$chromosome, ignore.case = TRUE)
  if(length(m_locale) > 0){
    message("Removing markers on the Mitochondrial chromosome\n")
    data_obj$chromosome <- data_obj$chromosome[-m_locale]
    data_obj$marker_location <- data_obj$marker_location[-m_locale]
    data_obj$marker_num <- data_obj$marker_num[-m_locale]	
    data_obj$geno_names[[3]] <- data_obj$geno_names[[3]][-m_locale]
  }
  
  
  #take out markers with no variation
  gene <- get_geno(data_obj, geno_obj)
  if(verbose){cat("Checking for invariant markers.\n")}
  allelic_variation <- function(one_gene){
  	allele_var <- apply(one_gene, 2, function(x) var(x, na.rm = TRUE))
  	if(all(allele_var == 0)){
  		return(0)
  	}else{
  		return(1)
  	}
  }
  num_allele <- apply(gene, 3, allelic_variation)
  mono_allele <- which(num_allele == 0)
  if(length(mono_allele) > 0){
    message(paste("\nRemoving", length(mono_allele), "invariant markers:\n"))
    if(verbose){cat(data_obj$geno_names[[3]][mono_allele], sep = ", ")}
    data_obj$chromosome <- data_obj$chromosome[-mono_allele]
    data_obj$marker_location <- data_obj$marker_location[-mono_allele]
    data_obj$marker_num <- data_obj$marker_num[-mono_allele]
    data_obj$geno_names[[3]] <- data_obj$geno_names[[3]][-mono_allele]
  }
  
  return(data_obj)
  
}
