#' Get chromosome numbers for markers
#' 
#' Given a vector of marker names or numbers, this 
#' function returns the chromosome on which each 
#' marker lives.Covariates are assigned to chromosome 0.
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param markers A vector of marker names 
#' @param character_names A logical value indicating whether
#' the marker names are characters (TRUE) or numbers (FALSE)
#' 
#' @return A vector the same length as the input markers vector
#' indicating which chromosome each marker in markers lives on.
#'
#' @examples 
#' \dontrun{
#' marker_names <- dimnames(geno_obj)[[3]]
#' marker_chr <- get_marker_chr(data_obj, marker_names)
#' }
#' 
#' @keywords internal 
#' 
get_marker_chr <- function(data_obj, markers, character_names = TRUE){
  
  und_check <- grep("_", markers)
  if(length(und_check) > 0){
    markers <- sapply(strsplit(markers, "_"), function(x) x[1])
  }
  
  if(is.null(character_names)){
    is_char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))
  }
  if(character_names){
    is_char <- TRUE
  }
  if(!character_names){
    is_char <- FALSE
  }
  
  if(is_char){
    marker_chr <- data_obj$chromosome[match(markers, data_obj$geno_names[[3]])]
    na_locale <- which(is.na(marker_chr))
    if(length(na_locale) > 0){
      covar_info <- get_covar(data_obj)
      geno_covar <- which(covar_info$covar_type == "g")
      if(length(geno_covar) > 0){
        unassigned <- markers[na_locale]
        unassigned_locale <- match(unassigned, data_obj$g_covar[1,])
        geno_covar_chr <- data_obj$g_covar[2,unassigned_locale]
        marker_chr[na_locale] <- geno_covar_chr
      }
    }
    #The rest are phenotypic covariates which are assigned chr 0
    na_locale <- which(is.na(marker_chr))
    if(length(na_locale) > 0){
      marker_chr[na_locale] <- 0
    }
    return(marker_chr)
  }
  
  if(!is_char){
    marker_chr <- data_obj$chromosome[match(markers, data_obj$marker_num)]
    na_locale <- which(is.na(marker_chr))
    if(length(na_locale) > 0){
      covar_info <- get_covar(data_obj)
      geno_covar <- which(covar_info$covar_type == "g")
      if(length(geno_covar) > 0){
        unassigned <- markers[na_locale]
        unassigned_locale <- match(unassigned, colnames(data_obj$g_covar))
        geno_covar_chr <- data_obj$g_covar[2,unassigned_locale]
        marker_chr[na_locale] <- geno_covar_chr
      }
    }
    #The rest are phenotypic covariates which are assigned chr 0
    na_locale <- which(is.na(marker_chr))
    if(length(na_locale) > 0){
      marker_chr[na_locale] <- 0
    }
    return(marker_chr)
  }
}
