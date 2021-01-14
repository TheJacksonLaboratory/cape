#' Get marker genomic position
#' 
#' Given a vector of marker names or numbers, this 
#' function returns the genomic coordinates for
#' each marker, not including the chromosome number,
#' which is retrieved using \code{\link{get_marker_chr}}.
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param markers A vector of marker names 
#' 
#' @return A vector the same length as the input markers vector
#' indicating the genomic coordinate of each marker.
#'
#' @examples 
#' \dontrun{
#' marker_names <- dimnames(geno_obj)[[3]]
#' marker_loc <- get_marker_location(data_obj, marker_names)
#' }
#'
#' @export

get_marker_location <- function(data_obj, markers){
  
  
  markers_with_na <- markers
  not_na_pos <- which(!is.na(markers_with_na))
  markers <- markers[not_na_pos]
  
  
  und_check <- grep("_", markers)
  if(length(und_check) > 0){
    markers <- sapply(strsplit(markers, "_"), function(x) x[1])
  }
  
  is_char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))
  
  if(is_char){
    marker_loc <- data_obj$marker_location[match(markers, data_obj$geno_names[[3]])]
    na_locale <- which(is.na(marker_loc))
    if(length(na_locale) > 0){
      covar_info <- get_covar(data_obj)
      geno_covar <- which(covar_info$covar_type == "g")
      if(length(geno_covar) > 0){
        geno_covar_locale <- match(markers[na_locale], covar_info$covar_names[geno_covar])
        geno_covar_loc <- data_obj$g_covar[3,geno_covar_locale]
        marker_loc[na_locale] <- geno_covar_loc
      }
    }
    na_locale <- which(is.na(marker_loc))
    #if there are still missing values, these are 
    #phenotypic covariates. They get dummy positions
    if(length(na_locale) > 0){ 
      pheno_covar <- which(covar_info$covar_type == "p")
      if(length(pheno_covar) > 0){
        marker_loc[na_locale] <- 1:length(na_locale)
      }
    } #end case for if there are still missing markers
    final_marker_pos <- rep(NA, length(markers_with_na))
    final_marker_pos[not_na_pos] <- marker_loc
    return(final_marker_pos)
  }
  
  
  if(!is_char){
    marker_loc <- data_obj$marker_location[match(markers, data_obj$marker_num)]
    na_locale <- which(is.na(marker_loc))
    if(length(na_locale) > 0){
      covar_info <- get_covar(data_obj)
      covar_locale <- which(!is.na(match(colnames(covar_info$covar_table), markers)))
      if(length(covar_locale) == 0){
        marker_loc[na_locale] <- match(markers[na_locale], covar_info$covar_names)
      }else{
        marker_loc[na_locale] <- match(markers[na_locale], colnames(covar_info$covar_table))
      }
    }
    final_marker_pos <- rep(NA, length(markers_with_na))
    final_marker_pos[not_na_pos] <- marker_loc
    return(final_marker_pos)
  }
  
}
