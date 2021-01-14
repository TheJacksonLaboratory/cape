#' Creates a covariate from a genetic marker
#' 
#' Occasionally, researchers may want to condition 
#' marker effects on another genetic marker. For example,
#' the HLA locus in humans has very strong effects on 
#' immune phenotypes, and can swamp smaller effects from
#' other markers. It can be helpful to condition on markers
#' in the HLA region to find genetic modifiers of these 
#' markers.
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object
#' @param singlescan_obj It is possible to automatically identify
#' markers to use as covariates based on their large main effects.
#' If this is desired, a singlescan object is required.
#' @param covar_thresh If designating markers as covariates based
#' on their main effect size is desired, the covar_thresh indicates
#' the main effect size above which a marker is designated as a 
#' covariate.
#' @param markers Marker covariates can also be designated manually.
#' markers takes in a vector of marker names or numbers and assigns
#' the designated markers as covariates.
#'
#' @return This function returns the data object with additional 
#' information specifying which markers are to be used as covariates.
#' this information can be retrieved with \code{\link{get_covar}}.
#' 
#' @examples 
#' \dontrun{
#' #convert markers with effect sizes greater than 6 to covariates.
#' #this requires a singlescan_obj
#' data_obj <- marker2covar(data_obj, geno_obj, singlescan_obj, covar_thresh = 6)
#' 
#' #convert the first marker to a covariate
#' #this does not require a singlescan_obj
#' marker_name <- dimnames(geno_obj)[[3]][1]
#' data_obj <- marker2covar(data_obj, geno_obj, markers = marker_name)
#' }
#' @seealso \code{\link{get_covar}}
#' @export

marker2covar <- function(data_obj, geno_obj, singlescan_obj = NULL, covar_thresh = NULL, markers = NULL){
  
  if(!is.null(covar_thresh)){
    oneD <- singlescan_obj$singlescan.results
    if(is.null(oneD)){stop("singlescan_obj is required if setting covariates by a t threshold.\n")}
  }
  
  geno_mat <- get_geno(data_obj, geno_obj)
  
  #if the user has specified a t threshold for covariate specification
  if(!is.null(covar_thresh)){
    
    marker_names <- data_obj$geno_names[[3]]
    
    singlescan_obj$covar_thresh <- covar_thresh
    
    covar_which <- lapply(oneD, function(x) which(x[,"t_stat"] >= covar_thresh))
    covar_names <- unique(unlist(lapply(covar_which, function(x) names(x))))
    new_covar_locale <- get_marker_idx(data_obj, covar_names)
    new_covar <- geno_mat[,new_covar_locale,drop=FALSE]
    dimnames(new_covar)[[3]] <- covar_names
    
    snp_names <- get_marker_name(data_obj, covar_names)
    
    g_covar_info <- rbind(snp_names, data_obj$chromosome[new_covar_locale], data_obj$marker_location[new_covar_locale])
    colnames(g_covar_info) <- data_obj$marker_num[new_covar_locale]
    rownames(g_covar_info) <- c("name", "chromosome", "position")
    
    data_obj <- remove_markers(data_obj, markers_to_remove = snp_names)
    data_obj$g_covar_table <- new_covar
    data_obj$g_covar <- g_covar_info
    return(data_obj)		
  } #end case for setting covariates by a threshold
  
  
  if(!is.null(markers)){
    marker_names <- data_obj$geno_names[[3]]
    
    marker_locale <- get_marker_idx(data_obj, markers)		
    new_covar <- geno_mat[,,marker_locale,drop=FALSE]
    ref_locale <- which(dimnames(new_covar)[[2]] == data_obj$ref_allele)
    the_rest <- setdiff(1:dim(new_covar)[[2]], ref_locale)
    covar_mat <- new_covar[,the_rest,]
    data_obj$g_covar_table <- covar_mat
    
    g_covar_info <- rbind(marker_names[marker_locale], data_obj$chromosome[marker_locale], data_obj$marker_location[marker_locale])
    colnames(g_covar_info) <- data_obj$marker_num[marker_locale]
    rownames(g_covar_info) <- c("name", "chromosome", "position")
    
    data_obj$g_covar <- g_covar_info
    data_obj <- remove_markers(data_obj, markers)
    
  }
  
  return(data_obj)
  
}