#' Write significant cape interactions to a csv file
#' 
#' This function takes in the final data object and 
#' writes the variant influences that are at or below 
#' the specified significance level.
#' 
#' The columns of the output file are the following:
#'  Source: The marker that is the source of the directed interaction
#'  Chr: The chromosome on which the source marker lives
#'  Position: The genomic position of the source marker
#'  Target: If the effect is an interaction, this column 
#'    lists the marker that is the target of the directed interaction.
#'    If the effect is a main effect, this column lists the 
#'    trait that is the target of the main effect.
#'  Chr: The chromosome on which the target marker lives.
#'    If the effect is a main effect, this is listed as 0.
#'  Position: The genomic position of the target marker. If
#'    the effect is a main effect, this is listed as 1.
#'  Conditioning: If the effect is a main effect, this column identifies
#'    the marker on which the main effect marker was conditioned when
#'    it had it's largest main effect.
#'  Chr: If the effect is a main effect, this column lists the chromosome
#'    on which the conditioning marker lives
#'  Position: If the effect is a main effect, this column lists the
#'    genomic position of the conditioning marker.
#'  Effect: The effect size of the effect, either main effect or interaction.
#'  SE: The standard error of the effect, either main effect or interaction.
#'  |Effect|/SE: The standardized effect
#'  P_empirical: The empirical p value calculated from permutations
#'  p_adjusted: The p value adjusted by the method specified in the parameter file.
#' @param data_obj a \code{\link{Cape}} object
#' @param p_or_q A threshold indicating the maximum adjusted p value considered 
#'   significant. If an FDR method has been used to correct for multiple testing, 
#'   this value specifies the maximum q value considered significant.
#' @param include_main_effects Whether to include main effects (TRUE) or only
#'    interaction effects (FALSE) in the output table.
#' @param filename A character vector specifying the name of the file.
#' @param delim A character string indicating the delimiter in the data file. 
#'    The default indicates a comma-separated file (",").
#' @param mark_covar A logical value. If TRUE, an asterisk is appended the 
#'  names of markers used as covariates in the pair scan.
#' @param write_file A logical value indicating whether the table should be 
#'   written to a file or simply returned.
#'
#' @return If write_file is TRUE, this function writes the results table
#' to a file and invisibly returns the table. If write_file is FALSE, the
#' function returns the results table without writing to file.
#' 
#' @importFrom utils write.table
#' 
#' @examples 
#' \dontrun{
#' inf_table <- write_variant_influences(data_obj)
#' }
#' 
#' @export
write_variant_influences <- function(data_obj, p_or_q = 0.05, include_main_effects = TRUE, 
                                   filename = "Variant.Influences.csv", delim = ",", 
                                   mark_covar = FALSE, write_file = TRUE){
  
  
  var_influences <- data_obj$var_to_var_p_val
  pheno_results <- data_obj$max_var_to_pheno_influence
  pheno_names <- names(pheno_results)
  
  if(is.null(var_influences)){
    stop("calc_p() must be run to calculate variant-to-variant influences.")
  }
  
  
  if(is.null(pheno_results)){
    stop("direct_influence() must be run to calculate variant-to-trait influences.")
  }
  
  var_sig_col <- which(colnames(var_influences) == "p_adjusted")
  
  sig_var <- which(as.numeric(var_influences[, var_sig_col]) <= p_or_q)
  
  
  if(length(sig_var) > 0){
    var_table <- var_influences[sig_var,,drop=FALSE]
  }else{
    var_table <- NULL
  }
  
  
  add_name <- function(pheno_table, pheno_name){
    final_result <- cbind(pheno_table[,1,drop=FALSE], rep(pheno_name, nrow(pheno_table)), pheno_table[,2:ncol(pheno_table),drop=FALSE])
    return(final_result)
  }
  
  final_table = NULL
  if(include_main_effects){	
    #pull out the significan main effects and
    #add the phenotype name to the phenotype 
    #results, so we can put everything in one big table
    pheno_sig_col <- which(colnames(pheno_results[[1]]) == "p_adjusted")
    sig_pheno <- lapply(pheno_results, function(x) x[which(as.numeric(x[,pheno_sig_col]) <= p_or_q),,drop=FALSE])
    
    #add phenotype names
    sig_pheno_labels <- mapply(add_name, sig_pheno, pheno_names, SIMPLIFY = FALSE)
    pheno_table <- Reduce("rbind", sig_pheno_labels)
    orig_names <- colnames(pheno_results[[1]])
    new_names <- append(orig_names, "target", after = 1)
    colnames(pheno_table) <- new_names
    
    #take out the raw t statistic column
    remove_col <- which(colnames(pheno_table) == "t_stat")
    pheno_table <- pheno_table[,-remove_col,drop=FALSE]
    colnames(pheno_table) <- c("Source", "Target", "conditioning_marker", "Effect", "SE", "|Effect|/SE", "P_empirical", "p_adjusted")
    
    
    #add a column of NAs to the var_table to account for 
    #the conditioning_marker column in the pheno_table	
    if(length(var_table) > 0){
      conditioning_marker <- rep(NA, nrow(var_table))
      var_table <- cbind(var_table[,1:2,drop=FALSE], conditioning_marker, var_table[,3:ncol(var_table),drop=FALSE])
    }
    
    final_table <- rbind(var_table, pheno_table)
  }else{
    
    if(length(var_table) > 0){
      conditioning_marker <- rep(NA, nrow(var_table))
      final_table <- cbind(var_table[,1:2,drop=FALSE], conditioning_marker, var_table[,3:ncol(var_table),drop=FALSE])
    }
  }
  
  if(is.null(final_table)){
    full_table <- final_table <- matrix(NA, nrow = 1, ncol = 14)
    colnames(full_table) <- c("Source",	"Chr",	"Position",	"Target",	"Chr",	
    "Position",	"conditioning_marker",	"Chr",	"Position",	"Effect",	"SE",
    "|Effect|/SE",	"P_empirical", "p_adjusted")
    }else{
    final_table <- final_table[order(as.numeric(final_table[,"|Effect|/SE"]), decreasing = TRUE),,drop=FALSE]
    
    source_chr <- get_marker_chr(data_obj, final_table[,1])
    target_chr <- get_marker_chr(data_obj, final_table[,2])
    source_loc <- get_marker_location(data_obj, final_table[,1])
    target_loc <- get_marker_location(data_obj, final_table[,2])
    if(include_main_effects){
      cond_chr <- get_marker_chr(data_obj, final_table[,3])
      cond_loc <- get_marker_location(data_obj, markers = final_table[,3])
    }else{
      cond_chr <- rep(NA, nrow(final_table))
      cond_loc <- rep(NA, nrow(final_table))
    }
    
    full_table <- cbind(final_table[,1,drop=FALSE], source_chr, source_loc, final_table[,2,drop=FALSE], target_chr, target_loc, final_table[,3,drop=FALSE], cond_chr, cond_loc, final_table[,4:ncol(final_table),drop=FALSE])
  }
  
  if(mark_covar){
    covar_info <- get_covar(data_obj)
    covar_names <- covar_info$covar_names
    covar_source_locale <- which(full_table[,1] %in% covar_names)
    covar_target_locale <- which(full_table[,4] %in% covar_names)
    if(length(covar_source_locale) > 0){
      full_table[covar_source_locale,1] <- paste(full_table[covar_source_locale,1], "*", sep = "")
    }
    if(length(covar_target_locale) > 0){
      full_table[covar_target_locale,4] <- paste(full_table[covar_target_locale,4], "*", sep = "")
    }
  }
  
  colnames(full_table)[1:9] <- c("Source", "Chr", "Position", "Target", "Chr", "Position", "conditioning_marker", "Chr", "Position")
  if(write_file){
    write.table(full_table, file = filename, quote = FALSE, sep = delim, row.names = FALSE)	
    invisible(full_table)
  }else{
    return(full_table)
  }
}