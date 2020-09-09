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
#'  Chr: The chomosome on which the target marker lives.
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
#'  p.adjusted: The p value adjusted by the method specified in the parameter file.
#' @param data.obj a \code{\link{Cape}} object
#' @param p.or.q A threshold indicating the maximum adjusted p value considered 
#'   significant. If an FDR method has been used to correct for multiple testing, 
#'   this value specifies the maximum q value considered significant.
#' @param include.main.effects Whether to include main effects (TRUE) or only
#'    interaction effects (FALSE) in the output table.
#' @param filename A character vector specifying the name of the file.
#' @param delim A character string indicating the delimeter in the data file. 
#'    The default indicates a comma-separated file (",").
#' @param mark.covar A logical value. If TRUE, an asterisk is appended the 
#'  names of markers used as covariates in the pair scan.
#' @param write.file A logical value indicating whether the table should be 
#'   written to a file or simply returned.
#'
#' @return If write.file is TRUE, this function writes the results table
#' to a file and invisibly returns the table. If write.file is FALSE, the
#' function returns the results table without writing to file.
#' 
#' @export
#' 
writeVariantInfluences <- function(data.obj, p.or.q = 0.05, include.main.effects = TRUE, 
                                   filename = "Variant.Influences.csv", delim = ",", 
                                   mark.covar = FALSE, write.file = TRUE){
  
  
  var.influences <- data.obj$var_to_var_p_val
  pheno.results <- data.obj$max_var_to_pheno_influence
  pheno.names <- names(pheno.results)
  
  if(is.null(var.influences)){
    stop("calc.p() must be run to calculate variant-to-variant influences.")
  }
  
  
  if(is.null(pheno.results)){
    stop("direct.influence() must be run to calculate variant-to-trait influences.")
  }
  
  var.sig.col <- which(colnames(var.influences) == "p.adjusted")
  
  sig.var <- which(as.numeric(var.influences[, var.sig.col]) <= p.or.q)
  
  
  if(length(sig.var) > 0){
    var.table <- var.influences[sig.var,,drop=FALSE]
  }else{
    var.table <- NULL
  }
  
  
  add.name <- function(pheno.table, pheno.name){
    final.result <- cbind(pheno.table[,1,drop=FALSE], rep(pheno.name, nrow(pheno.table)), pheno.table[,2:ncol(pheno.table),drop=FALSE])
    return(final.result)
  }
  
  final.table = NULL
  if(include.main.effects){	
    #pull out the significan main effects and
    #add the phenotype name to the phenotype 
    #results, so we can put everything in one big table
    pheno.sig.col <- which(colnames(pheno.results[[1]]) == "p.adjusted")
    sig.pheno <- lapply(pheno.results, function(x) x[which(as.numeric(x[,pheno.sig.col]) <= p.or.q),,drop=FALSE])
    
    #add phenotype names
    sig.pheno.labels <- mapply(add.name, sig.pheno, pheno.names, SIMPLIFY = FALSE)
    pheno.table <- Reduce("rbind", sig.pheno.labels)
    orig.names <- colnames(pheno.results[[1]])
    new.names <- append(orig.names, "target", after = 1)
    colnames(pheno.table) <- new.names
    
    #take out the raw t statistic column
    remove.col <- which(colnames(pheno.table) == "t.stat")
    pheno.table <- pheno.table[,-remove.col,drop=FALSE]
    colnames(pheno.table) <- c("Source", "Target", "Conditioning.Marker", "Effect", "SE", "|Effect|/SE", "P_empirical", "p.adjusted")
    
    
    #add a column of NAs to the var.table to account for 
    #the conditioning.marker column in the pheno.table	
    if(length(var.table) > 0){
      conditioning.marker <- rep(NA, nrow(var.table))
      var.table <- cbind(var.table[,1:2,drop=FALSE], conditioning.marker, var.table[,3:ncol(var.table),drop=FALSE])
    }
    
    final.table <- rbind(var.table, pheno.table)
  }else{
    
    if(length(var.table) > 0){
      conditioning.marker <- rep(NA, nrow(var.table))
      final.table <- cbind(var.table[,1:2,drop=FALSE], conditioning.marker, var.table[,3:ncol(var.table),drop=FALSE])
    }
  }
  
  if(is.null(final.table)){
    full.table <- final.table <- matrix(NA, nrow = 1, ncol = 14)
    colnames(full.table) <- c("Source",	"Chr",	"Position",	"Target",	"Chr",	
    "Position",	"Conditioning.Marker",	"Chr",	"Position",	"Effect",	"SE",
    "|Effect|/SE",	"P_empirical", "p.adjusted")
    }else{
    final.table <- final.table[order(as.numeric(final.table[,"|Effect|/SE"]), decreasing = TRUE),,drop=FALSE]
    
    source.chr <- get.marker.chr(data.obj, final.table[,1])
    target.chr <- get.marker.chr(data.obj, final.table[,2])
    source.loc <- get.marker.location(data.obj, final.table[,1])
    target.loc <- get.marker.location(data.obj, final.table[,2])
    if(include.main.effects){
      cond.chr <- get.marker.chr(data.obj, final.table[,3])
      cond.loc <- get.marker.location(data.obj, markers = final.table[,3])
    }else{
      cond.chr <- rep(NA, nrow(final.table))
      cond.loc <- rep(NA, nrow(final.table))
    }
    
    full.table <- cbind(final.table[,1,drop=FALSE], source.chr, source.loc, final.table[,2,drop=FALSE], target.chr, target.loc, final.table[,3,drop=FALSE], cond.chr, cond.loc, final.table[,4:ncol(final.table),drop=FALSE])
  }
  
  if(mark.covar){
    covar.info <- get.covar(data.obj)
    covar.names <- covar.info$covar.names
    covar.source.locale <- which(full.table[,1] %in% covar.names)
    covar.target.locale <- which(full.table[,4] %in% covar.names)
    if(length(covar.source.locale) > 0){
      full.table[covar.source.locale,1] <- paste(full.table[covar.source.locale,1], "*", sep = "")
    }
    if(length(covar.target.locale) > 0){
      full.table[covar.target.locale,4] <- paste(full.table[covar.target.locale,4], "*", sep = "")
    }
  }
  
  colnames(full.table)[1:9] <- c("Source", "Chr", "Position", "Target", "Chr", "Position", "Conditioning.Marker", "Chr", "Position")
  if(write.file){
    write.table(full.table, file = filename, quote = FALSE, sep = delim, row.names = FALSE)	
    invisible(full.table)
  }else{
    return(full.table)
  }
}