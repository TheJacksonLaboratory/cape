#' Read the parameter file, add missing entries
#' 
#' This function returns reads in the YAML file
#' and checks for any parameters that might not
#' be included. This may not matter for the given
#' run, but it's handy to be able to check for 
#' any and all potential variables.
#'
#' @param filename full path to the .yml file holding CAPE parameters
#'
#' @return Returns a named list with all possible options
#' 
#' @export
read.parameters <- function(filename = "cape.parameters.yml"){
  
  parameter.table <- yaml::read_yaml(filename)
  
  #================================================
  # general parameters
  #================================================
  gen.param <- c("traits", "covariates", "marker.covariates", "traits.scaled", "traits.normalized", "scan.what", "eig.which", "use.kinship", "kinship.type", "locus", "pop", "pval.correction")
  
  #================================================
  #single scan parameters
  #================================================
  single.param <- c("ref.allele", "singlescan.perm")
  
  #================================================
  # marker selection
  #================================================
  marker.param <- c("marker.selection.method", "SNPfile", "peak.density", "tolerance", "window.size", "num.alleles.in.pairscan", "bp.buffer", "organism")
  
  #================================================
  # pair scan
  #================================================
  pair.param <- c("max.pair.cor", "min.per.geno", "pairscan.null.size")
  
  for ( param.name in c(gen.param, single.param, marker.param, pair.param) ) {
    if ( !is.element(param.name, names(parameter.table)) ){
      parameter.table[param.name] = NULL
    }
  }
  
  return(parameter.table)
  
}
