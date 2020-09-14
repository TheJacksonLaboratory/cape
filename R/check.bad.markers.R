#' Checks for unused markers
#' 
#' This function checks for markers that aren't used
#' in cape. For example, markers on sex chromosomes,
#' or mitochonrdial markers. It also removes any 
#' markers that are invariant across all individuals.
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj
#'


check.bad.markers <- function(data.obj, geno.obj){
  
  remove.flag = 0
  gene <- get.geno(data.obj, geno.obj)
  
  #we do not scan markers on the sex chromosomes
  #take these out here.
  x.locale <- grep("X", data.obj$chromosome, ignore.case = TRUE)
  y.locale <- grep("Y", data.obj$chromosome, ignore.case = TRUE)
  
  
  #take out markers with missing data
  # num.allele <- apply(gene, locus.dim, function(x) colSums(x, na.rm = TRUE))
  num.allele <- apply(gene, 3, function(x) sum(x, na.rm = TRUE))
  mono.allele <- which(num.allele == 0)

  if(length(x.locale)){
    remove.flag = 1
  }
  
  if(length(y.locale) > 0){
    remove.flag = 1
  }
  if(length(mono.allele) > 0){
    remove.flag = 1
  }

  if(remove.flag == 1){
    cross.obj <- remove.unused.markers(data.obj, geno.obj)
  }else{
    cross.obj <- list(data.obj, geno.obj)
  }

  return(cross.obj)
}