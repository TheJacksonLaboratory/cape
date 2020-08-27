#' Save the cross data in R/qtl CSV format
#' 
#' This function writes out a cape object
#' in a csv format readable both by \code{\link{read.population}}
#' in Cape or by read.cross in R/qtl.
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object
#' @param ref.allele a character, e.g., "A", that represents the reference allele in the data object
#' @param na either NA or a character used to represent a missing data value in the output
#' @param filename absolute or relative path to the output file's destination
#'
#' @references Broman et al. (2003) R/qtl: QTL mapping in experimental crosses. 
#' Bioinformatics 19:889-890 doi:10.1093/bioinformatics/btg112
#' 
#' @return Writes a file to the destination path
#'
#' @export

writePopulation <- function(data.obj, geno.obj, ref.allele = "A", na = NA, 
filename = "capeData.csv"){
  
  if(class(geno.obj)=="array"){
    geno.obj<-geno.obj}
  else{geno.obj<-geno.obj$geno}
  
  if(dim(geno.obj)[2] > 2){
    stop("writePopulation only works for two-parent populations.")
  }
  
  ref.allele <- ref.allele[1]
  alt.allele <- setdiff(colnames(geno.obj), ref.allele)
  geno <- get.geno(data.obj, geno.obj)
  geno.mat <- geno[,which(colnames(geno) == alt.allele),]
  pheno <- data.obj$pheno
  chr <- data.obj$chromosome
  loc <- data.obj$marker_location
  marker.names <- data.obj$geno_names[[3]]
  
  #build a matrix out of the object so we can write it out to a csv file
  final.geno <- rbind(chr, loc, geno.mat)
  final.geno[which(is.na(final.geno))] <- na
  colnames(final.geno) <- marker.names
  
  pheno[which(is.na(pheno))] <- na
  pheno.padding <- matrix(NA, nrow = 2, ncol = dim(pheno)[2])
  final.pheno <- rbind(pheno.padding, pheno)
  
  final.data <- cbind(final.pheno, final.geno)
  write.table(final.data, file = filename, sep = ",", quote = FALSE, row.names = FALSE, na = as.character(na))
  
}
