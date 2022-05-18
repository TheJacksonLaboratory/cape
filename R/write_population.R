#' Save the cross data in R/qtl CSV format
#' 
#' This function writes out a cape object
#' in a csv format readable both by \code{\link{read_population}}
#' in Cape or by read.cross in R/qtl.
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object
#' @param ref_allele a character, e.g., "A", that represents the reference allele in the data object
#' @param na either NA or a character used to represent a missing data value in the output
#' @param filename absolute or relative path to the output file's destination
#'
#' @references Broman et al. (2003) R/qtl: QTL mapping in experimental crosses. 
#' Bioinformatics 19:889-890 doi:10.1093/bioinformatics/btg112
#' 
#' @return Writes a file to the destination path
#' 
#' @importFrom utils write.table
#' 
#' @examples 
#' \dontrun{
#' write_population(data_obj, geno_obj)
#' }
#' 
#' @export

write_population <- function(data_obj, geno_obj, ref_allele = "A", na = NA, 
  filename = "capeData.csv"){
  
  class_geno <- class(geno_obj)
  if(class_geno =="array"){
    geno_obj<-geno_obj}
  else{geno_obj<-geno_obj$geno}
  
  if(dim(geno_obj)[2] > 2){
    stop("write_population only works for two-parent populations.")
  }
  
  ref_allele <- ref_allele[1]
  alt_allele <- setdiff(colnames(geno_obj), ref_allele)
  geno <- get_geno(data_obj, geno_obj)
  geno_mat <- geno[,which(colnames(geno) == alt_allele),]
  pheno <- data_obj$pheno
  chr <- data_obj$chromosome
  loc <- data_obj$marker_location
  marker_names <- data_obj$geno_names[[3]]
  
  #build a matrix out of the object so we can write it out to a csv file
  final_geno <- rbind(chr, loc, geno_mat)
  final_geno[which(is.na(final_geno))] <- na
  colnames(final_geno) <- marker_names
  
  pheno[which(is.na(pheno))] <- na
  pheno_padding <- matrix(NA, nrow = 2, ncol = dim(pheno)[2])
  final_pheno <- rbind(pheno_padding, pheno)
  
  final_data <- cbind(final_pheno, final_geno)
  write.table(final_data, file = filename, sep = ",", quote = FALSE, row.names = FALSE, na = as.character(na))
  
}
