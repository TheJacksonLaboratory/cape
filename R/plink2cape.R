#' Convert plink2 files to cape format
#'
#' @param ped full path to the ped file
#' @param map full path to the map file
#' @param pheno full path to the pheno file
#' @param out full path to the output file
#' @param missing_genotype default is "0"
#' @param no_fid boolean, default is FALSE
#' @param no_parents boolean, default is FALSE
#' @param no_sex boolean, default is FALSE
#' @param no_pheno boolean, default is FALSE
#' @param verbose boolean, default is FALSE, gives some happy little progress messages
#' @param overwrite boolean, default is FALSE, will only remove the existing file if this is set to TRUE
#'
#' @details For further information about PLINK and its file formats,
#' see \url{https://zzz.bwh.harvard.edu/plink/}
#' 
#' @return A list with two elements: data_obj and geno_obj
#' These objects are formatted for use in cape and must then
#' be separated to use in \code{\link{run_cape}}.
#'
#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, 
#' Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007) 
#' PLINK: a toolset for whole-genome association and population-based 
#' linkage analysis. American Journal of Human Genetics, 81.
#' 
#' @examples 
#' \dontrun{
#' #convert files with default names to a data_obj
#' data_obj <- plink2cape()
#' }
#' @export
plink2cape <- function(ped = "test.ped", map = "test.map", pheno = "test.pheno", 
                       out = "out.csv", missing_genotype = "0", no_fid = FALSE,
                       no_parents = FALSE, no_sex = FALSE, no_pheno = FALSE,
                       verbose = FALSE, overwrite = FALSE){
  
  # first check to see if the output file already exists
  if (file.exists(out) & overwrite == TRUE) {
    file.remove(out)
  } else if (file.exists(out) & overwrite == FALSE) {
    stop(paste('The file', out, 'already exists. Please rename it or set overwrite to TRUE.', collapse = ' '))
  }
  
  pheno_table <- as.matrix(read.table(pheno, header = TRUE))
  # the map file should contain the following columns (see: https://www.cog-genomics.org/plink2/formats#map):
  # 1. Chromosome code. PLINK 1.9 also permits contig names here, but most older programs do not.
  # 2. Variant identifier
  # 3. Position in morgans or centimorgans (optional; also safe to use dummy value of '0')
  # 4. Base-pair coordinate
  # All lines must have the same number of columns 
  # (so either no lines contain the morgans/centimorgans column, or all of them do).
  mapdata <- read.table(map, colClasses=c("character"))
  colnames(mapdata) <- c("Chr", "ID", "cM", "BP")
  
  # create column names for the SNPs
  SNPcolnames <- paste0(unlist(lapply(mapdata[,"ID"],rep,2)), c(".A",".B"))
  
  # The ped file should contain the following (see: https://www.cog-genomics.org/plink2/formats#ped):
  # The ped file contains pedigree and genotype information
  # It has no header line, and one line per sample with 2V+6 fields where V is the number of variants. 
  # 
  # 1. Family ID ('FID')
  # 2. Within-family ID ('IID'; cannot be '0')
  # 3. Within-family ID of father ('0' if father isn't in dataset)
  # 4. Within-family ID of mother ('0' if mother isn't in dataset)
  # 5. Sex code ('1' = male, '2' = female, '0' = unknown)
  # 6. Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
  # 7 & 8. allele calls for the first variant in the .map file ('0' = no call)
  # 9 & 10. allele calls for the second variant in the .map file ('0' = no call)
  # 
  peddata <- scan(ped, what=character(), na.strings=c("-9"))
  columnNames <- NULL
  columnNames <- c(columnNames, "IID")
  if(!no_fid) columnNames <- c(columnNames, "FID")
  if(!no_parents) columnNames <- c(columnNames, "PID", "MID")
  if(!no_sex) columnNames <- c(columnNames, "Sex")
  if(!no_pheno) columnNames <- c(columnNames, "Pheno")
  
  peddata <- matrix(peddata, ncol=length(c(columnNames, SNPcolnames)), byrow = TRUE)
  
  colnames(peddata) <- c(columnNames, SNPcolnames)
  
  # If there is no phenotype, create a random one
  if(no_pheno) {
    peddata <- cbind(peddata, Pheno=runif(nrow(peddata)))
  }
  
  # If there is no sex, then everyone is a female
  if(no_sex) {
    peddata <- cbind(peddata, Sex=rep(0, nrow(peddata)))
  }
  
  # R/qtl uses m and f for males and females, so let's map those
  peddata[peddata[,"Sex"] == 1, "Sex"] <- "m"; peddata[peddata[,"Sex"] == 2, "Sex"] <- "f"
  
  # Start with an empty genotype matrix
  genotypes <- matrix(NA, length(mapdata[,"ID"]), nrow(peddata))
  rownames(genotypes) <- mapdata[,"ID"]
  column <- length(columnNames)+1
  for (snp in mapdata[,"ID"]) {
    # Get the columns associated with this SNP
    cols <- column:(column+1)
    # The SNP alleles
    snpalleles <- sort(unique(unlist(as.character(peddata[,cols]))))
    # Missing data should not count as an allele
    if (missing_genotype %in% snpalleles) {
      snpalleles <- snpalleles[-which(snpalleles == missing_genotype)]
    }
    if (length(snpalleles) > 2) {
      warning("WARNING", snp, "found multi allelic marker:", snpalleles, ", passed as all missing\n")
      genotype <- rep(NA, nrow(peddata))
    } else {
      # a bit of debugging info if you want some
      if (verbose) {
        cat((column - (length(columnNames)+1)) / 2,"/", length(mapdata[,"ID"]), snp,"found", snpalleles,"\n")
      }
      genotype <- apply(peddata[,cols], 1, function(x) {
        # if missing genotype data
        if (x[1] == missing_genotype) return(NA)
        # if heterozygous
        if (x[1] != x[2]) return(2)
        # if homozygous, allele 1/A
        if (x[1] == x[2] && x[1] == snpalleles[1]) return(1)
        # if homozygous, allele 2/B
        if (x[1] == x[2] && x[1] == snpalleles[2]) return(3)
      })
    }
    column <- column + 2
    genotypes[snp,] <- genotype
  }
  
  #create a table in the the qtl csv format
  geno_info <- t(as.matrix(cbind(mapdata[,c("ID","Chr","BP")])))
  geno_mat <- rbind(geno_info, t(genotypes)-1)
  pheno_padding <- matrix(NA, nrow = 3, ncol = ncol(pheno_table)-2)
  pheno_padding[1,] <- colnames(pheno_table)[3:ncol(pheno_table)]
  pheno_mat <- rbind(pheno_padding, pheno_table[,3:ncol(pheno_table)])
  final_table <- cbind(pheno_mat, geno_mat)
  write.table(final_table, out, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE,
  na = "-")
  
	cross_obj <- read_population(out, verbose = verbose)
 	new_obj <- cape2mpp(cross_obj)
 	return(new_obj)

}

