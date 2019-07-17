#' Select SNPs for pairscan based on a gene list
#'
#' This function selects SNPs to test for the CAPE
#' pairscan based on an ordered gene list
#' It 	goes down list of genes from top to bottom 
#' pulling out SNPs from around the named genes 
#' until we accumulate the target number of SNPs
#' This function also adds the distribution of 
#' effect sizes to the data.obj, so you can select
#' by this distribution for the pairscan null
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object
#' @param gene.list an array of Ensembl gene names
#' @param bp.buffer default = 1000
#' @param organism "mouse" or "human"
#' 
select.markers.for.pairscan.by.gene <- function(data.obj, geno.obj, gene.list, bp.buffer = 1000, organism = c("mouse", "human")){
  
  # install biomaRt if it's not installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  if (!("biomaRt" %in% (.packages()))) {
    BiocManager::install("biomaRt", update=FALSE, ask=FALSE)
    library(biomaRt)
  }
  
  ref.allele <- data.obj$ref_allele
  organism <- organism[1]
  
  data.obj$marker.selection.method <- "by.gene"
  
  if(length(data.obj$chromosome) != length(data.obj$geno_names[[3]])){
    stop("There are unequal lengths in the marker data\n")
  }
  
  #different archives work at different times. Use the most up to date that's actually working
  all.var <- ls(globalenv())
  if(length(which(all.var == "lib")) == 0){
    if(organism[1] == "mouse"){		
      # lib <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
      lib <<- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "may2017.archive.ensembl.org")
    }else{
      # lib <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")	
      lib <<- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "may2017.archive.ensembl.org")
    }
  }
  
  gene.info <- biomaRt::getBM(c("external_gene_name", "chromosome_name", "start_position", "end_position"), "external_gene_name", gene.list, lib)	
  #take out the scaffold chromosomes
  gene.info <- gene.info[which(!suppressWarnings(is.na(as.numeric(gene.info[,"chromosome_name"])))),]
  
  snp.list <- NULL
  for (gene.name in gene.list) {
    gene.locale <- which(tolower(gene.info[,1]) == tolower(gene.name))
    gene.chr <- gene.info[gene.locale,"chromosome_name"]
    gene.start <- gene.info[gene.locale,"start_position"] - bp.buffer
    gene.end <- gene.info[gene.locale,"end_position"] + bp.buffer
    
    chr.locale <- which(data.obj$chromosome == gene.chr)
    after.start <- which(data.obj$marker_location >= gene.start)
    before.end <- which(data.obj$marker_location <= gene.end)
    snp.pos <- Reduce("intersect", list(chr.locale, after.start, before.end))
    
    snp.list <- c(snp.list, data.obj$geno_names[[3]][snp.pos])
  }
  
  if (length(snp.list) == 0) {
    stop("The genes in the supplied gene list are not found either in the data.obj or in Ensembl.\n")
  }
  
  geno <- get.geno(data.obj, geno.obj)
  snp.locale <- match(snp.list, dimnames(geno)[[3]])
  geno.for.pairscan <- geno[,2,snp.locale]
  alt.allele <- setdiff(dimnames(geno)[[2]], ref.allele)
  colnames(geno.for.pairscan) <- paste0(colnames(geno.for.pairscan), "_", alt.allele)
  data.obj$geno_for_pairscan <- geno.for.pairscan
  
  data.obj$organism <- organism
  data.obj$bp_buffer <- bp.buffer
  data.obj$marker_selection_method = "by.gene"	
  return(data.obj)
  
}