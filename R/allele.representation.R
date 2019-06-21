#' Finds alleles from each founder were seleted for the pairscan
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param color.scheme string "CC/DO" or "other"
#' @param pdf.label string
#'
#' @export
allele.representation <- function(data.obj, color.scheme = c("CC/DO", "other"), pdf.label = "Allele.Representation.pdf"){
  
  
  get.alleles <- function(set){
    split.alleles <- strsplit(set, "_")
    alleles <- rep(NA, length(split.alleles))
    for(i in 1:length(alleles)){
      if(length(split.alleles[[i]]) == 1){
        alleles[i] <- NA
      }else{
        start.pos = 2
        alleles[i] <- paste(split.alleles[[i]][start.pos:length(split.alleles[[i]])], collapse = ",")	
      }
    }
    return(alleles)
  }
  
  selected.alleles <- colnames(data.obj$geno_for_pairscan)
  just.markers <- sapply(strsplit(selected.alleles, "_"), function(x) x[1])
  alleles <- get.alleles(selected.alleles)
  
  allele.chr <- get.marker.chr(data.obj, just.markers)
  
  allele.colors <- get.allele.colors(color.scheme, unique(alleles))
  
  allele.labels <- allele.colors[,1]
  allele.abbr <- allele.colors[,2]
  
  #count number of chromosomes per allele
  chr.count <- sort(sapply(allele.abbr, function(x) length(unique(allele.chr[which(alleles == x)]))), decreasing = TRUE)
  
  allele.count <- sort(table(alleles), decreasing = TRUE)
  
  # TODO move this into the CAPE methods
  pdf(pdf.label, width = 10, height = 5)
  
  allele.locale <- match(names(allele.count), allele.abbr)
  
  a <- barplot(allele.count, col = allele.colors[allele.locale,3], names = allele.labels[allele.locale], ylim = c(0, max(allele.count)+14), cex.names = 2)
  text(x = a, y = allele.count+8, labels = allele.count, cex = 2)
  
  a <- barplot(chr.count, col = allele.colors[allele.locale,3], names = allele.labels[allele.locale], ylim = c(0, max(chr.count)+2), cex.names = 2)
  text(x = a, y = chr.count+1, labels = chr.count, cex = 2)
  
  dev.off()
  
}