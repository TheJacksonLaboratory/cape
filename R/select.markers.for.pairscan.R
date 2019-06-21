#' This function takes in a singlescan object and bins the markers into bins based on the  either the LOD score, or the curve formed by the maximum allele t score at each point
#' 
#' This function takes in a singlescan object
#' and bins the markers into bins based on the 
#' either the LOD score, or the curve formed by 
#' the maximum allele t score at each point
#' num.alleles is the total number of alleles desired
#' peak.density is the density at which alleles will
#' be sampled within a peak. A value of 0.2 means about
#' 20% of the alleles will be selected, while a value of 
#' 0.5 means that 50% of alleles will be selected.
#' they are sampled uniformly at random
#' This function sets a cutoff such that you get the 
#' if you are specifying markers, you do not need to 
#' include a singlescan object.
#' Makes effects plots for multi-allelic 1D scans.
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param singlescan.obj a singlescan object
#' @param geno.obj a genotype object
#' @param specific.markers
#' @param num.alleles
#' @param peak.density default = 0.5
#' @param window.size
#' @param tolerance
#' @param plot.peaks
#' @param verbose
#' @param pdf.filename default = "Peak.Plots.pdf"
#' 
#' @export
select.markers.for.pairscan <- function(data.obj, singlescan.obj, geno.obj, specific.markers = NULL, 
                                        num.alleles = 1500, peak.density = 0.5, window.size = NULL, 
                                        tolerance = 5, plot.peaks = FALSE, verbose = FALSE, 
                                        pdf.filename = "Peak.Plots.pdf"){
  
  require(abind)
  chr <- unique(data.obj$chromosome)
  
  if(is.null(specific.markers)){
    data.obj$marker_selection_method <- "top.effects"
  }else{
    data.obj$marker_selection_method <- "from.list"	
  }
  
  cols <- c("grey", "white")	
  
  geno <- get.geno(data.obj, geno.obj)
  alleles <- dimnames(geno)[[2]]
  n.alleles <- length(alleles)
  
  if(class(singlescan.obj) == "list"){
    ref.allele <- singlescan.obj$ref.allele
    data.obj$ref_allele <- ref.allele
  }else{
    ref.allele <- data.obj$ref_allele
  }
  
  if(is.null(ref.allele)){
    allele.text <- paste(alleles, collapse = ", ")
    ref.allele <- readline(prompt = paste("Which allele do you want to use as the reference?\n", allele.text, "\n"))
    data.obj$ref_allele <- ref.allele
  }
  
  #if we are asking for more markers than there are in the dataset
  #just take all of them.
  if(num.alleles >= dim(geno)[3]*(dim(geno)[2]-1)){
    num.alleles = dim(geno)[3]
    data.obj$marker_selection_method <- "from.list"	
    alt.alleles <- setdiff(dimnames(geno)[[2]], ref.allele)
    specific.markers <- paste(dimnames(geno)[[3]], alt.alleles, sep = "_")
  }
  
  
  if(!is.null(specific.markers)){
    if(n.alleles == 2){
      ref.allele.locale <- which(data.obj$geno_names[[2]] == ref.allele)
      other.allele <- setdiff(1:2, ref.allele.locale)
      split.markers <- strsplit(as.character(specific.markers), "_")
      just.markers <- sapply(split.markers, function(x) x[1])
      geno.for.pairscan <- geno[,other.allele,just.markers]
      colnames(geno.for.pairscan) <- paste(colnames(geno.for.pairscan), data.obj$geno_names[[2]][other.allele], sep = "_")
      data.obj$geno_for_pairscan <- geno.for.pairscan
    }else{
      split.markers <- strsplit(specific.markers, "_")
      just.markers <- sapply(split.markers, function(x) x[1])
      just.alleles <- sapply(split.markers, function(x) x[2])
      geno.for.pairscan <- matrix(nrow = nrow(data.obj$pheno), ncol = length(just.markers))
      colnames(geno.for.pairscan) <- specific.markers
      for(i in 1:length(just.markers)){
        geno.for.pairscan[,i] <- geno[,just.alleles[i], just.markers[i]]
      }
    }
    
    if(verbose){cat("Removing markers that are not linearly independent...\n")}
    data.obj$geno_for_pairscan <- geno.for.pairscan
    geno.ind <- get.linearly.independent(data.obj, verbose = verbose)
    if(verbose){
      cat(length(geno.ind[[2]]), "allele(s) rejected.\n")
      cat("Final alleles selected:", "\t", ncol(geno.ind$independent.markers), "\n")
    }
    #we still need to specify values for selecting markers for the pairscan null distribution
    #so just use default value
    # cat("Generating the null distribution for the pairscan requires values for peak.density, window.size, and tolerance.\nSetting default values...")
    # data.obj$peak.density <- peak.density
    # data.obj$window.size <- window.size
    # data.obj$tolerance <- tolerance
    return(data.obj)
  }
  
  if(class(singlescan.obj) == "list"){ 
    results <- abs(singlescan.obj$singlescan.t.stats) #an actual singlescan object
  }else{
    results <- abs(singlescan.obj) #a singlescan matrix for calculating pairscan null distribution
  }
  
  
  filtered.results <- results
  
  covar.info <- get.covar(data.obj)
  results.no.covar <- results[which(!rownames(results) %in% covar.info$covar.names),,,drop=FALSE]
  result.chr <- get.marker.chr(data.obj, markers = rownames(results.no.covar), character.names = TRUE)
  
  
  sorted.results <- sort(abs(as.vector(results.no.covar)))
  #start with a cutoff that might be near the number of 
  #alleles desired
  guess.point <- num.alleles
  min.effect.size = min(tail(sorted.results, guess.point))
  
  #===============================================================
  #internal functions
  #===============================================================
  
  #how may peaks are above a given cutoff?
  num.peaks <- function(allele.curves, bins, cutoff){
    filtered.bins <- bins
    filtered.bins[which(abs(allele.curves) < cutoff)] <- NA
    num.peaks <- apply(filtered.bins, 2, function(x) length(unique(x))-1)
    return(num.peaks)
  }
  
  #how many alleles are above a given t stat cutoff
  num.markers <- function(allele.curves, cutoff){
    filtered.curves <- allele.curves
    filtered.curves[which(abs(allele.curves) < cutoff)] <- NA
    num.markers <- apply(filtered.curves, 2, function(x) length(which(!is.na(x))))
    return(num.markers)
  }
  
  #how many markers are in each peak at a given cutoff?
  markers.per.peak <- function(allele.curves, bins, cutoff){
    #make a results mat with enough columns for each allele
    #and enough rows for each bin. Each cell will count the
    #number of markers in the bin at the designated cutoff
    result.mat <- matrix(0, ncol = dim(allele.curves)[[3]], nrow = max(bins, na.rm = TRUE))
    rownames(result.mat) <- 1:nrow(result.mat)
    colnames(result.mat) <- dimnames(allele.curves)[[3]]
    
    #delete all effects that are less than the cutoff
    filtered.bins <- bins
    filtered.bins[which(abs(allele.curves) < cutoff)] <- NA
    #image(filtered.bins)
    #count the number of markers in each bin for each allele
    for(i in 1:ncol(filtered.bins)){
      counts <- table(filtered.bins[,i])
      result.mat[names(counts),i] <- counts
    }
    return(result.mat)
  }
  
  sample.peaks <- function(pheno.results, num.per.peak, bins){
    #sample markers based on peaks across all alleles
    sampled.markers <- vector(mode = "list", length = ncol(pheno.results))
    names(sampled.markers) <- colnames(pheno.results)
    for(i in 1:ncol(pheno.results)){
      #figure out which peaks we will sample from
      allele.markers <- NULL
      peaks.which <- which(num.per.peak[,i] > 0)
      if(length(peaks.which) > 0){
        for(j in 1:length(peaks.which)){
          #in each peak, pick the max, and sample the rest
          marker.locale <- which(bins[,i] == peaks.which[j])
          allele.markers <- c(allele.markers, marker.locale[which.max(abs(pheno.results[marker.locale,i]))])
          num.to.sample <- num.per.peak[peaks.which[j],i] - 1 #take off the maximum marker
          if(num.to.sample > 0){ #if there are still markers to get after grabbing the max
            unif.markers <- round(runif(num.to.sample, min = min(marker.locale), max = max(marker.locale)))
            allele.markers <- c(allele.markers, unif.markers)
          }#end case for sampling peak uniformly
        }#end looping through peaks for one allele
        sampled.markers[[i]] <- rownames(pheno.results)[sort(allele.markers)]
      }#end looping through alleles
    }
    return(sampled.markers)
  }
  
  allele.coin.flip <- function(alleles.per.bin, peak.density){
    one.locale <- which(alleles.per.bin == 1)
    coin.flip <- runif(length(one.locale), 0, 1)
    alleles.per.bin[which(coin.flip < peak.density)] <- 0
    alleles.per.bin[which(coin.flip >= peak.density)] <- 1/peak.density
    return(alleles.per.bin)
  }
  
  #===============================================================
  
  #===============================================================
  #group markers into bins based on their effect size profiles
  #===============================================================	
  num.pheno <- dim(results)[[2]]
  allele.bins <- vector(mode = "list", length = num.pheno)
  for(ph in 1:num.pheno){
    if(verbose){cat("\nBinning markers for", colnames(filtered.results)[ph], "\n")}
    pheno.results <- results.no.covar[,ph,,drop=FALSE]
    
    if(plot.peaks){				
      pdf(pdf.filename, width = nrow(results.no.covar)*0.5, height = 15)
      layout.mat <- get.layout.mat(ncol(pheno.results), "upright")
      # quartz(width = 15, height = 15)
      layout(layout.mat)
      par(mar = c(2,2,2,2))
    }
    
    #bin markers by chromosome
    for(ch in 1:length(chr)){
      if(verbose){
        report.progress(ch, length(chr))
      }
      chr.locale <- which(result.chr == chr[ch])
      chr.results <- pheno.results[chr.locale,,,drop=FALSE]
      
      #bin each allele curve by chromosome, so we don't have bins overlapping 
      #chromosome breaks. Make sure the first bin on one chromosome is 1+ the max
      #bin on the last chromosome
      chr.bins <- apply(chr.results, 3, function(x) bin.curve(x, plot.peaks = plot.peaks, window.size = window.size)$bins)
      if(!is.null(allele.bins[[ph]])){ 
        max.bin <- apply(allele.bins[[ph]], 2, max)
      }else{
        max.bin <- rep(0, ncol(chr.bins))#if we don't have any bins yet, initialize at 0
      }
      total.bins <- Reduce("cbind", lapply(1:ncol(chr.bins), function(x) chr.bins[,x]+max.bin[x]))
      total.bins <- matrix(total.bins, ncol = ncol(chr.bins))
      allele.bins[[ph]] <- rbind(allele.bins[[ph]], total.bins)
    }
    
  }
  #===============================================================
  
  
  #===============================================================			
  #find an effect size cutoff that gives us the approximate number 
  #of alleles requested
  #===============================================================
  if(verbose){cat("\nFinding effect size threshold...\n")}
  
  total.alleles <- 0
  alleles.checked <- NULL
  repeats <- 0 #This checks for bouncing around the same numbers over
  #and over. If we start to see repeat numbers without 
  #getting close to the desired, exit the loop anyway.
  #This prevents infinite loops when the tolerance is 
  #set too small.
  while((total.alleles < (num.alleles - tolerance) || total.alleles > (num.alleles + tolerance)) && repeats == 0){
    ph.alleles <- vector(mode = "list", length = num.pheno)
    
    #adjust the guess point based on how much we need to shift
    #if we have many fewer than we need lower the min.effect.size
    #a lot, if we need only a few more, don't lower it too much
    guess.point <- round(guess.point + num.alleles - total.alleles)
    min.effect.size = min(tail(sorted.results, guess.point))
    
    for(ph in 1:num.pheno){
      pheno.results <- results.no.covar[,ph,,drop=FALSE]
      ph.alleles[[ph]] <- markers.per.peak(pheno.results, allele.bins[[ph]], min.effect.size)
      total.alleles <- round(sum(unlist(ph.alleles))*peak.density)
    }
    if(verbose){
      cat(signif(min.effect.size,3), "\t\t", total.alleles, "\n")
    }
    
    already.seen <- which(alleles.checked == total.alleles)
    if(length(already.seen) > 0){
      repeats = 1
    }
    
    alleles.checked <- c(alleles.checked, total.alleles)
  }
  
  
  if(plot.peaks){
    for(ph in 1:dim(results.no.covar)[[2]]){
      # quartz(width = 11, height = 7)
      pheno.results <- results.no.covar[,ph,,drop=FALSE]
      throw.out <- apply(pheno.results, 2, function(x) bin.curve(x, plot.peaks = plot.peaks, window.size = window.size)$bins)
    }
  }
  
  
  #for each bin that has only 1 allele, flip a biased coin
  #to decide whether to take it.
  ph.alleles <- lapply(ph.alleles, function(x) allele.coin.flip(x, peak.density))
  
  #sample the peaks that reach above the min.effect.size
  #multiply the number of markers from each bin by the peak density
  #to determine how many markers from each bin we should take
  alleles.per.peak <- lapply(ph.alleles, function(x) round(x*peak.density))
  sampled.markers <- vector(mode = "list", length = num.pheno)
  for(ph in 1:num.pheno){	
    #phenotype results across all alleles
    #should be a matrix in two dimensions with 
    #individuals in rows and alleles in columns
    pheno.results <- results.no.covar[,ph,,drop=FALSE]
    pheno.results <- adrop(pheno.results, drop = 2)
    sampled.markers[[ph]] <- sample.peaks(pheno.results = pheno.results, num.per.peak = alleles.per.peak[[ph]], bins = allele.bins[[ph]])			
  }
  
  #pull out the sampled alleled for all phenotypes and build a genotype matrix
  num.parents <- ncol(pheno.results)
  markers.by.parent <- vector(mode = "list", length = num.parents)
  names(markers.by.parent) <- colnames(pheno.results)
  for(p in 1:num.parents){
    markers.by.parent[[p]] <- unique(unlist(lapply(sampled.markers, function(x) x[[colnames(pheno.results)[p]]])))
  }
  
  total.markers <- length(unlist(markers.by.parent))
  if(verbose){cat("total unique alleles:", "\t", total.markers, "\n")}
  
  
  geno.for.pairscan <- matrix(NA, nrow = nrow(data.obj$pheno), ncol = total.markers)
  colnames(geno.for.pairscan) <- 1:ncol(geno.for.pairscan)
  start.allele <- 1
  for(p in 1:length(markers.by.parent)){
    all.allele.markers <- markers.by.parent[[p]]
    if(length(all.allele.markers) > 0){
      geno.section <- geno[,colnames(pheno.results)[p],all.allele.markers,drop=FALSE]
      geno.for.pairscan[,start.allele:(start.allele+length(all.allele.markers)-1)] <- geno.section
      colnames(geno.for.pairscan)[start.allele:(start.allele+length(all.allele.markers)-1)] <- paste(all.allele.markers, colnames(pheno.results)[p], sep = "_")
      start.allele <- start.allele+length(all.allele.markers)
    }
  }
  
  #now order the matrix first by marker then by allele
  marker.split <- strsplit(colnames(geno.for.pairscan), "_")
  just.markers <- unlist(lapply(marker.split, function(x) x[1]))
  just.allele <- unlist(lapply(marker.split, function(x) x[2]))
  marker.table <- cbind(get.marker.num(data.obj, just.markers), just.allele)
  sorted.table <- sortByThenBy(marker.table, col.type = c("n", "c"), return.order = TRUE)
  
  for(i in 1:ncol(sorted.table)){
    geno.for.pairscan <- geno.for.pairscan[,sorted.table[,i]]
  }
  
  
  if(verbose){cat("Checking for linear independence...\n")}
  data.obj$geno_for_pairscan <- geno.for.pairscan
  geno.ind <- get.linearly.independent(data.obj)
  
  rownames(geno.ind$independent.markers) <- rownames(data.obj$pheno)
  data.obj$geno_for_pairscan <- geno.ind$independent.markers
  data.obj$effect_size_cutoff <- min.effect.size
  data.obj$peak_density = peak.density
  data.obj$window_size = window.size
  data.obj$tolerance = tolerance
  
  if(verbose){
    cat(length(geno.ind[[2]]), "allele(s) rejected.\n")
    cat("Final alleles selected:", "\t", ncol(geno.ind$independent.markers), "\n")
  }
  data.obj$marker_selection_method = "top.effects"
  if(plot.peaks){
    dev.off()
  }
  
  return(data.obj)
  
  
}