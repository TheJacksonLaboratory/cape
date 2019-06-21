#' Makes effects plots for multi-allelic 1D scans.
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param singlescan.obj a singlescan object
#' @param chr an array of chromosome names to include in the plot
#' @param traits an array of traits to include in the plot
#' @param alpha the alpha significance level
#' @param standardized boolean; singlescan.t.stats if TRUE, else singlescan.effects
#' @param color.scheme "DO/CC" or "other"
#' @param allele.labels an array of character labels
#' @param include.covars boolean
#' @param show.selected boolean; if TRUE, run select.markers.for.pairscan() first
#' @param line.type 1-character string giving the type of plot desired. The following values are possible, for details, see plot: "p" for points, "l" for lines, "b" for both points and lines, "c" for empty points joined by lines, "o" for overplotted points and lines, "s" and "S" for stair steps and "h" for histogram-like vertical lines. Finally, "n" does not produce any points or lines.
#' @param lwd line width, default is 1
#' @param pch see the "points()" R function. Default is 16 (a point).
#' @param cex see the "points()" R function. Default is 1.
#' @param covar.label.size default is 0.7
#' 
#' @export
plotSinglescan <- function(data.obj, singlescan.obj, chr = NULL, traits = NULL, alpha = c(0.01, 0.05),
                           standardized = TRUE, color.scheme = c("DO/CC","other"), allele.labels = NULL,
                           include.covars = TRUE, show.selected = FALSE, line.type = "l", lwd = 1,
                           pch = 16, cex = 1, covar.label.size = 0.7){
  
  geno.names <- data.obj$geno.names
  
  if(is.null(chr)){
    chr <- sort(as.numeric(unique(data.obj$chromosome)))
  }
  
  
  calc.alpha <- singlescan.obj$alpha	
  calc.alpha.locale <- which(calc.alpha %in% alpha)
  if(length(calc.alpha.locale) > 0){
    alpha.to.use = calc.alpha[calc.alpha.locale]
    thresh.to.use = unlist(singlescan.obj$alpha.thresh)[calc.alpha.locale]
  }else{
    alpha.to.use = NULL
    thresh.to.use = NULL
  }
  
  covar.info <- get.covar(data.obj)
  covar.names <- covar.info$covar.names
  
  
  #Get the dimension names to minimize confusion	
  mouse.dim <- which(names(geno.names) == "mouse")
  locus.dim <- which(names(geno.names) == "locus")
  allele.dim <- which(names(geno.names) == "allele")
  
  
  all.chromosomes <- data.obj$chromosome
  lod.scores <- singlescan.obj$locus.score.scores
  
  if(!standardized){
    results <- singlescan.obj$singlescan.effects
    plot.type.label <- "beta"
  }else{
    results <- singlescan.obj$singlescan.t.stats
    plot.type.label <- "t.stat"
  }
  
  
  if(is.null(traits)){
    traits <- dimnames(results)[[2]]
  }
  results.el <- which(dimnames(results)[[2]] %in% traits)
  
  if(length(results.el) < length(traits)){
    if(length(results.el) > 0){
      not.found <- traits[-results.el]
    }else{
      not.found <- traits
    }
    message("I couldn't find the following traits:")
    cat(not.found, sep = "\n")
    return()
  }
  
  
  #subset the results based on which chromosomes
  #are desired.
  covar.locale <- which(rownames(results) %in% covar.names)
  chr.locale <- c(which(all.chromosomes %in% chr), covar.locale)
  sub.results <- results[chr.locale,,,drop=FALSE]
  # lod.scores <- lod.scores[chr.locale,,drop=FALSE]
  
  if(include.covars){
    covar.locale <- which(rownames(sub.results) %in% covar.names)
    non.covar.locale <- setdiff(1:nrow(sub.results), covar.locale)
    plot.length <- length(non.covar.locale)
    if(length(covar.locale) > 1){
      covar.x <- segment.region((plot.length+1), round(plot.length*1.1), length(covar.locale))
    }else{
      covar.x <- mean(c((plot.length+1), round(plot.length*1.1)))
    }
  }else{
    covar.locale <- NULL	
    covar.x <- NULL
    non.covar.locale <- which(!rownames(sub.results) %in% covar.names)
    sub.results <- sub.results[non.covar.locale,,,drop=FALSE]
    lod.scores <- lod.scores[non.covar.locale,,drop=FALSE]
  }
  
  max.x <- max(c(length(non.covar.locale), covar.x))
  
  
  
  if(length(results) == 0){
    stop("You must run singlescan.R before plotting effects.")
  }
  
  if(length(dim(results)) < 3){
    stop("This function is for plotting effects of multiple alleles.\nYou only have two alleles at each locus.")
  }
  
  #For each phenotype, we want to plot the effects of the
  #presence of each parent allele across the genome in its
  #own color
  
  num.loci <- dim(sub.results)[[1]]
  num.alleles <- dim(sub.results)[[3]]
  ref.allele <- singlescan.obj$ref.allele
  alleles <- geno.names[[allele.dim]]
  ref.allele.locale <- which(alleles == ref.allele)
  
  used.alleles <- alleles[-ref.allele.locale]
  allele.colors <- get.allele.colors(color.scheme, used.alleles)
  
  if(!is.null(allele.labels)){
    if(length(allele.labels) == length(geno.names[[allele.dim]])){
      used.alleles <- allele.labels[-ref.allele.locale]
    }else{
      used.alleles <- allele.labels	
    }
  }
  
  phenos.scanned  <- dimnames(results)[[2]]
  
  
  
  if(show.selected){
    ind.markers <- colnames(data.obj$geno.for.pairscan)
    if(is.null(ind.markers)){stop("select.markers.for.pairscan() must be run before showing selected markers")}
    ind.loci <- apply(matrix(ind.markers, ncol = 1), 1, function(x) strsplit(x, "_")[[1]][1]) 
    ind.alleles <- apply(matrix(ind.markers, ncol = 1), 1, function(x) strsplit(x, "_")[[1]][2]) 
    ind.locale <- which(dimnames(sub.results)[[1]] %in% ind.loci)
  }
  
  
  
  if(standardized){
    ylim <- c(min(c(min(abs(sub.results), na.rm = TRUE), thresh.to.use)), max(c(max(abs(sub.results), na.rm = TRUE), thresh.to.use)))
  }else{
    ylim <- c(min(c(min(sub.results, na.rm = TRUE))), max(c(max(sub.results, na.rm = TRUE))))	
  }
  yrange <- ylim[2]-ylim[1]
  
  
  t.layout.mat <- matrix(c(1,2), nrow = 2)
  eff.layout.mat <- matrix(c(1:3), nrow = 3)
  
  
  for(i in results.el){
    # dev.new(width = 15, height = 5)
    if(plot.type.label == "t.stat"){
      layout(t.layout.mat, heights = c(0.85, 0.15))
    }else{
      layout(eff.layout.mat, heights = c(0.45, 0.4, 0.15))	
    }
    
    if(plot.type.label == "t.stat"){
      if(show.selected){par(mar = c(5,5,7,2))}else{par(mar = c(4,5,7,2))}
      plot.new()
      plot.window(xlim = c(1,max.x), ylim = ylim)
    }else{
      if(show.selected){par(mar = c(5,5,7,2))}else{par(mar = c(4,5,5,2))}
      plot.new()
      plot.window(xlim = c(1,max.x), ylim = c(0, max(lod.scores, na.rm = TRUE)))
      points(non.covar.locale, lod.scores[non.covar.locale,i], type = line.type, pch = pch, cex = cex)
      if(length(covar.locale) > 0){
        points(covar.x, lod.scores[covar.locale,i], type = "h")
      }
      abline(h = 0)
      axis(2)
      mtext("F statistic", side = 2, line = 3)
      legend(0, (max(lod.scores)*1.2), legend = used.alleles, col = allele.colors[,3], lty = 1, lwd = 3, xpd = TRUE, horiz = TRUE)
      par(xpd = TRUE)
      text(covar.x, ylim[2]*-0.05, labels = covar.names, srt = 90, adj = 1, cex = covar.label.size)
      par(xpd = FALSE)
      if(show.selected){par(mar = c(5,5,0,2))}else{par(mar = c(4,5,0,2))}
      plot.new()
      # plot.window(xlim = c(1,num.loci), ylim = ylim)
      plot.window(xlim = c(1,max.x), ylim = ylim)
    }
    
    if(length(covar.locale) > 0){
      covar.effects <- as.vector(sub.results[covar.locale,i,1])
      if(plot.type.label == "t.stat"){
        points(covar.x, abs(covar.effects), col = "black", type = "h", lwd = lwd)
        par(xpd = TRUE)
        text(covar.x, ylim[2]*-0.05, labels = covar.names, srt = 90, adj = 1, cex = covar.label.size)
        par(xpd = FALSE)
        
      }else{
        points(covar.x, covar.effects, col = "black", type = "h", lwd = lwd)	
        par(xpd = TRUE)
        text(covar.x, ylim[2]*-0.05, labels = covar.names, srt = 90, adj = 1, cex = covar.label.size)
        par(xpd = FALSE)
      }
    }
    
    for(j in 1:num.alleles){
      #pull out the effects of the presence of
      #allele j on phenotype i
      allele.effects <- as.vector(sub.results[non.covar.locale,i,j])
      if(plot.type.label == "t.stat"){ #plot the absolute value of the t.statistics
        points(non.covar.locale, abs(allele.effects), col = allele.colors[j,3], type = line.type, lwd = lwd, pch = pch, cex = cex)
      }else{
        points(non.covar.locale, allele.effects, col = allele.colors[j,3], type = line.type, lwd = lwd, pch = pch, cex = cex)	
      }
      
      if(plot.type.label == "t.stat"){
        lines(x = c(1,num.loci), y = rep(data.obj$pairscan.thresh, 2), lty = 1, col = "darkgray")
        lines(x = c(1,num.loci), y = rep(data.obj$covar.thresh, 2), lty = 2, col = "darkgray")
        # abline(h = data.obj$covar.thresh, lty = 2, col = "darkgray")
        par(xpd = TRUE)
        if(length(alpha.to.use) > 0){
          for(a in 1:length(alpha.to.use)){
            text(x = num.loci*1.02, y = thresh.to.use[a], labels = paste("p =", alpha.to.use[a]), cex = 0.5, adj = 0)
            par(xpd = FALSE)
            abline(h = thresh.to.use[a], lty = a)
          }
        }
        par(xpd = FALSE)
      }
    } #end looping over alleles
    
    abline(h = 0)
    axis(2)
    # axis(1, labels = FALSE)
    mtext(paste("Effect Relative to Allele", ref.allele), side = 2, line = 2.5)
    mtext(phenos.scanned[i], outer = TRUE, line = -3, cex = 2)
    
    if(plot.type.label == "t.stat"){
      # TODO find out why the legend entries are missing
      browser()  # 'legend' is of length 0 
      legend(x = 0, y = ylim[2]+yrange*0.15, legend = used.alleles, col = allele.colors[,3], lty = 1, lwd = 3, xpd = TRUE, horiz = TRUE)
    }
    
    #put in lines for chromosome boundaries
    for(ch in 1:length(chr)){
      abline(v = max(which(all.chromosomes == chr[ch])) - min(which(all.chromosomes %in% chr)), lty = 3)
    }	
    
    if(show.selected){
      abline(h = data.obj$effect.size.cutoff)
      par(xpd = TRUE)
      allele.cols <- get.allele.colors(color.scheme, used.alleles)
      y.pos <- ylim[1] - ylim[2]*0.01
      for(cl in 1:nrow(allele.cols)){
        allele.locale <- which(ind.alleles == allele.cols[cl,2])
        if(length(allele.locale) > 0){
          points(ind.locale[allele.locale], rep(y.pos, length(allele.locale)), col = allele.cols[cl,3], pch = 16, cex = 0.5)
          y.pos <- y.pos - ylim[2]*0.01
        }
      }
      par(xpd = FALSE)
    }
    
    par(mar = c(0,5,0,2))	
    plot.new()
    plot.window(xlim = c(1,max.x), ylim = c(0,1))		
    #indicate where chromosomes are and rewrite the
    #phenotype for each, so we can see it on really
    #zoomed in plots
    for(ch in 1:length(chr)){
      #find the mean position of where each chromosome starts and stops
      mid.pos <- mean(which(all.chromosomes == chr[ch])) - min(which(all.chromosomes %in% chr))
      text(mid.pos, 0.7, labels = chr[ch], xpd = TRUE, cex = 1.5)
    }
    
    mtext("Chromosome", side = 1, line = -1.6, cex = 1.2)	
    
  } #end looping over phenotypes
  
  
  
}