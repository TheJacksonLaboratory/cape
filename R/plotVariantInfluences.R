#' This function plots m12 and m21 using myImagePlot
#' 
#' This function plots the reparameterized influences of variants on each other. 
#' The epistatic interactions from the pairwise scan are reparameteriezed to 
#' the terms \eqn{m_{12}} and \eqn{m_{21}}, where the subscripts indicate the 
#' source and target variants respectively. These terms are interpreted as the 
#' effect that the source variant exerts on the target variant when both are 
#' present. Negative influences represent suppression while positive influences 
#' represent enhancement.
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param p.or.q A threshold indicating the maximum adjusted p value considered 
#' significant. If an fdr method has been used to correct for multiple testing, 
#' this value specifies the maximum q value considered significant.
#' @param min.std.effect 
#' @param plot.all.vals 
#' @param all.markers
#' @param standardize
#' @param color.scheme
#' @param pos.col 
#' @param neg.col 
#' @param not.tested.col 
#' @param show.marker.labels 
#' @param show.chr 
#' @param label.chr 
#' @param show.alleles 
#' @param scale.effects 
#' @param pheno.width 
#' @param covar.width 
#' @param covar.labels 
#' @param phenotype.labels 
#' @param show.not.tested 
#' 
plotVariantInfluences <- function(data.obj, p.or.q = 0.05, min.std.effect = 0, plot.all.vals = FALSE, 
                                  all.markers = FALSE, standardize = TRUE, color.scheme = c("DO/CC", "other"),
                                  pos.col = "brown", neg.col = "blue", not.tested.col = "lightgray", 
                                  show.marker.labels = FALSE, show.chr = TRUE, label.chr = TRUE, show.alleles = TRUE, 
                                  scale.effects = c("log10", "sqrt", "none"), pheno.width = 11, covar.width = 11, 
                                  covar.labels = NULL, phenotype.labels = NULL, show.not.tested = TRUE){
  require(igraph)
  
  if(!show.not.tested){
    not.tested.col = FALSE
  }
  
  geno.names <- data.obj$geno_names
  marker.names <- geno.names[[3]]
  
  if(length(grep("n", scale.effects)) > 0){
    scale.effects <- "none"
  }
  if(length(scale.effects) == 1){
    if(scale.effects != "log10" & scale.effects != "sqrt" & scale.effects != "none"){
      stop("scale.effects must be 'log10', 'sqrt' or 'none.'")
    }
  }
  
  if(covar.width < 1 || !is.numeric(covar.width)){
    stop("covar.width must be a whole positive number.")
  }
  
  if(pheno.width < 1 || !is.numeric(pheno.width)){
    stop("pheno.width must be a whole positive number.")
  }
  
  var.influences <- data.obj$var_to_var_p_val
  
  pheno.inf <- data.obj$max_var_to_pheno_influence
  if(is.null(phenotype.labels)){
    pheno.names <- names(data.obj$max_var_to_pheno_influence)
  }else{
    pheno.names <- phenotype.labels
    if(length(pheno.names) != length(names(data.obj$max_var_to_pheno_influence))){
      stop("I am detecting the wrong number of phenotype labels for the phenotypes present.")
    }
  }
  num.pheno <- length(pheno.names)
  
  if(not.tested.col == TRUE){
    not.tested.col = "lightgray"
  }
  
  if(is.null(var.influences)){
    stop("calc.p() must be run to calculate variant-to-variant influences.")
  }
  
  if(is.null(pheno.inf)){
    stop("direct.influence() must be run to calculate variant-to-trait influences.")
  }
  
  #This function expands the given rows or columns of 
  #a matrix by a given amount
  expand.matrix <- function(mat, row.col.num, row.or.col, expansion.factor){
    if(row.or.col == "row"){
      row.labels <- 1:nrow(mat)
      for(i in 1:length(row.col.num)){
        mat.before <- mat[which(row.labels < row.col.num[i]),,drop=FALSE]
        mat.after <- mat[which(row.labels > row.col.num[i]),,drop=FALSE]
        row.to.expand <- mat[which(row.labels == row.col.num[i]),,drop=FALSE]
        mat.to.add <- matrix(row.to.expand, nrow = expansion.factor, ncol = ncol(mat), byrow = TRUE)
        label.locale <- which(row.labels == row.col.num[i])
        label <- rownames(mat)[label.locale]
        rowname.v <- rep("", expansion.factor)
        rowname.v[round(expansion.factor/2)] <- label
        rownames(mat.to.add) <- rowname.v
        row.labels <- c(row.labels[which(row.labels < row.col.num[i])], rep(row.col.num[i], expansion.factor), row.labels[which(row.labels > row.col.num[i])])
        mat <- rbind(mat.before, mat.to.add, mat.after)
      }
      return(mat)
    }
    
    
    if(row.or.col == "col"){
      col.labels <- 1:ncol(mat)
      for(i in 1:length(row.col.num)){
        mat.before <- mat[,which(col.labels < row.col.num[i]),drop=FALSE]
        mat.after <- mat[,which(col.labels > row.col.num[i]),drop=FALSE]
        col.to.expand <- mat[,which(col.labels == row.col.num[i]),drop=FALSE]
        mat.to.add <- matrix(col.to.expand, ncol = expansion.factor, nrow = nrow(mat), byrow = FALSE)
        
        label.locale <- which(col.labels == row.col.num[i])
        label <- colnames(mat)[label.locale]
        colname.v <- rep("", expansion.factor)
        colname.v[round(expansion.factor/2)] <- label
        colnames(mat.to.add) <- colname.v
        
        col.labels <- c(col.labels[which(col.labels < row.col.num[i])], rep(row.col.num[i], expansion.factor), col.labels[which(col.labels > row.col.num[i])])
        mat <- cbind(mat.before, mat.to.add, mat.after)
      }
      return(mat)
    }
  }
  
  #This function replaces a row or column in a 
  #matrix with another matrix. If adding rows,
  #the two matrices must match in the number 
  #of columns
  replace.row.col <- function(orig.matrix, replace.matrix, row.col.num, row.or.col){
    if(row.or.col == "row"){
      orig.before <- orig.matrix[1:(min(row.col.num)-1),]
      if(max(row.col.num) < nrow(orig.matrix)){
        orig.after <- orig.matrix[(max(row.col.num)+1):nrow(orig.matrix), ]
      }else{
        orig.after <- NULL
      }
      new.mat <- rbind(orig.before, replace.matrix, orig.after)
      return(new.mat)
    }
    if(row.or.col == "col"){
      orig.before <- orig.matrix[,1:(min(row.col.num)-1)]
      if(max(row.col.num) < nrow(orig.matrix)){
        orig.after <- orig.matrix[,(max(row.col.num)+1):nrow(orig.matrix)]
      }else{
        orig.after <- NULL
      }
      new.mat <- cbind(orig.before, replace.matrix, orig.after)
      return(new.mat)
    }
  }
  
  if(all.markers){
    unique.markers <- geno.names[[3]]
  }else{
    unique.markers <- unique(c(as.vector(var.influences[,"Source"]), as.vector(var.influences[,"Target"]), rownames(pheno.inf[[1]])))
  }
  
  just.markers <- sapply(strsplit(unique.markers, "_"), function(x) x[[1]][1])
  unique.marker.locale <- match(just.markers, marker.names)		
  sorted.markers <- unique.markers[order(unique.marker.locale)]
  
  if(show.alleles){	
    alleles <- unique(sapply(strsplit(colnames(data.obj$geno_for_pairscan), "_"), function(x) x[2]))
    allele.colors <- get.allele.colors(color.scheme, alleles)
    allele.names <- allele.colors[,1]
    allele.labels <- allele.colors[,2]
    all.alleles <- unlist(lapply(strsplit(sorted.markers, "_"), function(x) x[2]))
    allele.cols <- allele.colors[match(all.alleles, alleles),3]
  }else{
    allele.cols <- NULL
  }
  
  
  #update the markers based on the covariate width
  covar.info <- get.covar(data.obj)
  covar.names <- covar.info$covar.names
  if(length(covar.names) > 0){
    covar.markers <- covar.names
    covar.locale <- match(covar.markers, sorted.markers)
    sorted.markers[sort(covar.locale)] <- covar.names #make sure the names are in the right order
    new.covar.markers <- sort(rep(covar.markers, covar.width))
    just.markers <- sorted.markers[-covar.locale]
    expanded.markers <- c(just.markers, new.covar.markers)
    #extend allele cols with gray for covariates
    allele.cols <- c(allele.cols, rep("lightgray", length(new.covar.markers)))
  }else{
    expanded.markers <- sorted.markers
    covar.locale <- NULL	
  }
  
  
  
  #get coordinates of the chromosome boundaries
  if(show.chr){
    num.covar <- length(covar.names)
    orig.chromosomes <- get.marker.chr(data.obj, markers =  sapply(strsplit(sorted.markers, "_"), function(x) x[[1]]))
    chromosomes <- orig.chromosomes[which(orig.chromosomes != 0)]
    if(num.covar > 0){
      for(i in 1:length(covar.names)){
        chromosomes <- c(chromosomes, rep(covar.names[i], covar.width))
      }
    }
    u_chr <- unique(chromosomes[which(!is.na(chromosomes))])
    chr.boundaries <- apply(matrix(u_chr, ncol = 1), 1, function(x) max(which(chromosomes == x))) + 0.5
    chr.boundaries <- c(0, chr.boundaries)
    if(label.chr){
      #use only the first two characters of each chromosome
      chr.names <- unique(chromosomes)
      if(!is.null(covar.labels)){
        chr.names[(length(chr.names)-num.covar+1):length(chr.names)] <- covar.labels
      }
      
      # chr.names <- unlist(lapply(strsplit(unique(chromosomes), ""), function(x) paste(x[1:2], collapse = "")))
      # chr.names <- unlist(strsplit(chr.names, "NA"))
    }else{
      chr.names <- NULL
    }
  }else{
    chr.boundaries <- NULL
    chr.names <- NULL
  }
  
  #make the variant influence matrix from significant interactions
  var.sig.col <- which(colnames(var.influences) == "p.adjusted")
  
  if(standardize){
    edge.weights <- as.numeric(var.influences[,5])
  }else{
    edge.weights <- as.numeric(var.influences[,3])	
  }
  pvals <- as.numeric(var.influences[,var.sig.col])
  
  var.influence.net <- graph_from_edgelist(var.influences[,1:2], directed = TRUE)
  var.pval.net <- var.influence.net
  E(var.influence.net)$weight <- edge.weights
  E(var.pval.net)$weight <- pvals
  var.influence.mat <- as.matrix(as_adjacency_matrix(var.influence.net, attr = "weight"))
  var.pval.mat <- as.matrix(as_adjacency_matrix(var.pval.net, attr = "weight"))
  
  #turn the not-tested elements to NA
  var.influence.mat[which(var.influence.mat == 0)] <- NA
  var.pval.mat[which(var.pval.mat == 0)] <- NA
  
  pheno.influence.mat <- matrix(NA, nrow = length(unique.markers), ncol = num.pheno)
  pheno.pval.mat <- matrix(NA, nrow = length(unique.markers), ncol = num.pheno)
  colnames(pheno.influence.mat) <- colnames(pheno.pval.mat) <- pheno.names
  rownames(pheno.influence.mat) <- rownames(pheno.pval.mat) <- sorted.markers
  
  
  #expand the covariate rows and columns to the specified width
  if(length(covar.names) > 0){
    covar.locale <- which(sorted.markers %in% covar.names)
    
    new.var.inf <- expand.matrix(mat = var.influence.mat, row.col.num = covar.locale, row.or.col = "row", expansion.factor = covar.width)
    new.var.inf <- expand.matrix(new.var.inf, covar.locale, "col", covar.width)
    
    new.var.pval <- expand.matrix(var.pval.mat, covar.locale, "row", covar.width)
    new.var.pval <- expand.matrix(new.var.pval, covar.locale, "col", covar.width)		
    
    var.influence.mat <- new.var.inf
    var.pval.mat <- new.var.pval
  }
  
  
  
  #fill the variant-to-phenotype matrix with test statistics 
  #(still with sources in rows and targets in columns)
  #use phenotypes or eigentraits based on user input
  pheno.sig.col <- which(colnames(pheno.inf[[1]]) == "p.adjusted")
  for(i in 1:length(unique.markers)){
    for(j in 1:length(pheno.names)){
      marker.locale <- which(pheno.inf[[j]][,1] == unique.markers[i])
      if(length(marker.locale) > 0){
        if(standardize){	
          pheno.influence.mat[as.character(unique.markers[i]), pheno.names[j]] <- pheno.inf[[j]][marker.locale, "t.stat"]
        }else{
          pheno.influence.mat[as.character(unique.markers[i]), pheno.names[j]] <- pheno.inf[[j]][marker.locale, "coef"]
        }
        pheno.pval.mat[as.character(unique.markers[i]), pheno.names[j]] <- pheno.inf[[j]][marker.locale, pheno.sig.col]
      }else{
        pheno.influence.mat[as.character(unique.markers[i]), pheno.names[j]] <- NA
        pheno.pval.mat[unique.markers[i], pheno.names[j]] <- NA
      }
    }
  }
  
  #expand the phenotype influence matrix to give it more visual weight in the plot
  expanded.pheno.mat <- expand.matrix(pheno.influence.mat, 1:ncol(pheno.influence.mat), "col", pheno.width)
  expanded.pheno.pval.mat <- expand.matrix(pheno.pval.mat, 1:ncol(pheno.pval.mat), "col", pheno.width)
  
  #also expand the regions where the covariates are if there are covariates
  if(length(covar.names) > 0){
    expanded.pheno.mat <- expand.matrix(expanded.pheno.mat, covar.locale, "row", covar.width)
    expanded.pheno.pval.mat <- expand.matrix(expanded.pheno.pval.mat, covar.locale, "row", covar.width)
  }
  
  
  full.inf.mat <- cbind(var.influence.mat, expanded.pheno.mat)
  full.pval.mat <- cbind(var.pval.mat, expanded.pheno.pval.mat) 
  
  full.inf.mat.num <- apply(full.inf.mat, 2, as.numeric)
  rownames(full.inf.mat.num) <- rownames(full.inf.mat)
  colnames(full.inf.mat.num) <- colnames(full.inf.mat)
  
  full.pval.mat.num <- apply(full.pval.mat, 2, as.numeric)
  rownames(full.pval.mat.num) <- rownames(full.pval.mat)
  colnames(full.pval.mat.num) <- colnames(full.pval.mat)
  
  
  #get the coordinates for all pairs not tested
  # not.tested.locale <- which(is.na(rotate.mat(full.inf.mat)), arr.ind = TRUE)
  not.tested.locale <- which(is.na(rotate.mat(full.inf.mat.num)), arr.ind = TRUE)
  
  if(not.tested.col == FALSE || is.na(not.tested.col)){
    not.tested.locale <- NULL
  }
  
  #take out any values that aren't significant according
  #to the user cutoff, and do not have a high enough
  #effect size
  #if we are not plotting all value
  #use an extra color matrix to highlight significant interactions
  #if we are plotting all markers
  
  extra.col.mat <- NULL
  if(!plot.all.vals){
    full.inf.mat.num[which(full.pval.mat.num > p.or.q)] <- NA
    full.inf.mat.num[which(abs(full.inf.mat.num) < min.std.effect)] <- NA
  }else{
    extra.col.mat <- matrix(NA, nrow = nrow(full.pval.mat.num), ncol = ncol(full.pval.mat.num))
    
    min.effect <- which(abs(full.inf.mat.num) > min.std.effect)
    neg.effect <- intersect(min.effect, which(full.inf.mat.num < 0))
    pos.effect <- intersect(min.effect, which(full.inf.mat.num > 0))
    
    sig.neg <- intersect(neg.effect, which(full.pval.mat.num < p.or.q))
    sig.pos <- intersect(pos.effect, which(full.pval.mat.num < p.or.q))
    
    
    extra.col.mat[sig.neg] <- get.color(neg.col, light.dark = "d")[3]
    extra.col.mat[sig.pos] <- get.color(pos.col, light.dark = "d")[3]
  }
  main <- "Variant Influences"
  
  if(scale.effects == "log10"){
    neg.locale <- which(full.inf.mat.num < 0)
    scaled.effects <- log10(abs(full.inf.mat.num))
    scaled.effects[neg.locale] <- scaled.effects[neg.locale]*-1	
    full.inf.mat.num <- scaled.effects
    main <- "log10 Variant Influences"
  }
  if(scale.effects == "sqrt"){
    neg.locale <- which(full.inf.mat.num < 0)
    scaled.effects <- sqrt(abs(full.inf.mat.num))
    scaled.effects[neg.locale] <- scaled.effects[neg.locale]*-1	
    full.inf.mat.num <- scaled.effects
    main <- "Square Root of Variant Influences"
  }
  
  
  if(length(which(na.omit(as.vector(full.inf.mat.num)) != 0)) == 0){
    plot.new()
    plot.window(xlim = c(0,1), ylim = c(0,1))
    text(0.5, 0.5, "No Significant Interactions")
  }else{
    
    myImagePlot(x = full.inf.mat.num, min.x = min(full.inf.mat.num, na.rm = TRUE), max.x = max(full.inf.mat.num, na.rm = TRUE), main = main, xlab = "Target", ylab = "Source", mark.coords = not.tested.locale, mark.col = not.tested.col, show.labels = show.marker.labels, chromosome.coordinates = chr.boundaries, chr.names = chr.names, show.pheno.labels = TRUE, extra.col.mat = extra.col.mat, allele.cols = allele.cols)
    
    
    #add phenotype names
    if(!is.null(not.tested.locale)){
      legend("topright", legend = "not testable", col = not.tested.col, pch = 16)
    }
  }
  
  invisible(full.inf.mat.num)
  
  
}