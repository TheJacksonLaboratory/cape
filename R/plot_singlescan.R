#' Plot results of single-locus scans
#'
#' This function plots the results of \code{\link{singlescan}}
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param singlescan_obj a singlescan object from \code{\link{singlescan}}
#' @param chr a vector of chromosome names to include in the plot. Defaults to all chromosomes.
#' @param traits a vector of trait names to plot. Defaults to all traits.
#' @param alpha the alpha significance level. Lines for significance values will only
#' be plotted if n_perm > 0 when \code{\link{singlescan}} was run. And only alpha values
#' specified in \code{\link{singlescan}} can be plotted.
#' @param standardized If TRUE t statistics are plotted. If FALSE, effect sizes are plotted.
#' @param color_scheme A character value of either "DO/CC" or other indicating the 
#' color scheme of main effects. If "DO/CC" allele effects can be plotted with the
#' DO/CC colors.
#' @param allele_labels A vector of labels for the alleles if different than those
#' stored in the data_object.
#' @param include_covars Whether to include covariates in the plot.
#' @param show_selected If TRUE will indicate which markers were selected for the pairscan.
#' In order for these to be plotted, \code{\link{select_markers_for_pairscan}} must be run first.
#' @param line_type as defined in plot
#' @param lwd line width, default is 1
#' @param pch see the "points()" R function. Default is 16 (a point).
#' @param cex see the "points()" R function. Default is 1.
#' @param covar_label_size default is 0.7
#' 
#' @importFrom graphics lines
#' 
#' @export

plot_singlescan <- function(data_obj, singlescan_obj, chr = NULL, traits = NULL, 
  alpha = c(0.01, 0.05), standardized = TRUE, color_scheme = c("DO/CC","other"), 
  allele_labels = NULL, include_covars = TRUE, show_selected = FALSE, 
  line_type = "l", lwd = 1, pch = 16, cex = 1, covar_label_size = 0.7){

  oldPar <- par(no.readonly = TRUE)
  on.exit(oldPar)

  geno_names <- data_obj$geno_names
  
  if(is.null(chr)){
    chr <- sort(as.numeric(unique(data_obj$chromosome)))
  }
  
  
  calc_alpha <- singlescan_obj$alpha	
  calc_alpha_locale <- which(calc_alpha %in% alpha)
  if(length(calc_alpha_locale) > 0){
    alpha_to_use = calc_alpha[calc_alpha_locale]
    thresh_to_use = unlist(singlescan_obj$alpha_thresh)[calc_alpha_locale]
  }else{
    alpha_to_use = NULL
    thresh_to_use = NULL
  }
  
  covar_info <- get_covar(data_obj)
  covar_names <- covar_info$covar_names
  
  #Get the dimension names to minimize confusion	
  allele_dim <- which(names(geno_names) == "allele")
  
  all_chromosomes <- data_obj$chromosome
  lod_scores <- singlescan_obj$locus_score_scores
  
  if(!standardized){
    results <- singlescan_obj$singlescan_effects
    plot_type_label <- "beta"
  }else{
    results <- singlescan_obj$singlescan_t_stats
    plot_type_label <- "t_stat"
  }
  
  
  if(is.null(traits)){
    traits <- dimnames(results)[[2]]
  }
  results_el <- which(dimnames(results)[[2]] %in% traits)
  
  if(length(results_el) < length(traits)){
    if(length(results_el) > 0){
      not_found <- traits[-results_el]
    }else{
      not_found <- traits
    }
    warning("The following traits could not be found: ", paste(not_found, collapse = ", "))
    stop()
  }
  
  
  #subset the results based on which chromosomes
  #are desired.
  covar_locale <- which(rownames(results) %in% covar_names)
  chr_locale <- c(which(all_chromosomes %in% chr), covar_locale)
  sub_results <- results[chr_locale,,,drop=FALSE]
  # lod_scores <- lod_scores[chr_locale,,drop=FALSE]
  
  if(include_covars){
    covar_locale <- which(rownames(sub_results) %in% covar_names)
    non_covar_locale <- setdiff(1:nrow(sub_results), covar_locale)
    plot_length <- length(non_covar_locale)
    if(length(covar_locale) > 1){
      covar_x <- segment_region((plot_length+1), round(plot_length*1.1), length(covar_locale))
    }else{
      covar_x <- mean(c((plot_length+1), round(plot_length*1.1)))
    }
  }else{
    covar_locale <- NULL	
    covar_x <- NULL
    non_covar_locale <- which(!rownames(sub_results) %in% covar_names)
    sub_results <- sub_results[non_covar_locale,,,drop=FALSE]
    lod_scores <- lod_scores[non_covar_locale,,drop=FALSE]
  }
  
  max_x <- max(c(length(non_covar_locale), covar_x))

  if(length(results) == 0){
    stop("You must run singlescan.R before plotting effects.")
  }
  
  if(length(dim(results)) < 3){
    stop("This function is for plotting effects of multiple alleles.\nYou only have two alleles at each locus.")
  }
  
  #For each phenotype, we want to plot the effects of the
  #presence of each parent allele across the genome in its
  #own color
  
  num_loci <- dim(sub_results)[[1]]
  num_alleles <- dim(sub_results)[[3]]
  ref_allele <- singlescan_obj$ref_allele
  alleles <- geno_names[[allele_dim]]
  ref_allele_locale <- which(alleles == ref_allele)
  
  used_alleles <- alleles[-ref_allele_locale]
  allele_colors <- get_allele_colors(color_scheme, used_alleles)
  
  if(!is.null(allele_labels)){
    if(length(allele_labels) == length(geno_names[[allele_dim]])){
      used_alleles <- allele_labels[-ref_allele_locale]
    }else{
      used_alleles <- allele_labels	
    }
  }
  
  phenos_scanned  <- dimnames(results)[[2]]
  
  if(show_selected){
    ind_markers <- colnames(data_obj$geno_for_pairscan)
    if(is.null(ind_markers)){stop("select_markers_for_pairscan() must be run before showing selected markers")}
    ind_loci <- apply(matrix(ind_markers, ncol = 1), 1, function(x) strsplit(x, "_")[[1]][1]) 
    ind_alleles <- apply(matrix(ind_markers, ncol = 1), 1, function(x) strsplit(x, "_")[[1]][2]) 
    ind_locale <- which(dimnames(sub_results)[[1]] %in% ind_loci)
  }
  
  
  
  if(standardized){
    ylim <- c(min(c(min(abs(sub_results), na.rm = TRUE), thresh_to_use)), max(c(max(abs(sub_results), na.rm = TRUE), thresh_to_use)))
  }else{
    ylim <- c(min(c(min(sub_results, na.rm = TRUE))), max(c(max(sub_results, na.rm = TRUE))))	
  }
  yrange <- ylim[2]-ylim[1]
  
  
  t_layout_mat <- matrix(c(1,2), nrow = 2)
  eff_layout_mat <- matrix(c(1:3), nrow = 3)
  
  
  for(i in results_el){
    # dev.new(width = 15, height = 5)
    if(plot_type_label == "t_stat"){
      layout(t_layout_mat, heights = c(0.85, 0.15))
    }else{
      layout(eff_layout_mat, heights = c(0.45, 0.4, 0.15))	
    }
    
    if(plot_type_label == "t_stat"){
      if(show_selected){par(mar = c(5,5,7,2))}else{par(mar = c(4,5,7,2))}
      plot.new()
      plot.window(xlim = c(1,max_x), ylim = ylim)
    }else{
      if(show_selected){par(mar = c(5,5,7,2))}else{par(mar = c(4,5,5,2))}
      plot.new()
      plot.window(xlim = c(1,max_x), ylim = c(0, max(lod_scores, na.rm = TRUE)))
      points(non_covar_locale, lod_scores[non_covar_locale,i], type = line_type, pch = pch, cex = cex)
      if(length(covar_locale) > 0){
        points(covar_x, lod_scores[covar_locale,i], type = "h")
      }
      abline(h = 0)
      axis(2)
      mtext("F statistic", side = 2, line = 3)
      legend(0, (max(lod_scores)*1.2), legend = used_alleles, col = allele_colors[,3], lty = 1, lwd = 3, xpd = TRUE, horiz = TRUE)
      par(xpd = TRUE)
      text(covar_x, ylim[2]*-0.05, labels = covar_names, srt = 90, adj = 1, cex = covar_label_size)
      par(xpd = FALSE)
      if(show_selected){par(mar = c(5,5,0,2))}else{par(mar = c(4,5,0,2))}
      plot.new()
      # plot.window(xlim = c(1,num_loci), ylim = ylim)
      plot.window(xlim = c(1,max_x), ylim = ylim)
    }
    
    if(length(covar_locale) > 0){
      covar_effects <- as.vector(sub_results[covar_locale,i,1])
      if(plot_type_label == "t_stat"){
        points(covar_x, abs(covar_effects), col = "black", type = "h", lwd = lwd)
        par(xpd = TRUE)
        text(covar_x, ylim[2]*-0.05, labels = covar_names, srt = 90, adj = 1, cex = covar_label_size)
        par(xpd = FALSE)
        
      }else{
        points(covar_x, covar_effects, col = "black", type = "h", lwd = lwd)	
        par(xpd = TRUE)
        text(covar_x, ylim[2]*-0.05, labels = covar_names, srt = 90, adj = 1, cex = covar_label_size)
        par(xpd = FALSE)
      }
    }
    
    for(j in 1:num_alleles){
      #pull out the effects of the presence of
      #allele j on phenotype i
      allele_effects <- as.vector(sub_results[non_covar_locale,i,j])
      if(plot_type_label == "t_stat"){ #plot the absolute value of the t_statistics
        points(non_covar_locale, abs(allele_effects), col = allele_colors[j,3], type = line_type, lwd = lwd, pch = pch, cex = cex)
      }else{
        points(non_covar_locale, allele_effects, col = allele_colors[j,3], type = line_type, lwd = lwd, pch = pch, cex = cex)	
      }
      
      if(plot_type_label == "t_stat"){
        lines(x = c(1,num_loci), y = rep(data_obj$pairscan_thresh, 2), lty = 1, col = "darkgray")
        lines(x = c(1,num_loci), y = rep(data_obj$covar_thresh, 2), lty = 2, col = "darkgray")
        par(xpd = TRUE)
        if(length(alpha_to_use) > 0){
          for(a in 1:length(alpha_to_use)){
            text(x = num_loci*1.02, y = thresh_to_use[a], labels = paste("p =", alpha_to_use[a]), cex = 0.5, adj = 0)
            par(xpd = FALSE)
            abline(h = thresh_to_use[a], lty = a)
          }
        }
        par(xpd = FALSE)
      }
    } #end looping over alleles
    
    abline(h = 0)
    axis(2)
    # axis(1, labels = FALSE)
    mtext(paste("Effect Relative to Allele", ref_allele), side = 2, line = 2.5)
    mtext(phenos_scanned[i], outer = TRUE, line = -3, cex = 2)
    
    if(plot_type_label == "t_stat"){
      legend(x = 0, y = ylim[2]+yrange*0.15, legend = used_alleles, col = allele_colors[,3], lty = 1, lwd = 3, xpd = TRUE, horiz = TRUE)
    }
    
    #put in lines for chromosome boundaries
    for(ch in 1:length(chr)){
      abline(v = max(which(all_chromosomes == chr[ch])) - min(which(all_chromosomes %in% chr)), lty = 3)
    }	
    
    if(show_selected){
      abline(h = data_obj$effect_size_cutoff)
      par(xpd = TRUE)
      allele_cols <- get_allele_colors(color_scheme, used_alleles)
      y_pos <- ylim[1] - ylim[2]*0.01
      for(cl in 1:nrow(allele_cols)){
        allele_locale <- which(ind_alleles == allele_cols[cl,2])
        if(length(allele_locale) > 0){
          points(ind_locale[allele_locale], rep(y_pos, length(allele_locale)), col = allele_cols[cl,3], pch = 16, cex = 0.5)
          y_pos <- y_pos - ylim[2]*0.01
        }
      }
      par(xpd = FALSE)
    }
    
    par(mar = c(0,5,0,2))	
    plot.new()
    plot.window(xlim = c(1,max_x), ylim = c(0,1))		
    #indicate where chromosomes are and rewrite the
    #phenotype for each, so we can see it on really
    #zoomed in plots
    for(ch in 1:length(chr)){
      #find the mean position of where each chromosome starts and stops
      mid_pos <- mean(which(all_chromosomes == chr[ch])) - min(which(all_chromosomes %in% chr))
      text(mid_pos, 0.7, labels = chr[ch], xpd = TRUE, cex = 1.5)
    }
    
    mtext("Chromosome", side = 1, line = -1.6, cex = 1.2)	
    
  } #end looping over phenotypes
  
}