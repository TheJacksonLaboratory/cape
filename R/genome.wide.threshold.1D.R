#This function does a permutation test on the 1D scan
#to calculate the genome-wide significance threshold
#of effect.
#A threshold function must be specified for both the 
#selection of alleles to use in the pair scan and for
#the alleles to be used as covariates. Both thresholds
#are added to the data object along with their function
#call. To fit the p values to the extreme value distribution,
#the previous default, use p.val.thresh() with your desired
#alpha as the threshold parameter.
#an additional option is the t.stat.thresh.sd in which the 
#threshold is calculated based on standard deviations away
#from the mean of the t statistic.
#additional thresholding functions that use the permuted 
#t statistics and a single parameter can easily be added
#to the arsenal.

#' Calculate a genome-wide significance threshold for the single-variant scan
#' 
#' This function uses permtuation testing to calculate a genome-wide significance
#' threshold for the single-variant scan. Two user-defined thresholds are 
#' calculated: the lower threshold (higher alpha) can be used to determine which 
#' variants will be used in the pairwise scan. The higher threshold (lower alpha) 
#' is used to determine which variants are used as covariates in the pairwise scan. 
#' In each permutation, the phenotype or eigentrait is shuffled and all markers 
#' are retested with the permuted phenotype. The regression coefficients are 
#' collected from each permutation and the extreme value distribution is used to 
#' determine thresholds for the user-defined alpha values.
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object. If this is not supplied then it is 
#' generated here.
#' @param n.perm The number of permutations to perform. The default is 100.
#' @param scan.what A character value that uniquely specifies whether the 
#' eigentraits or phenotypes should be scanned. Options are "eigentraits"
#' or "raw.traits".
#' @param ref.allele \code{\link{Cape}} requires that one of the alleles 
#' in the population be selected as a reference allele. The effects of all 
#' alleles are then reported as the effects relative to the reference allele. 
#' In the DO population B represents the Black6 mouse. Because this strain is 
#' considered a very standard strain, the default reference allele is B. 
#' If another allele has a positive effect, it means that it increases the given 
#' phenotype relative to that in the Black6 mouse. Any other allele can be selected 
#' as the reference allele simply by specifying it with this argument. 
#' @param alpha The alpha value used to calculate the threshold for which variants 
#' will be used as covariates in the pairwise scan.
#' @param model.family Indicates the model family of the phenotypes. This can be 
#' either "gaussian" or "binomial".
#' @param run.parallel boolean, default = TRUE
#' @param n.cores integer number of cores to use if running in parallel
#' 
genome.wide.threshold.1D <- function(data.obj, geno.obj = NULL, n.perm = 100, 
                                     scan.what = c("eigentraits", "raw.traits"), 
                                     ref.allele = NULL, alpha = c(0.01, 0.05), 
                                     model.family, run.parallel = TRUE, n.cores = 4
                                     verbose = verbose){
  
  if(!run.parallel){n.cores = 1}
  
  require("evd")
  
  if(n.perm < 2){
    stop("You must do more than one permutation.")
  }
  
  gene <- get.geno(data.obj, geno.obj)
  
  #Get the dimension names to minimize confusion
  locus.dim <- which(names(dimnames(gene)) == "locus")
  allele.dim <- which(names(dimnames(gene)) == "allele")
  
  #check to make sure a reference allele has been chosen
  #and that we can find it in the array
  if(length(ref.allele) != 1){
    stop("You must pick one reference allele")
  }else{
    ref.col <- which(dimnames(gene)[[allele.dim]] == ref.allele)
    if(length(ref.col) < 1){
      stop("I couldn't find reference allele: ", ref.allele)
    }
  }
  
  
  #calculate the numbers of markers, phenotypes and samples
  # TODO remove the lines below if not needed
  # n.gene <- dim(gene)[locus.dim]
  # n.allele <- dim(gene)[allele.dim]-1 #find the number of alleles we will be regressing on
  
  #pull out genotype and phenotype data for
  #clarity of code later.
  #If the user has not specified a scan.what,
  #from the outer function (singlescan.R),
  #default to eigentraits, basically, if eigen,
  #et, or ET are anywhere in the string, use the
  #eigentrais, otherwise, use raw phenotypes
  use.eigentraits <- length(c(grep("eigen", scan.what), grep("ET", scan.what), grep("et", scan.what)))
  if(use.eigentraits){
    pheno <- data.obj$ET
  }else{
    pheno <- data.obj$pheno
  }
  
  
  num.samples <- dim(pheno)[1]
  
  if(length(model.family) == 1){
    model.family <- rep(model.family, ncol(pheno))
  }
  
  
  #============================================================	
  #functions
  #============================================================	
  #do the regression at a single locus and return the model
  multi.allele.regress <- function(locus.mat, phenotype, model.family){
    #remove the reference allele from the matrix
    regress.mat <- locus.mat[,-ref.col]
    model <- glm(phenotype~regress.mat, family = model.family)
    model.coef <- coef(summary(model))
    #gather the t statistics for all alleles
    t.stats <- model.coef[2:dim(model.coef)[1], 3]
    return(as.vector(t.stats))
  }
  
  get.s <- function(evd.result, alpha){
    s <- qgev(1-alpha,loc=evd.result$estimate[1], scale=evd.result$estimate[2], shape=evd.result$estimate[3], lower.tail = TRUE)
    return(s)
  }
  
  
  #This funtion runs a single permutation
  #it can be parallelized to speed things up
  one.perm <- function(){
    #shuffle the vector of individuals
    #This keeps all correlation structure
    #between phenotypes and genotypes the
    #same, but reassigns which phenotypes
    #are associated with which genotypes
    sampled.vector <- sample(1:num.samples)
    
    #permute the individuals
    gene_perm <- gene[sampled.vector,,]
    
    #for each phenotype, do the regression on all loci
    #at once and pull out the t statistic
    
    stat.mat <- array(NA, dim = c(dim(gene_perm)[locus.dim], ncol(gene_perm)-1, ncol(pheno)))
    
    # ptm <- proc.time()
    for(et in 1:length(pheno[1,])){	 	  			
      stat.mat[,,et] <- Reduce("rbind", lapply(1:dim(gene_perm)[locus.dim], function(x) multi.allele.regress(gene_perm[,,x], phenotype = pheno[,et], model.family = model.family[et])))
    }
    
    #find the maximum t statistic for the permutation
    #and record it.
    max.stat <- apply(stat.mat, 2, max)
    return(max.stat)        				
  }
  
  #====================================================
  
  # TODO remove the two lines below if not needed
  #create a matrix to hold permutation results
  #perm.max <- array(NA, dim = c(n.perm, ncol(gene)-1, ncol(pheno)))
  
  
  if (run.parallel) {
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    max.stat <- foreach::foreach(p = 1:n.perm, .combine = "rbind")  %dopar% {
      one.perm()
    }
    parallel::stopCluster(cl)
    
  } else {
    
    max.stat <- c()
    index <- 1:n.perm
    for (p in index) {
      if(verbose){report.progress(p, n.perm)}
      max.stat <- rbind(max.stat, one.perm())
    }
    
  }
  # proc.time() - ptm
  
  #====================================================
  # calculate the pairscan and covar thresholds and 
  # add them to the data object.
  #====================================================
  
  #apply the extreme value distribution to the results
  evd <- apply(max.stat, 2, function(x) fgev(x, std.err = FALSE))
  data.obj$save_rds(max.stat, "singlescan.permutations.RData")
  
  s <- vector(mode = "list", length = length(alpha))
  for(a in 1:length(alpha)){
    s[[a]] <- as.vector(sapply(evd, function(x) get.s(x, alpha[a])))
  }
  
  #calculate one threshold over all phenotypes
  thresholds <- lapply(s, mean)
  names(thresholds) <- alpha
  
  if(verbose(cat("\n")))

  return(thresholds)
  
  
}
