#This function adds a pair interacting QTL to a data object
#specifically for testing kinship corrections
#pop.assig is a numeric vector indicating the population each
#individual belongs to
#allele.freq is a matrix with two rows, one for each SNP and
#the same number of columns as populations. It indicates the
#causal allele frequency (0-1) for each SNP in each population
#main.effect.size is a vector of length 2 indicating the main
#effect of each SNP in the interacting pair
#int.effect.size is a number stating the effect size of the 
#interaction
#pheno is a number indicating which phenotype the QTL relates to
#qtl position is a vector of length 2 indicating where the interacting
#qtl should be placed. If NULL, positions are chosen at random.

sim.qtl.interaction <- function(data.obj, geno.obj, pop.assig, allele.freq, main.effect.size = c(1,1), int.effect.size = 2, pheno = 1, qtl.position = NULL)
{
  
  geno <- get.geno(data.obj, geno.obj)
  
  #break the relationship between all phenotypes and genotypes,
  #but retain phenotype correlations
  data.obj$pheno <- data.obj$pheno[sample(nrow(data.obj$pheno)),]
  
  all.sig <- 0
  all.effects <- c(main.effect.size, int.effect.size)
  num.sig <- length(which(all.effects != 0))
  
  while(all.sig == 0){ #while we haven't made significant QTL for main effects and the interaction
    #simulate genotypes for the different populations
    #in the data object
    num.pop <- length(unique(pop.assig))	
    if(is.null(dim(allele.freq))){
      allele.freq <- as.matrix(allele.freq, nrow = 1)
    }
    if(is.null(dim(main.effect.size))){
      main.effect.size <- as.matrix(main.effect.size, ncol = 1)
    }
    if(dim(allele.freq)[2] != num.pop){
      stop("allele.freq must be a vector with length equal to the number of populations")
    }
  }
  
  #simulate alleles based on the causal allele effect
  #size in all populations
  allele.mat <- matrix(NA, nrow = length(pop.assig), ncol = dim(allele.freq)[1])
  for(p in 1:num.pop){ #for each distinct population in the data set
    pop.locale <- which(pop.assig == p)
    for(m in 1:dim(allele.freq)[1]){ #for each marker we are simulating
      allele.mat[pop.locale,m] <- sim.marker(length(pop.locale), allele.freq[m,p])
    }
  }
  
  #Now we alter the phenotype to be affected by the genotype
  new.pheno <- data.obj$pheno[,ph]
  
  #adjust the phenotype based on the main effects and the interaction effect	
  qtl.pheno <- new.pheno + allele.mat %*% main.effect.size + (allele.mat[,1,drop=FALSE]*allele.mat[,2,drop=FALSE]) %*% int.effect.size
  model <- lm(qtl.pheno~allele.mat[,1]*allele.mat[,2])
  
  
  sig <- which(summary(model)$coefficients[2:4,4] <= 0.05)
  if(length(sig) >= num.sig){
    all.sig <- 1
  }			
  
  print(summary(model))
  
  #pick random places to insert the markers
  if(is.null(qtl.position)){
    rnd.pos <- sort(sample(1:dim(geno)[3], dim(allele.mat)[2], replace = FALSE))
  }else{
    rnd.pos <- qtl.position	
  }
  for(i in 1:length(rnd.pos)){
    geno[,,rnd.pos[i]] <- allele.mat[,i]
  }
  data.obj$pheno[,pheno] <- qtl.pheno
  
  updated.objects <- list("data.obj" = data.obj, "geno.obj" = geno)	
  
  return(updated.objects)
  
}