#This function simulates genotypes for individuals
#for a marker with given allele frequency

sim.marker <- function(n.ind, maf){
  
  genotypes <- c(0, 0.5, 1)
  # genotypes <- c(-1, 0, 1)
  
  p <- 1-maf
  q <- maf
  maj.hom <- p^2
  het <- 2*p*q
  min.hom <- q^2
  
  num.maj.hom = round(maj.hom*n.ind)
  num.het <- round(het*n.ind)
  num.min.hom <- round(min.hom*n.ind)
  
  genos <- c(rep(genotypes[1], num.maj.hom), rep(genotypes[2], num.het), rep(genotypes[3], num.min.hom))
  samp.genos <- sample(genos, n.ind, replace = TRUE)
  
  return(samp.genos)
}

