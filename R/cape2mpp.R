#This function converts an old cape object to
#capeDO objects. The genotype object should ideally
#be separate from the rest of the data

cape2mpp <- function(data.obj, geno.obj = NULL){
	
	geno.locale <- which(names(data.obj) == "geno")
	if(length(geno.locale) > 0){
		geno <- data.obj$geno
		}else{
		geno <- geno.obj$geno	
		}
	
	geno.array <- array(NA, dim = c(nrow(geno),2,ncol(geno)))
	geno.array[,1,] <- 1-geno
	geno.array[,2,] <- geno	
	
	dimnames(geno.array) <- list(rownames(geno), c("A", "B"), colnames(geno))
	names(dimnames(geno.array)) <- c("mouse", "allele" ,"locus")
	
	geno.obj$geno <- geno.array
	
	data.obj$geno <- NULL
	data.obj$marker.names <- NULL

	geno.names <- list(rownames(data.obj$pheno), c("A", "B"), colnames(geno))
	names(geno.names) <- c("mouse", "allele", "locus")
	data.obj$geno.names <- geno.names


	results <- list("data.obj" = data.obj, "geno.obj" = geno.obj)
	return(results)

}