#' Convert qtl2 object to cape format
#'
#' This function converts a data object contructed by qtl2 using the read_cross()
#' function to cape format. It returns a list in which the first element is the
#' cape data object, and the second element is the cape genotype object.
#'
#' @references Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use
#'   of pleiotropy to model genetic interactions in a population. PLoS genetics,
#'   8(10), e1003010. doi:10.1371/journal.pgen.1003010
#'
#' @references Broman, Karl W., Daniel M. Gatti, Petr Simecek, Nicholas A. Furlotte, 
#' Pjotr Prins, Śaunak Sen, Brian S. Yandell, and Gary A. Churchill. "R/qtl2: software 
#' for mapping quantitative trait loci with high-dimensional data and multiparent populations." 
#' Genetics 211, no. 2 (2019): 495-502.
#'
#' @param cross a cross object created by the R/qtl2 function read_cross()
#'
#' @return This function returns a list of two elements. The first element is a cape data
#' object. The second element is a cape genotype object.
#'
#' @export

qtl2_to_cape <- function(cross){

	phenotype.matrix = as.matrix(cross$pheno)
	genoprobs = cross$geno
	map = cross$pmap
	covar = as.matrix(cross$covar)
    covar.ind <- rownames(covar)

    #only take numeric covariates
    covar.num <- apply(covar, 2, function(x) suppressWarnings(as.numeric(x)))
    covar.which <- which(apply(covar.num, 2, function(x) !all(is.na(x))))

    if(length(covar.which) < ncol(covar)){
        cat("Removing non-numeric covariates\n")
    }
    covar <- covar.num[,covar.which]
    rownames(covar) <- covar.ind
    
    chr.to.add <- setdiff(names(genoprobs), c("1", "X", "Y", "M"))

    cat("Converting genoprobs to array...\n")
    n.dim.geno <- length(dim(genoprobs[[1]]))
    if(n.dim.geno == 2){
    	    geno <- abind(genoprobs[[1]], 1-genoprobs[[1]], along = 3)
		 for(i in chr.to.add){
		 	a.geno <- abind(genoprobs[[i]], 1-genoprobs[[i]], along = 3)
	        geno <- abind(geno, a.geno, along = 2)
    		}
		geno <- aperm(geno, c(1,3,2))
    }else{
        geno <- abind(genoprobs[[1]], along = 3)
		for(i in chr.to.add){
	        geno <- abind(geno, genoprobs[[i]], along = 3)
    		}
    }
    colnames(geno) <- LETTERS[1:ncol(geno)]

    geno.names <- dimnames(geno)
	
    common.ind <- Reduce("intersect", list(rownames(phenotype.matrix), dimnames(geno)[[1]], rownames(covar)))
    common.pheno.locale <- match(common.ind, rownames(phenotype.matrix))
    common.geno.locale <- match(common.ind, rownames(geno))
    common.covar.locale <- match(common.ind, rownames(covar))
    
    pheno <- as.matrix(phenotype.matrix[common.pheno.locale,])
    rownames(pheno) <- common.ind
    geno <- geno[common.geno.locale,,]

    geno.names <- dimnames(geno)
    names(geno.names) <- c("mouse", "allele" ,"locus")
    
    data.obj <- list()
    data.obj$pheno <- pheno
    data.obj$geno_names <- geno.names
    names(data.obj$geno_names) <- c("mouse", "allele", "locus")
    data.obj$marker_num <- 1:dim(geno)[[3]]
    
    chr.used <- setdiff(names(genoprobs), c("X", "Y", "M"))
    un_map <- unlist(map[chr.used])
    data.obj$marker_location <- as.numeric(un_map)
    split.marker <- strsplit(names(un_map), "\\.")
    chr <- sapply(split.marker, function(x) x[1])
    data.obj$chromosome <- as.numeric(chr)
    data.obj$pheno <- cbind(data.obj$pheno, covar[common.covar.locale,])

    result <- list("data.obj" = data.obj, "geno.obj" = geno)    
    
}