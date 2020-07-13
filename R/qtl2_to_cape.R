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
#' @genoprobs an optional argument for providing previously calculated genoprobs.
#' if this parameter is missing, genoprobs are calculated by qtl_to_cape.
#' @map The qtl2 map. This can be omitted if the map is included in the cross 
#' object as either pmap or gmap. By default the physical map (pmap) is used. 
#' If it is missing, the genetic map is used. A provided map will be used 
#' preferrentially over a map included in the cross object.
#'
#' @return This function returns a list of two elements. The first element is a cape data
#' object. The second element is a cape genotype object.
#'
#' @export

qtl2_to_cape <- function(cross, genoprobs = NULL, map = NULL){

	phenotype.matrix = as.matrix(cross$pheno)
	
	if(is.null(map)){
		map = cross$pmap
		}
	
	if(is.null(map)){
		map <- cross$gmap
	}
	
	crosstype = cross$crosstype
	mpp.types <- c("do", "riself4", "riself8", "riself16", "magic19")
	is.mpp <- as.logical(length(which(mpp.types == crosstype)))
	
	if(is.null(genoprobs)){
		geno <- cross$geno
		geno.ind <- rownames(geno[[1]])

		if(is.mpp){
		    genoprobs<-qtl2convert::probs_doqtl_to_qtl2(geno, map = map, pos_column = "pos")
		    genoprobs<-qtl2::genoprob_to_alleleprob(genoprobs)
			}else{
			genoprobs <- calc_genoprob(cross, map = map)
			}
		}else{
			geno.ind <- rownames(genoprobs[[1]])
		}
		
	
	if(!is.null(cross$covar)){
		covar = as.matrix(cross$covar)
		num.covar <- matrix(NA, nrow = nrow(covar),  ncol = ncol(covar))
		dimnames(num.covar) <- dimnames(covar)
    		covar.ind <- rownames(covar)

    		#convert non-numeric covariates to numeric
    		for(i in 1:ncol(covar)){
    			as.num <- suppressWarnings(as.numeric(covar[,i]))
    			if(all(is.na(as.num))){
    				cat("Converting", colnames(covar)[i], "to numeric.\n")
    				new.covar <- as.numeric(as.factor(covar[,i])) - 1
    				num.covar[,i] <- new.covar
    				}
    		}
    		covar <- num.covar
    	}else{
    		covar <- NULL
    		covar.ind <- rownames(phenotype.matrix)
    	}

    
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

	if(!is.mpp && dim(geno)[2] == 3){ #convert to bi-allelic probabilities
		to_biallelic <- function(marker.mat){
			allele1 <- marker.mat[,1]
			het <- marker.mat[,2]
			allele2 <- marker.mat[,3]
			allele1 <- allele1 + het/2
			allele2 <- allele2 + het/2
			#test <- cbind(allele1, allele2)
			biallele.mat <- cbind(allele1, allele2)
			return(biallele.mat)
		}
		geno.temp <- apply(geno, 3, to_biallelic)
		geno <- array(geno.temp, dim = c(nrow(phenotype.matrix), 2, dim(geno)[3]))
		rownames(geno) <- geno.ind
	}

    colnames(geno) <- LETTERS[1:ncol(geno)]
	rownames(geno) <- geno.ind

    common.ind <- Reduce("intersect", list(rownames(phenotype.matrix), geno.ind, covar.ind))
    common.pheno.locale <- match(common.ind, rownames(phenotype.matrix))
    common.geno.locale <- match(common.ind, rownames(geno))
    if(!is.null(covar)){
	    common.covar.locale <- match(common.ind, rownames(covar))
	    }
    
    pheno <- as.matrix(phenotype.matrix[common.pheno.locale,])
    rownames(pheno) <- common.ind
    geno <- geno[common.geno.locale,,]


    chr.used <- setdiff(names(genoprobs), c("X", "Y", "M"))
    un_map <- unlist(map[chr.used])
 	marker.location <- as.numeric(un_map)
 
    split.marker <- strsplit(names(un_map), "\\.")
    chr <- as.numeric(sapply(split.marker, function(x) x[1]))
    marker_names <- sapply(split.marker, function(x) x[2])
    dimnames(geno)[[3]] <- marker_names

    geno.names <- dimnames(geno)
    names(geno.names) <- c("mouse", "allele" ,"locus")


    data.obj <- list()
    data.obj$pheno <- pheno
    data.obj$geno_names <- geno.names
    data.obj$chromosome <- chr
    data.obj$marker_num <- 1:length(chr)
	data.obj$marker_location <- marker.location
	
	
    if(!is.null(covar)){
	    data.obj$pheno <- cbind(data.obj$pheno, covar[common.covar.locale,])
	    }

    result <- list("data.obj" = data.obj, "geno.obj" = geno)    
    
}