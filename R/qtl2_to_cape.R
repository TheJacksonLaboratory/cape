#' Convert qtl2 object to cape format
#'
#' This function converts a data object constructed by qtl2 using the read_cross()
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
#' @param genoprobs an optional argument for providing previously calculated genoprobs.
#' if this parameter is missing, genoprobs are calculated by qtl_to_cape.
#' @param map The qtl2 map. This can be omitted if the map is included in the cross 
#' object as either pmap or gmap. By default the physical map (pmap) is used. 
#' If it is missing, the genetic map is used. A provided map will be used 
#' preferentially over a map included in the cross object.
#' @param covar Optional matrix of any covariates to be included in the analysis.
#' @param verbose A logical value indicating whether to print progress to the screen.
#' Defaults to TRUE.
#'
#' @return This function returns a list of two elements. The first element is a cape data
#' object. The second element is a cape genotype object.
#'
#' @import qtl2convert
#' @importFrom qtl2 genoprob_to_alleleprob calc_genoprob
#' 
#' @examples 
#' \dontrun{
#' data_obj <- qtl2_to_cape(cross_obj, genoprobs, map, covar, verbose = TRUE)
#' }
#' 
#'
#' @export

qtl2_to_cape <- function(cross, genoprobs = NULL, map = NULL, covar = NULL, verbose = TRUE){

	phenotype_matrix = as.matrix(cross$pheno)
	
	if(is.null(map)){
		map = cross$pmap
		}
	
	if(is.null(map)){
		map <- cross$gmap
	}
	
	crosstype = cross$crosstype
	mpp_types <- c("do", "riself4", "riself8", "riself16", "magic19")
	is_mpp <- as.logical(length(which(mpp_types == crosstype)))
	
	if(is.null(genoprobs)){
		geno <- cross$geno
		geno_ind <- rownames(geno[[1]])

		if(is_mpp){
		    genoprobs<-probs_doqtl_to_qtl2(geno, map = map, pos_column = "pos")
		    genoprobs<-genoprob_to_alleleprob(genoprobs)
		}else{
			genoprobs <- calc_genoprob(cross, map = map)
		}
	}else{
		geno_ind <- rownames(genoprobs[[1]])
	}
		
	if(is.null(covar)){
		if(!is.null(cross$covar)){
			covar = as.matrix(cross$covar)
			num_covar <- matrix(NA, nrow = nrow(covar),  ncol = ncol(covar))
			dimnames(num_covar) <- dimnames(covar)
				covar_ind <- rownames(covar)

				#convert non-numeric covariates to numeric
				for(i in 1:ncol(covar)){
					as_num <- suppressWarnings(as.numeric(covar[,i]))
					if(all(is.na(as_num))){
						if(verbose){cat("Converting", colnames(covar)[i], "to numeric.\n")}
						new_covar <- as.numeric(as.factor(covar[,i])) - 1
						num_covar[,i] <- new_covar
					}
				}
				covar <- num_covar
			}else{
				covar <- NULL
				covar_ind <- rownames(phenotype_matrix)
			}
	}else{
		covar_ind <- rownames(covar)
	}
	
    chr_to_add <- setdiff(names(genoprobs), c("1", "X", "Y", "M"))

    if(verbose){cat("Converting genoprobs to array...\n")}
    n_dim_geno <- length(dim(genoprobs[[1]]))
    if(n_dim_geno == 2){
    	    geno <- abind(genoprobs[[1]], 1-genoprobs[[1]], along = 3)
		 for(i in chr_to_add){
		 	a_geno <- abind(genoprobs[[i]], 1-genoprobs[[i]], along = 3)
	        geno <- abind(geno, a_geno, along = 2)
    	}
		geno <- aperm(geno, c(1,3,2))
    }else{
        geno <- abind(genoprobs[[1]], along = 3)
		for(i in chr_to_add){
	        geno <- abind(geno, genoprobs[[i]], along = 3)
    	}
    }

	if(!is_mpp && dim(geno)[2] == 3){ #convert to bi-allelic probabilities
		to_biallelic <- function(marker_mat){
			allele1 <- marker_mat[,1]
			het <- marker_mat[,2]
			allele2 <- marker_mat[,3]
			allele1 <- allele1 + het/2
			allele2 <- allele2 + het/2
			#test <- cbind(allele1, allele2)
			biallele_mat <- cbind(allele1, allele2)
			return(biallele_mat)
		}
		geno_temp <- apply(geno, 3, to_biallelic)
		geno <- array(geno_temp, dim = c(nrow(phenotype_matrix), 2, dim(geno)[3]))
		rownames(geno) <- geno_ind
	}

    colnames(geno) <- LETTERS[1:ncol(geno)]
	rownames(geno) <- geno_ind

    common_ind <- Reduce("intersect", list(rownames(phenotype_matrix), geno_ind, covar_ind))
    common_pheno_locale <- match(common_ind, rownames(phenotype_matrix))
    common_geno_locale <- match(common_ind, rownames(geno))
    if(!is.null(covar)){
	    common_covar_locale <- match(common_ind, rownames(covar))
	}
    
    pheno <- as.matrix(phenotype_matrix[common_pheno_locale,])
    rownames(pheno) <- common_ind
    geno <- geno[common_geno_locale,,]


    chr_used <- setdiff(names(genoprobs), c("X", "Y", "M"))
    un_map <- unlist(map[chr_used])
 	marker_location <- as.numeric(un_map)
 
    split_marker <- strsplit(names(un_map), "\\.")
    chr <- as.numeric(sapply(split_marker, function(x) x[1]))
    marker_names <- sapply(split_marker, function(x) x[2])
    dimnames(geno)[[3]] <- marker_names

    geno_names <- dimnames(geno)
    names(geno_names) <- c("mouse", "allele" ,"locus")


    data_obj <- list()
    data_obj$pheno <- pheno
    data_obj$geno_names <- geno_names
    data_obj$chromosome <- chr
    data_obj$marker_num <- 1:length(chr)
	data_obj$marker_location <- marker_location
	
	
    if(!is.null(covar)){
	    data_obj$pheno <- cbind(data_obj$pheno, covar[common_covar_locale,])
	}

    result <- list("data_obj" = data_obj, "geno_obj" = geno)    
    
}