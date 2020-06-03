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
#' @param genoprobs a genoprobs object created by the R/qtl2 function calc_genoprobs()
#' @param map the genetic or physical map from an R/qtl2 object.
#'
#' @return This function returns a list of two elements. The first element is a cape data
#' object. The second element is a cape genotype object.
#'
#' @export

qtl2_to_cape <- function(cross, genoprobs, map){

    if(class(cross)[1] == "cross2"){
        along = 3
    }else{
        along = 2
    }

    chr.to.add <- setdiff(names(genoprobs), c("1", "X", "Y", "M"))
    geno <- genoprobs[[1]]
    for(i in chr.to.add){
        geno <- abind(geno, genoprobs[[i]], along = along)
    }

    if(class(cross)[1] != "cross2"){
        #switch the second and third dimensions
        geno <- aperm(geno, c(1,3,2))
        rownames(geno) <- rownames(cross$pheno)
    }
    
    geno.names <- dimnames(geno)
	
    common.ind <- intersect(rownames(cross$pheno), dimnames(geno)[[1]])
    common.pheno.locale <- match(common.ind, rownames(cross$pheno))
    common.geno.locale <- match(common.ind, rownames(geno))
    
    pheno <- as.matrix(cross$pheno[common.pheno.locale,])
    rownames(pheno) <- common.ind
    geno <- geno[common.geno.locale,,]

    geno.names <- dimnames(geno)
    names(geno.names) <- c("mouse", "allele" ,"locus")
    
    data.obj <- list()
    data.obj$pheno <- pheno
    data.obj$geno_names <- geno.names
    data.obj$marker_num <- 1:dim(geno)[[3]]
    
    chr.used <- setdiff(names(genoprobs), c("X", "Y", "M"))
    un_map <- unlist(map[chr.used])
    data.obj$marker_location <- as.numeric(un_map)
    split.marker <- strsplit(names(un_map), "\\.")
    chr <- sapply(split.marker, function(x) x[1])
    data.obj$chromosome <- as.numeric(chr)
    data.obj$p_covar_table <- as.matrix(cross$covar)
    rownames(data.obj$p_covar_table) <- common.ind
    data.obj$p_covar <- colnames(cross$covar)

    result <- list("data.obj" = data.obj, "geno.obj" = geno)    
    
}