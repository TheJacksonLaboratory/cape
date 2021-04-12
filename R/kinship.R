#' Calculate the kinship matrix
#' 
#' This function produces a realized relationship matrix 
#' (kinship matrix) for use in adjusting for the effect 
#' of inbred relatedness. We use the R/qtl2 function
#' calc_kinship.
#' 
#' Broman KW, Gatti DM, Simecek P, Furlotte NA, Prins P, 
#' Sen Ś, Yandell BS, Churchill GA (2018) R/qtl2: software
#' for mapping quantitative trait loci with high-dimensional 
#' data and multi-parent populations. Genetics
#' 211:495-502 doi:10.1534/genetics.118.301595
#' 
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object
#' @param type type of  kinship correction. Default is overall.
#' @param pop population type, "MPP" (multi-parental population), 
#' "2PP" (2 parents), "RIL" (recombinant inbred line)
#' @param n_cores The number of cores. Defaults to 4.
#' 
#' @details This uses the function probs_doqtl_to_qtl2 
#' from qtl2convert:
#' Karl W Broman (2019). qtl2convert: Convert Data 
#' among R/qtl2, R/qtl, and DOQTL. 
#' \url{https://kbroman.org/qtl2/},
#' \url{https://github.com/rqtl/qtl2convert/}.
#' And genoprob_to_alleleprob from qtl2.
#'
#' @return This function returns an n by n matrix, where 
#' n is the number of individuals in the test population. 
#' The entries of the matrix represent the level of relatedness
#' between pairs of individuals. For more information see
#' Kang, H. M. et al. Efficient control of population 
#' structure in model organism association mapping. Genetics 
#' 178, 1709–1723 (2008).
#' 
#' @import qtl
#' @import qtl2convert
#' @importFrom qtl2 calc_kinship genoprob_to_alleleprob
#' @importFrom utils combn
#' 
#' @examples 
#' \dontrun{
#' K <- kinship(data_obj, geno_obj)
#' }
#'
#' @export
kinship <- function(data_obj, geno_obj, type=c("overall"), n_cores=4, 
  pop=c("MPP","2PP","RIL")){
  #file input could be geno_obj or genoprobs
  
  
  ##############################################################
  #                                                            #
  #        Determine if locus is numerical or character        #
  #                                                            #
  ##############################################################
  num <- "1"
  snp <- data_obj$geno_names$locus[[1]]
  
  if(grepl(num,snp)==TRUE) {
    locus<-"num"
  } else {
    locus<-"char"
  }
  
  #################################################################
  #                                                               #
  # Create probability and map file if locus is numerical for MPP #
  #                                                               #
  #################################################################

  ##check to see if genoprobs have been calculated, if not calculate genotype probablities
  if(!("calc_genoprob" %in% class(geno_obj)) & locus=="num" & pop=="MPP"){
    
    ### create map and genoprobs using geno file
    
    map <- data.frame(marker=dimnames(geno_obj)[[3]],chr=dimnames(geno_obj)[[3]],pos=dimnames(geno_obj)[[3]],stringsAsFactors = F)
    
    temp <- strsplit(map$chr,":") #splits column into two separate columns
    
    map$chr <- sapply(temp,"[",1) #pulls first column and makes it a list
    
    map$pos <- as.numeric(sapply(temp,"[",2)) #pulls second column and makes it a list
    
    #data_obj$save_rds(map,"map.RDS")
    
    genoprobs <- probs_doqtl_to_qtl2(geno_obj,map=map,pos_column = "pos") #creates genotype probabilities from DOqtl...can only be used for DO genotype file
    
    genoprobs <- genoprob_to_alleleprob(genoprobs)
  }
  
  
  #################################################################
  #                                                               #
  # Create probability and map file if locus is character for MPP #
  #                                                               #
  #################################################################
  
  if(!("calc_genoprob" %in% class(geno_obj)) & locus=="char" & pop=="MPP"){
    
    ### create map and genoprobs using geno file
    map <- data.frame(marker = data_obj$geno_names[[3]],chr = data_obj$geno_names[[3]], 
    		pos = data_obj$geno_names[[3]], stringsAsFactors = F)
    
    map$marker <- as.list(data_obj$geno_names$locus)
    
    map$chr <- as.list(data_obj$chromosome)
    
    map$pos <- as.list(data_obj$marker_location)
    
    #data_obj$save_rds(map,"map.RDS")
    
    genoprobs <- probs_doqtl_to_qtl2(geno_obj,map=map,pos_column = "pos") #creates genotype probabilities from DOqtl...can only be used for DO genotype file
    
    genoprobs <- genoprob_to_alleleprob(genoprobs)
    
  }
  
  # The QTL format file is a temporary object for transferring population data to R.qtl
  qtl_file <- "QTL_format.csv"
  qtl_path <- data_obj$results_path
  
  ##############################################################
  #                                                            #
  #          Create probability and map file if RIL            #
  #                                                            #
  ##############################################################
  if(!("calc_genoprob" %in% class(geno_obj)) & pop=="RIL"){
    write_population(data_obj, geno_obj, filename = file.path(qtl_path, qtl_file), na = "")
    cross <- read.cross(format="csv", dir = qtl_path, qtl_file, genotypes=c(0,.5,1))
    unlink(file.path(qtl_path, qtl_file)) #delete the file
    cross <- convert2risib(cross)
    cross <- jittermap(cross)
    qtlprobs <- calc_genoprob(cross)
    probs <- probs_qtl_to_qtl2(qtlprobs)
    genoprobs <- probs$probs
    map <- probs$map
  }
  
  ##############################################################
  #                                                            #
  #           Create probability and map file if 2PP           #
  #                                                            #
  ##############################################################
  if(!("calc_genoprob" %in% class(geno_obj)) & pop=="2PP"){
    write_population(data_obj, geno_obj, filename = file.path(qtl_path, qtl_file), 
    na = "")
    cross <- read.cross(format="csv", dir = qtl_path, qtl_file, genotypes=c(0,.5,1))
    unlink(file.path(qtl_path, qtl_file)) #delete the file
    qtlprobs <- calc_genoprob(cross)
    probs <- probs_qtl_to_qtl2(qtlprobs)
    genoprobs <- probs$probs
    map <- probs$map
  }
  
  ##############################################################
  #                                                            #
  #            Create overall kinship matrix if MPP            #
  #                                                            #
  ##############################################################
  if(type=="overall"){
    
    if(pop=="MPP"){
      
      ## calculate kinship matrix using genotype or allele probabilities
      if(type=="chr"){
        stop("Must be type overall")
      }
      else if(type=="loco"){
        stop("Must be type overall")
      }

      map <- map_df_to_list(map,pos_column = "pos")
      
      #data_obj$save_rds(map,"map.RDS")
      
      K <- calc_kinship(probs=genoprobs,type=type, cores=n_cores)
    }
    
    
    ##############################################################
    #                                                            #
    #        Create overall kinship matrix if RIL or 2PP         #
    #                                                            #
    ##############################################################
    
    if(pop== "RIL"|| pop == "2PP"){
      
      if(type=="chr"){
        stop("Must be type overall")
      }
      else if(type=="loco"){
        stop("Must be type overall")
      }
      
      kinship <- calc_kinship(probs=genoprobs,type=type, cores=n_cores)
      rownames(kinship) <- colnames(kinship) <- rownames(data_obj$pheno)
      K <- list(kinship)
      names(K)[1] <- "overall"

    } 
    
    # TODO these calls to the class() function don't work; it just returns the object type
    # TODO the class(file) call just returns "function"
    if ("calc_genoprob" %in% class(file)){
      
      ## Convert to allele probabilities if it isn't a 2PP 
      if(pop=="MPP"){ genoprobs<-genoprob_to_alleleprob(genoprobs)}
      
      if(type=="chr"){
        stop("Must be type overall")
      }
      else if(type=="chr"){
        stop("Must be type overall")
      }
      
      kinship <- calc_kinship(probs=genoprobs,type=type, cores=n_cores)
      rownames(kinship) <- colnames(kinship) <- rownames(data_obj$pheno)
      K <- list(kinship)
      names(K)[1] <- "overall"
      
    }
  }
  
  ##############################################################
  #                                                            #
  #       Create Leave two chromosome out kinship matrix       #
  #       For now this block of code will never run as we      #
  #       have removed the ltco option, because of some        #
  #       instability                                          #
  #                                                            #
  ##############################################################
  
  if(type=="ltco"){
    # create the list of chromosome pairs
    exclude_none <- matrix(c("-0", "-0"), nrow=2, ncol=1, byrow=TRUE)
    exclude_one <- matrix(names(genoprobs), nrow=2, ncol=length(genoprobs), byrow=TRUE)
    exclude_two <- combn(names(genoprobs), 2)
    
    # for all 2-chromosome combinations, filter them out of the geno object
    excludes <- cbind(exclude_none, exclude_one, exclude_two)
    
    # assign names to elements like "1,1", "1,2" ...
    chr_names <- function(x){return(paste(x, collapse=","))}
    colnames(excludes) <- apply(excludes, 2, chr_names)
    
    # a place to hold all the kinship matrices
    K = list()
    
    # get a pair to exclude, e.g.,  excludes[,"18,19"] gives c("18, "19")
    for(i in colnames(excludes)) {
      
      # calculate kinship matrix with some chromosomes excluded
      gp <- genoprobs[,!names(genoprobs) %in% excludes[,i]]
      kinship <- calc_kinship(probs=gp,type="overall", cores=n_cores)
      
      # assign row,column names using names from the pheno object
      rownames(kinship) <- colnames(kinship) <- rownames(data_obj$pheno)
      
      K[[i]] <- kinship
    }
    
  }
  # rename the "-0,-0" element to "overall"
  names(K)[1] <- "overall"
  if(length(K) == 1){
    K <- K[[1]]
  }

  return(K)
}
