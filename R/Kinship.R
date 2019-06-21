#' Calculate the kinship correction matrix
#' 
#' This function produces a realized relationship matrix (kinship matrix)
#' for use in adjusting for the effect of inbred relatedness 
#' plot.adj.mat is used to plot the covariance matrices
#' vs. the positive definite matrices
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object
#' @param type correction type, must be either "overall" or "ltco" (leave two chromosomes out)
#' @param pop population type, "MPP" (multi-parental population), "2PP" (2 parents), "RIL" (recombinant inbred line)
#' @param n.cores integer, default 4
#'
#' @export
Kinship<-function(data.obj, geno.obj, type=c("overall","ltco"), n.cores=4, pop=c("MPP","2PP","RIL")){
  #file input could be geno.obj or genoprobs
  
  
  ##############################################################
  #                                                            #
  #        Determine if locus is numerical or character        #
  #                                                            #
  ##############################################################
  num<-"1"
  snp<-data.obj$geno_names$locus[[1]]
  
  if(grepl(num,snp)==TRUE)
  {locus<-"num"
  }else {locus<-"char"}
  
  #################################################################
  #                                                               #
  # Create probability and map file if locus is numerical for MPP #
  #                                                               #
  #################################################################

  
  ##check to see if genoprobs have been calculated, if not calculate genotype probablities
  if(!("calc_genoprob" %in% class(geno.obj)) & locus=="num" & pop=="MPP"){
    
    ### create map and genoprobs using geno file
    
    map<-data.frame(marker=dimnames(geno.obj)[[3]],chr=dimnames(geno.obj)[[3]],pos=dimnames(geno.obj)[[3]],stringsAsFactors = F)
    
    temp<-strsplit(map$chr,":") #splits column into two separate columns
    
    map$chr<-sapply(temp,"[",1) #pulls first column and makes it a list
    
    map$pos<-as.numeric(sapply(temp,"[",2)) #pulls second column and makes it a list
    
    saveRDS(map,"map.RData")
    
    genoprobs<-probs_doqtl_to_qtl2(geno.obj,map=map,pos_column = "pos") #creates genotype probabilities from DOqtl...can only be used for DO genotype file
    
    genoprobs<-genoprob_to_alleleprob(genoprobs)
  }
  
  
  #################################################################
  #                                                               #
  # Create probabiltiy and map file if locus is character for MPP #
  #                                                               #
  #################################################################
  
  if(!("calc_genoprob" %in% class(geno.obj)) & locus=="char" & pop=="MPP"){
    
    ### create map and genoprobs using geno file
    map<-data.frame(marker=dimnames(geno.obj)[[3]],chr=dimnames(geno.obj)[[3]],pos=dimnames(geno.obj)[[3]],stringsAsFactors = F)
    
    map$marker<-as.list(data.obj$geno_names$locus)
    
    map$chr<-as.list(data.obj$chromosome)
    
    map$pos<-as.list(data.obj$marker_location)
    
    saveRDS(map,"map.RData")
    
    genoprobs<-probs_doqtl_to_qtl2(geno.obj,map=map,pos_column = "pos") #creates genotype probabilities from DOqtl...can only be used for DO genotype file
    
    genoprobs<-genoprob_to_alleleprob(genoprobs)
    
  }
  
  ##############################################################
  #                                                            #
  #          Create probability and map file if RIL            #
  #                                                            #
  ##############################################################
  if(!("calc_genoprob" %in% class(geno.obj)) & pop=="RIL"){
    writePopulation(data.obj,geno.obj,filename = "QTL_format.csv",na = "")
    cross<-read.cross(format="csv",".","QTL_format.csv", genotypes=c(0,.5,1))
    unlink("QTL_format.csv") #delete the file
    cross<-convert2risib(cross)
    cross<-jittermap(cross)
    qtlprobs<-calc.genoprob(cross)
    probs<-probs_qtl_to_qtl2(qtlprobs)
    genoprobs<-probs$probs
    map<-probs$map
  }
  
  ##############################################################
  #                                                            #
  #           Create probability and map file if 2PP           #
  #                                                            #
  ##############################################################
  if(!("calc_genoprob" %in% class(geno.obj)) & pop=="2PP"){
    writePopulation(data.obj,geno.obj,filename = "QTL_format.csv",na = "")
    cross<-read.cross(format="csv",".","QTL_format.csv", genotypes=c(0,.5,1))
    unlink("QTL_format.csv") #delete the file
    qtlprobs<-calc.genoprob(cross)
    probs<-probs_qtl_to_qtl2(qtlprobs)
    genoprobs<-probs$probs
    map<-probs$map
  }
  
  ##############################################################
  #                                                            #
  #            Create overall Kinship matrix if MPP            #
  #                                                            #
  ##############################################################
  if(type=="overall"){
    
    if(pop=="MPP"){
      
      ## calculate kinship matrix using genotype or allele probabilities
      if(type=="chr"){
        stop("Must be type overall or ltco")
      }
      else if(type=="loco"){
        stop("Must be type overall or ltco")}
      
      map<-map_df_to_list(map,pos_column = "pos")
      
      saveRDS(map,"map.RData")
      
      K<-calc_kinship(probs=genoprobs,type=type, cores=n.cores)}
    
    
    ##############################################################
    #                                                            #
    #        Create overall Kinship matrix if RIL or 2PP         #
    #                                                            #
    ##############################################################
    
    if(pop== "RIL"|| pop == "2PP"){
      
      if(type=="chr"){
        stop("Must be type overall or ltco")
      }
      else if(type=="loco"){
        stop("Must be type overall or ltco")}
      
      kinship <- calc_kinship(probs=genoprobs,type=type, cores=n.cores)
      rownames(kinship) <- colnames(kinship) <- rownames(data.obj$pheno)
      K <- list(kinship)
      names(K)[1] <- "overall"
      
    } 
    
    # TODO these calls to the class() function don't work; it just returns the object type
    # TODO the class(file) call just returns "function"
    if ("calc_genoprob" %in% class(file)){
      
      ## Convert to allele probabilities if it isn't a 2PP 
      if(pop=="MPP"){ genoprobs<-genoprob_to_alleleprob(genoprobs)}
      
      if(type=="chr"){
        stop("Must be type overall, or ltco")
      }
      else if(type=="chr"){
        stop("Must be type overall or ltco")}
      
      kinship <- calc_kinship(probs=genoprobs,type=type, cores=n.cores)
      rownames(kinship) <- colnames(kinship) <- rownames(data.obj$pheno)
      K <- list(kinship)
      names(K)[1] <- "overall"
      
    }
  }
  
  ##############################################################
  #                                                            #
  #       Create Leave two chromosome out Kinship matrix       #
  #                                                            #
  ##############################################################
  
  if(type=="ltco"){
    # create the list of chromosome pairs
    exclude.none <- matrix(c("-0", "-0"), nrow=2, ncol=1, byrow=TRUE)
    exclude.one <- matrix(names(genoprobs), nrow=2, ncol=length(genoprobs), byrow=TRUE)
    exclude.two <- combn(names(genoprobs), 2)
    
    # for all 2-chromosome combinations, filter them out of the geno object
    excludes <- cbind(exclude.none, exclude.one, exclude.two)
    
    # assign names to elements like "1,1", "1,2" ...
    chr.names <- function(x){return(paste(x, collapse=","))}
    colnames(excludes) <- apply(excludes, 2, chr.names)
    
    # a place to hold all the kinship matrices
    K = list()
    
    # get a pair to exclude, e.g.,  excludes[,"18,19"] gives c("18, "19")
    for(i in colnames(excludes)) {
      
      # calculate kinship matrix with some chromosomes excluded
      gp <- genoprobs[,!names(genoprobs) %in% excludes[,i]]
      kinship <- calc_kinship(probs=gp,type="overall", cores=n.cores)
      
      # assign row,column names using names from the pheno object
      rownames(kinship) <- colnames(kinship) <- rownames(data.obj$pheno)
      
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
