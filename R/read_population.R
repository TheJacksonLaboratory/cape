#' Reads in data in the R/qtl csv format
#' 
#' This function reads in a data file in the r/qtl format
#' It converts letter genotypes to numbers if required.
#' It parses the data into a data object.
#' if filename is left empty, the script will ask the
#' use to choose a file.
#' phenotypes can be specified with a vector of column 
#' numbers or character strings. For each phenotype
#' specified with a name, the script will find its location. 
#' 
#' @param filename The name of the file to read in
#' @param pheno_col Column numbers of desired traits. The default
#' behavior is to read in all traits.
#' @param geno_col Column numbers of desired markers. The default
#' behavior is to read in all markers.
#' @param id_col The column number of an ID column. This is helpful to
#' specify if the individual IDs are strings. Strings are only
#' allowed in the ID column. All other trait data must be numeric.
#' @param delim column delimiter for the file, default is ","
#' @param na_strings a character string indicating how NA values are specified, default is "-"
#' @param check_chr_order boolean, default is TRUE
#' @param verbose A  logical value indicating whether to print progress 
#' and cross information to the screen. Defaults to TRUE.
#' 
#' @references Broman et al. (2003) R/qtl: QTL mapping in experimental crosses. 
#' Bioinformatics 19:889-890 doi:10.1093/bioinformatics/btg112
#'
#' @return This function returns a cape object in a former cape format.
#' It must be updated using \code{\link{cape2mpp}}
#' 
#' @examples 
#' \dontrun{
#' cape_obj <- read_population("cross.csv")
#' combined_obj <- cape2mpp(cape_obj)
#' data_obj <- combined_obj$data_obj
#' geno_obj <- combined_obj$geno_obj
#' }
#'
#' @export
read_population <- function(filename = NULL, pheno_col = NULL, geno_col = NULL, id_col = NULL, 
	delim = ",", na_strings = "-", check_chr_order = TRUE, verbose = TRUE) {

		if(is.null(filename)){
			filename <- file.choose()
		}
			
		cross_data <- read.table(filename, na.strings = na_strings, stringsAsFactors = FALSE, 
		sep = delim, header = TRUE)

		if(!is.null(id_col)){
			ind_names <- cross_data[3:nrow(cross_data),id_col]
			cross_data <- cross_data[,-id_col]
		}else{
			ind_names <- 1:(nrow(cross_data)-2)
		}

	    #determine where phenotypes end and genotypes begin by blanks in row 1
		beginGeno = match(FALSE, is.na(suppressWarnings(as.numeric(cross_data[1,]))))
	
		#if no phenotypes are specified, just take all phenotypes
		if(is.null(pheno_col)){
			pheno_col <- 1:(beginGeno-1)
		}


		#add a check for non-numeric phenotypes
		pheno_classes <- NULL
		for(i in pheno_col){
			pheno_classes <- c(pheno_classes, class(cross_data[,i]))
		}
			
		char_pheno <- which(pheno_classes == "character")
	
		if(length(char_pheno) > 0){
			warning("All phenotypes must be numeric.")
			warning("The following phenotype columns have character values:", paste(colnames(cross_data)[char_pheno], collapse = ", "))
			message("This error can occur if NA's are coded with multiple characters, or if na_strings is mis-specified. Make sure NA coding is consistent throughout the data set and specified correctly with na_strings.")
			stop()
		}
		
		chr <- as.vector(as.matrix(cross_data[1,beginGeno:dim(cross_data)[2]]))
		
		if(check_chr_order){
			x_locale <- grep("x", chr, ignore.case = TRUE)
			y_locale <- grep("y", chr, ignore.case = TRUE)
			m_locale <- grep("m", chr, ignore.case = TRUE)
			just_num_chr <- setdiff(1:length(chr), c(x_locale, y_locale, m_locale))
			consec_chr <- consec_pairs(as.numeric(chr[just_num_chr]))
			order_check <- apply(consec_chr, 1, function(x) x[2] - x[1])
			if(length(which(order_check < 0)) > 0){
				warning("The chromosomes appear to be out of order.\nIt is best to sort the chromosomes before beginning the analysis.")
			}
		}
			
		marker_loc <- as.numeric(cross_data[2,beginGeno:dim(cross_data)[2]])
	
		#take out the genotype matrix
		#It begins in the 3rd row after chromosome numbers and marker locations
		geno <- as.matrix(cross_data[3:dim(cross_data)[1],beginGeno:dim(cross_data)[2]])
	
			
		#if phenotypes are specified as characters, find their
		#locations
		pheno_columns <- get_col_num(cross_data, pheno_col)
			
		#take out the phenotype matrix
		#It begins in the third row and includes the columns
		#specified by the user
		pheno <- as.matrix(cross_data[3:dim(cross_data)[1],pheno_columns])
	
		#convert pheno into a numeric matrix, so we can do 
		#matrix algebra on it later. If there are any non-numeric
		#phenotypes, convert them to numeric 
		pheno <- matrix(apply(pheno, 2, as.numeric), ncol = dim(pheno)[2], byrow = FALSE)
		colnames(pheno) <- colnames(cross_data)[pheno_columns]
		rownames(pheno) <- ind_names
	   
	   if(is.null(geno_col)){
		   	geno_col <- 1:dim(geno)[2]
	   	}
	   	geno_columns <- get_col_num(geno, geno_col)
	  	geno <- geno[,geno_columns] 
	  	chr <- chr[geno_columns]
		marker_loc <- marker_loc[geno_columns]
	   
	    #run a check to see how the genotypes are stored.
		#genotypes can be stored as (0,1,2), ("A","H","B")
		#or as probabilities between 0 and 1
		#if the genotypes are stored as (0,1,2) or 
		#letters, we need to convert them to probabilities
	    
	    found_genotype_mode <- 0 #create a flag to determine whether we have figured out how the genotypes are encoded
		is_multiparent <- 0 #create a flag to identify a multi-parent cross	    

	    genotype_class <- class(cross_data[,beginGeno])
	    if(genotype_class != "character"){
		    all_genotypes <- sort(unique(na.omit(as.numeric(as.matrix(geno))))) #get a vector of all genotypes used
 			if(length(all_genotypes) > 3){
				found_genotype_mode <- 1
				is_multiparent <- 1
				
				fill_array <- function(genotypes){
					geno_mat <- matrix(0, nrow = nrow(geno), ncol = length(all_genotypes))
					for(al in 1:ncol(geno_mat)){
						geno_mat[which(genotypes == all_genotypes[al]),al] <- 1
						}
					return(geno_mat)
					}

				if(verbose){message("I have detected a multi-parent cross\nConverting to array...\n")}
				geno_list <- lapply(1:dim(geno)[2], function(x) fill_array(geno[,x]))

				#fun <- function(x,y) abind(x,y, along = 3)
				#geno_array <- Reduce(fun, geno_list)

				geno_array <- array(NA, dim = c(nrow(geno), length(all_genotypes), ncol(geno)))
				rownames(geno_array) <- ind_names
				colnames(geno_array) <- sort(all_genotypes)
				dimnames(geno_array)[[3]] <- colnames(geno)
				for(m in 1:length(geno_list)){
					geno_array[,,m] <- geno_list[[m]]
					}
		   	}
		}else{
		    all_genotypes <- sort(unique(na.omit(as.vector(as.matrix(geno)))))
		}
		
		#check to see if the genotypes are encoded as letters
	  	if(genotype_class == "character"){
	  		het_present <- grep("H", all_genotypes)
	  		if(length(all_genotypes) > 2 && length(het_present) == 0){
	  			warning("I am detecting more than 2 genotypes: ", paste(all_genotypes, collapse = ", "))
	  			warning("But no H")
	  			stop("Heterozygotes must be coded by H")
	  		}
	  			
	 		
	  		found_genotype_mode <- 1 #indicate that we have found the genotype mode for the file
			#assign baseGeno and notBaseGeno
		    #the baseGeno is assigned a numeric value
		    #of 0. The notBaseGeno is assigned 1
		    #by default we make the first letter
		    #alphabetically the base genotype
		    if(length(all_genotypes) == 3){
				baseGeno <- all_genotypes[all_genotypes != "H"][1]
				notBaseGeno <- all_genotypes[all_genotypes != "H"][2]
				if(verbose){cat("The genotypes are encoded as ", baseGeno, ", H, ", notBaseGeno, "\nConverting to 0, 0.5, 1.\n", sep = "")}
			}else{
				# baseGeno <- all_genotypes[all_genotypes != "H"][1]
				# notBaseGeno <- "H"
				baseGeno <- sort(all_genotypes)[1]
				notBaseGeno <- sort(all_genotypes)[2]
				if(verbose){cat("The genotypes are encoded as ", baseGeno, " and ", notBaseGeno, "\nConverting to 0 and 1.\n", sep = "")}
		 	}
				
		    #turn baseGeno, H, and notBaseGeno to 0, .5, and 1 respectively
			#This function takes in a vector and converts the letters to
			#the appropriate numbers
		
			convert_geno_letter <- function(genotypes){
				genotypes[which(as.character(genotypes) == baseGeno)] <- 0
		    	genotypes[which(as.character(genotypes) == notBaseGeno)] <- 1
		    	if(length(all_genotypes) == 3){
		    		genotypes[which(as.character(genotypes) == "H")] <- 0.5
		    	}else{
		    		genotypes[which(as.character(genotypes) == "H")] <- 1
		    	}
		    	return(as.numeric(genotypes))
		    }
		 	geno <- apply(geno, 2, convert_geno_letter) 
	  	}
	    
		if(!as.logical(is_multiparent)){		
	 		#check to see if the genotypes are encoded as (0, 1, 2)
	 		numeric_test <- which(all_genotypes == 2) #check for 2, since 2 is unique to this encoding
	 		if(length(numeric_test) > 0){
	 			outside_upper_bound <- which(all_genotypes > 2)
	 			outside_lower_bound <- which(all_genotypes < 0)
	 			if(length(outside_upper_bound) > 0 || length(outside_lower_bound) > 0){
	 				stop("Assuming (0,1,2) coding, but I detected genotypes greater than 2 or less than 0.")
	 			}
				if(verbose){cat("The genotypes are encoded as 0, 1, 2.\nConverting to 0, 0.5, 1.\n")}
				found_genotype_mode <- 1 #set the flag indicating we've figured out the encoding
				#turn 0, 1, 2 into 0, 0.5 and 1 respectively
				convert_geno_number <- function(genotypes){
					genotypes[which(as.numeric(genotypes) == 1)] <- 0.5
					genotypes[which(as.numeric(genotypes) == 2)] <- 1
					return(as.numeric(genotypes))
				}
			
				geno <- apply(geno, 2, convert_geno_number) 
		
			}	
		}

	 	#if we still haven't found the genotype mode yet
	 	#check to see if the genotypes are encoded as probabilities
	 	if(found_genotype_mode == 0){
	 		min_geno <- min(all_genotypes); max_geno <- max(all_genotypes) #find the max and min values for genotypes
	 		if(min_geno >= 0 && max_geno <= 1){ #and make sure they are bounded as probabilities are
	 			found_genotype_mode <- 1 #set the flag to indicate we have found the genotype mode
	 			if(verbose){cat("The genotypes are encoded as probabilities.\nNo conversion needed.\n")}
	 			
				#we still need to conver the data frame to a numeric matrix for later compatibility
	 			convert_geno_prob <- function(genotypes){
			        return(as.numeric(genotypes))
				}
		
			 	geno <- apply(geno, 2, convert_geno_prob) 
	
	 		}
	 	}
	 	
	 	
	 	#If after all this, we haven't found the genotype encoding, stop and warn the user
	 	if(found_genotype_mode == 0){
	 		stop("\nGenotypes must be encoded as (0, 1, 2), (A,H,B), or probabilities.\n")
	 	}
	
	
		#take out the sex chromosomes and invariant markers
		x_locale <- grep("X", chr, ignore.case = TRUE)
		if(length(x_locale) > 0){
			if(verbose){cat("\nRemoving markers on the X chromosome")}
			if(is_multiparent){
				geno_array <- geno_array[,,-x_locale]
			}else{
				geno <- geno[,-x_locale]
			}
		chr <- chr[-x_locale]
		marker_loc <- marker_loc[-x_locale]
		}
			
		y_locale <- grep("Y", chr, ignore.case = TRUE)
		if(length(y_locale) > 0){
			if(verbose){cat("\nRemoving markers on the Y chromosome")}
			if(is_multiparent){
				geno_array <- geno_array[,,-y_locale]
			}else{
				geno <- geno[,-y_locale]
			}
			chr <- chr[-y_locale]
			marker_loc <- marker_loc[-y_locale]
		}

		m_locale <- grep("M", chr, ignore.case = TRUE)
		if(length(m_locale) > 0){
			if(verbose){cat("\nRemoving markers on the mitochondrial chromosome")}
			if(is_multiparent){
				geno_array <- geno_array[,,-m_locale]
			}else{
				geno <- geno[,-m_locale]
			}
			chr <- chr[-m_locale]
			marker_loc <- marker_loc[-m_locale]
		}

		#take out markers with only 1 allele
		if(!is_multiparent){
			num_allele <- apply(geno, 2, function(x) length(unique(x)))
			mono_allele <- which(num_allele == 1)
			if(length(mono_allele) > 0){
				if(verbose){cat("\nRemoving invariant markers.\n")}
				geno <- geno[,-mono_allele]
				chr <- chr[-mono_allele]
				marker_loc <- marker_loc[-mono_allele]
			}
		}
	
		na_locale <- which(is.na(geno))
		if(length(na_locale) > 0){
			message("Missing values detected in the genotype matrix.\n\tIf you are planning to use the kinship correction, please use impute.geno() to impute the genotype data.\n")
		}
	
		marker_names <- colnames(geno)
		# colnames(geno) <- 1:dim(geno)[2]
		rownames(geno) <- rownames(pheno)
	
		if(!is_multiparent){
			#scale the genotypes to lie between 0 and 1
			#even if it's a backcross
			geno <- geno/max(geno, na.rm = TRUE)
		}
		
		marker_num <- 1:dim(geno)[2]
	
		if(is_multiparent){
			final_data <- list(pheno, geno_array, chr, marker_names, marker_num, marker_loc)
		}else{
			final_data <- list(pheno, geno, chr, marker_names, marker_num, marker_loc)
		}
		names(final_data) <- c("pheno", "geno", "chromosome", "marker_names", "marker_num", "marker_location")

		if(verbose){     
			cat("Read in the following data:\n")
			cat("\t-", dim(pheno)[1], "individuals -\n")
			cat("\t-", dim(geno)[2], "markers -\n")
			cat("\t-", dim(pheno)[2], "phenotypes -\n")
		}
	
	
	    return(final_data)
}
