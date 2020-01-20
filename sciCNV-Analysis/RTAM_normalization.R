

##################################################################################################
######                           RTAM Nomalization Function                                ####### 
######                               (for scRNA-seq data)                                  #######
######     Tiedemann Lab - Princess Margaret Cancer Centre, University of Toronto          #######
######                        copyright@AliMahdipourShirayeh                               #######
##################################################################################################

# Please see the reference and supplemental materials described in the README file.

# The raw data matrix has genes as rows and cells as columns
# Genes are annotated with their genomic location (start and end locations)
# Order_Matrix is the original order of genes in the raw data
# MinNumGenes is a threshold used to exclude poor QC cells that have few detectable genes (default= minimum of 250 genes).
# OptNumGenes is the number of genes that when employed for RTAM normalization will optimally minimize the variance across cells of gene-set average expression


RTAM_normalization <- function(mat,           # Raw scRNA-seq data 
                               method,        # either RTAM1 or RTAM 2
                               Min_nGn,       # minimum number of genes per cell, below which cells are excluded. 
                                              # Defines the smallest least complex allowable cell. (default = 250 genes per cell)
                               Optimizing     # True or FALSE options specify whether to run an optional optimization code
){
  
  # argument validations
  #if( is.na(mat) ){
  #  stop("Error, the given matrix is not in the accurate format")
  #}
  
  if( is.na(method) ){
    method <- c("RTAM2")
    flog.info(sprintf("The considered normalization method is %s",method))
  }
  
  
  if(! (method %in% c("RTAM1", "RTAM2")) ){
    stop( "The chosen normalization method is not valid.")
  }
  
  #if( (Optimizing=="TRUE") & (! is.na(Nt)) 
  #    ){
  #  stop("Error, Nt (the number of top-ranked genes) can be specified only if Optimizing = FALSE")
  #}
  
  if( is.na(Optimizing) ){
    Optimizing = "FALSE"
  }
  
  if( Optimizing == "FALSE" ){
    Nt <- Min_nGn
  } 
  
  #if ( (!is.na(Nt) ) &  ( Nt > Min_nGn ) ){
  #  stop("Error, Nt is required to be not greater than Min_nGn")
  #}
  
  general <- as.matrix(mat) #[ , -1])
  rownames(general) <- rownames(mat) #mat[ , 1]
  
  # -------------------------------------------------------------------
  # Initial "TPM" normalization
  # -------------------------------------------------------------------
  
  matrix <- as.matrix(general)
  norm_mat <- matrix %*% diag(1/colSums(matrix))
  norm_mat <- norm_mat*1e5
  normlog_mat <- log2(norm_mat+1) 
  rownames(normlog_mat) <- rownames(general) 
  colnames(normlog_mat) <- colnames(general)
  
  # -------------------------------------------------------------------
  # Data Organization
  # -------------------------------------------------------------------
  
  ## Saving the original order of genes
  Order_Matrix <- matrix(0, ncol = ncol(normlog_mat), nrow = nrow(normlog_mat))
  for(l in 1:ncol(normlog_mat)){
    Order_Matrix[ , l] <- as.matrix(order(normlog_mat[ , l] , decreasing=TRUE) )
  }
  
  ## Ranking genes by expression values in each cell
  nGene <- matrix(0, ncol=ncol(normlog_mat), nrow= 1)
  nonzero <- function(x) sum(x != 0)
  
  
  ## Ranked matrix is called ranked_mat
  ranked_mat <- matrix(0, ncol=ncol(normlog_mat), nrow= nrow(normlog_mat))
  for(w in 1:ncol(normlog_mat)){
    
    SS <-  as.matrix(sort(normlog_mat[ ,w], decreasing=TRUE))
    LS <- dim(SS)[1]
    for(i in 1:LS){
      ranked_mat[ i,w] <- SS[i,1]
      
    }
    
  }
  colnames(ranked_mat) <- colnames(normlog_mat)
  
  # -------------------------------------------------------------------
  # RTAM normalization
  # -------------------------------------------------------------------
  
  for(i in 1:ncol(ranked_mat)){ nGene[ ,i] <- nonzero(ranked_mat[,i])  }
  colnames(nGene) <- colnames(normlog_mat)
  
  if ( method == c("RTAM2")){
    G_mtrx <- as.matrix(norm_mat)
  }
  
  if( is.na(Min_nGn) ){
    Min_nGene <- 250
  } else {
    Min_nGene <- max(250, min(nGene), Min_nGn) 
  }
  
  if( Min_nGene > min(nGene) ){
    
    KK <- t(as.matrix(which(t(as.matrix(nGene)) < Min_nGene) ))
    ##
    ranked_mat <- as.matrix( ranked_mat[ , -KK ] )
    colnames(ranked_mat) <- colnames(general[, -KK])
    ##
    normlog_mat <- as.matrix( normlog_mat[ , -KK ] )
    rownames(normlog_mat) <- rownames(general) 
    colnames(normlog_mat) <- colnames(general[, -KK])
    ##
    Order_Matrix <-  as.matrix( Order_Matrix[ , -KK ] )
    if ( method == c("RTAM2")){
      G_mtrx <- as.matrix( G_mtrx[ , -KK ] )
    }
  }
  
  ############## Ranked matrix of colum-sum normalized data for RTAM2 method
  if ( method == c("RTAM2")){
    
    ## Asigning GG matrix as a ranked matrix of colum-sum normalized data
    GG <- matrix(0, ncol=ncol(ranked_mat), nrow= nrow(ranked_mat))
    for(w in 1:ncol(G_mtrx)){
      SS <-  as.matrix(sort(G_mtrx[ ,w], decreasing=TRUE))
      LLS <- dim(SS)[1]
      for(i in 1:LLS){
        GG[ i,w] <- SS[i,1]
      }
    }
    colnames(GG)<- colnames(ranked_mat)
    
    
    ############# Asigning R_mat to ranked matrix of ordered gene expressions in each cell
  } else if ( method == c("RTAM1")){
    
    
    R_mtrx <- matrix(0,ncol=ncol(ranked_mat),nrow=nrow(ranked_mat))
    for( j in 1:ncol(ranked_mat)){
      CC <- NULL
      LLL <- nGene[1,j]
      for(i in 1:LLL){
        if(!(i %in% CC) ){
          if(i == nrow(ranked_mat)){
            R_mtrx[i,j] <- i
          }else if(i < nrow(ranked_mat)){
            if(ranked_mat[i,j]>ranked_mat[i+1,j] ){
              R_mtrx[i,j] <- i
            }else if(ranked_mat[i,j]==ranked_mat[i+1,j]){ 
              Equal_list <- as.matrix(which( ranked_mat[,j] == ranked_mat[i,j] ) )
              Length <- length( as.matrix(Equal_list) )
              Equal_value <- sum( seq(min(Equal_list),max(Equal_list),1) )/Length
              R_mtrx[ Equal_list , j ] <- Equal_value
              if( (i + Length-1) <= nrow(ranked_mat)){
                CC <- seq(min(Equal_list),max(Equal_list),1)
              }
            }
          }
        }
      }
    }
    
    
  }  
  
  ## Finding the optimal number of genes for RTAM normalization (to yield minimum variance in average geneset expression across the matrix)
  if( Optimizing == "TRUE" ) {
      
      SEQ <- seq(Min_nGene-floor(Min_nGene/2), Min_nGene, 5)
      Dispersion <- rep("NA", length(SEQ))
    
        if ( method == "RTAM2"){
          Dispersion[1:length(SEQ)] <- lapply( 1:length(SEQ) , function(i){
                    Opt_MeanSD_RTAM2(Normalized_log = normlog_mat, 
                                     Order_Matrix = Order_Matrix, 
                                     G_mtrx =  G_mtrx, 
                                     nGene =  nGene, 
                                     Min_nGene = Min_nGene, 
                                     gene_cutoff = SEQ[i])
          })
          
        }else{
          
          Dispersion[1:length(SEQ)] <- lapply( 1:length(SEQ) , function(i){
                    Opt_MeanSD_RTAM1(Normalized_log = normlog_mat, 
                                     Order_Matrix = Order_Matrix, 
                                     nGene = nGene, 
                                     Min_nGene = Min_nGene, 
                                     gene_cutoff = SEQ[i])
          })
        
        }
      
      Dispersion <- as.numeric(Dispersion)
      Nt <- SEQ[which(Dispersion == min(Dispersion))]
      
  } else {
    
    Nt <- Min_nGene
    
  }
  
  
  # --------------------------------------------------------------------
  # RTAM2 - differential normalization main step
  # -------------------------------------------------------------------
  if ( method == c("RTAM2")){
    
    ## Defining top_mat as the sub-matrix of top-expressed genes per cell
    Y <- seq(1, Nt, 1)
    top_mat <- matrix(0, ncol = ncol(ranked_mat) , nrow= Nt )
    top_mat <- as.matrix(ranked_mat[Y, ])
    
    #------------------
    
    pSum <- matrix(0, ncol = ncol(ranked_mat), nrow = 1)
    ratio <- matrix(0, ncol = ncol(ranked_mat), nrow = 1)
    Scaled_Normalized_log <- matrix(0, ncol = ncol(ranked_mat), nrow = nrow(ranked_mat))
    HH <- matrix(0, ncol = ncol(ranked_mat), nrow = 1)
    
    pSum[1,] <- as.matrix( colSums(top_mat))
    
    MIN_intesnse <- mean( pSum)   
    
    for(j in 1:ncol(G_mtrx)){
      hh <- 0
      for( i in 1:Nt){ 
        hh <- hh + (GG[ 1,j] - GG[ i,j] +1)   
      }
      HH[1 , j] <- hh
    }
    
    
    h_mtrx <- matrix(0, ncol = ncol(ranked_mat), nrow = nrow(ranked_mat))
    
    for(j in 1:ncol(ranked_mat)){
      for(i in 1:nrow(as.matrix(GG[ , j][GG[ , j]>0]))){ 
        h_mtrx[i ,j] <- (MIN_intesnse - pSum[1,j])*(GG[ 1,j] - GG[ i,j]+1)/(HH[1,j] )
      }
    }
    
    
    #------------------------------------------
    ## Normalization limitations to prevent: 
    ## (i) changes in the expression-based ranking of expressed genes
    ## (ii) negative expression values
    
    for(j in 1:ncol(ranked_mat)){
      
      if( h_mtrx[ Nt, j] >= (GG[ 1,j] + 1 - GG[ Nt, j])/((GG[ 1,j] + 1)*log(2, base=exp(1)))){
        h_mtrx[,j] <- h_mtrx[ ,j]*(GG[ 1,j] + 1 - GG[ Nt, j])/((GG[ 1,j] + 1)*log(2, 
                                                                                base=exp(1))*h_mtrx[ Nt, j])
      }
      
      h_mtrx_min <- min(h_mtrx[ ,j])
      GG_min <- min(log2(GG[ , j ][GG[ , j ] > 0 ] + 1))
      
      if( h_mtrx_min <= - GG_min ){
        h_mtrx[,j] <- h_mtrx[ ,j]*(-GG_min/h_mtrx_min )
      }
      
    }
    
    
  } 
  # ----------------------------------------------------------------
  # RTAM1 -differential normalization main step
  # ----------------------------------------------------------------
  else if ( method == c("RTAM1")){
    
    ## Using the optimized or specified number of top-expressed genes
    
    
    Y <- seq(1, Nt, 1)
    top_mat <- matrix(0, ncol = ncol(ranked_mat) , nrow = Nt )
    top_mat <- as.matrix(ranked_mat[Y, ])
    
    #------------------
    pSum <- matrix(0, ncol = ncol(ranked_mat), nrow = 1)
    ratio <- matrix(0, ncol = ncol(ranked_mat), nrow = 1)
    Scaled_Normalized_log <- matrix(0, ncol = ncol(ranked_mat), nrow = nrow(ranked_mat))
    
    pSum[1,] <- as.matrix( colSums(top_mat))
    
    MIN_intesnse <- mean( pSum)  
    
    for(j in 1:ncol(ranked_mat)){
      ratio[1 ,j] <- ( MIN_intesnse - pSum[1,j] )/( (Nt + 3/2)*log2(Nt + 3/2) - Nt/log(2)-3/2*log2(3/2) )  
    }
    
    h_mtrx <- matrix(0, ncol = ncol(ranked_mat), nrow = nrow(ranked_mat))
    
    for(j in 1:ncol(ranked_mat)){
      for(i in 1:nrow(ranked_mat)){ 
        h_mtrx[i ,j] <- ratio[1 ,j]*log2( R_mtrx[i,j] + 1 )
        
      }
    }
    
    
    
  }
  
  ############################################
  
  ## Normalizing the data by adding the calculated corrections to each gene expression value
  LS <- ranked_mat + h_mtrx
  
  Scaled_Normalized_log <- matrix(0, ncol=ncol(ranked_mat), nrow=nrow(ranked_mat))
  
  for(j in 1:ncol(ranked_mat)){
    Scaled_Normalized_log[ as.matrix(Order_Matrix[, j]) , j] <- LS[  , j]
  }
  
  rownames(Scaled_Normalized_log) <- rownames(normlog_mat)
  colnames(Scaled_Normalized_log) <- colnames(ranked_mat)
  
  
  return(Scaled_Normalized_log)
  
  
}




