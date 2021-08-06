

#######################################################################################
######                      Optimizing Function in RTAM2                        ####### 
######(Detects the optimal number of top-expressed genes for scaling a data-set)#######
######  Tiedemann Lab - Princess Margaret Cancer centre, University of Toronto  #######
######                       copyright@AliMahdipourShirayeh                     #######
#######################################################################################

# Please refer to the reference and supplemental materials described in the README for additional details.
#
# Definitions:
# Normalized_log: the matrix of Log2(.+1)-normalized data 
# Order_Matrix: the matrix of original order of gene expressions per cell
# G_mtrx: TPM-normalized matrix of the data
# nGene: number of expressed genes
# Min_nGene: minimum number of expressed genes across entire population
# gene_cutoff: the cutoff to calculate the summation of top (nonzero) gene expressions

Opt_MeanSD_RTAM2 <- function(Normalized_log, 
                             Order_Matrix, 
                             G_mtrx, 
                             nGene, 
                             Min_nGene, 
                             gene_cutoff
){    
  

  L <- dim(Normalized_log)[1]
  W <- dim(Normalized_log)[2]
  
  AA <- matrix(0, ncol=W, nrow= L)
  
  
  for(w in 1:W){
    
    SS <-  as.matrix(sort(Normalized_log[ ,w], decreasing=TRUE))
    LLS <- dim(SS)[1]
    
    for(i in 1:LLS){
      AA[ i,w] <- SS[i,1]
    }
    
  }
  
  colnames(AA)<- colnames(Normalized_log)
  
  # ---------------------------------------------------------------------------
  # Step 2: RTAM2-Lin main step, adjusts cellular gene expression
  # ---------------------------------------------------------------------------
  
  #Asigning G_ij matrix as a ranked matrix of colum-sum normalized data
  GG <- matrix(0, ncol=ncol(AA), nrow= nrow(AA))
  
  for(w in 1:ncol(G_mtrx)){
    
    SS <-  as.matrix(sort(G_mtrx[ ,w], decreasing=TRUE))
    
    LLS <- dim(SS)[1]
    for(i in 1:LLS){
      
      GG[ i,w] <- SS[i,1]
      
    }
  }
  
  colnames(GG)<- colnames(AA)
  

  #Asigning G_ij matrix as a ranked matrix of colum-sum normalized data
  GG <- matrix(0, ncol=ncol(AA), nrow= nrow(AA))
  
  for(w in 1:ncol(G_mtrx)){
    
    SS <-  as.matrix(sort(G_mtrx[ ,w], decreasing=TRUE))
    
    LLS <- dim(SS)[1]
    for(i in 1:LLS){
      
      GG[ i,w] <- SS[i,1]
      
    }
  }
  
  colnames(GG)<- colnames(AA)
  
  
 ############################################# 

  
  Y <- seq(1, gene_cutoff, 1)
  BB <- matrix(0, ncol= ncol(AA) , nrow= gene_cutoff )
  BB <- as.matrix(AA[Y, ])
  
  #------------------
  pSum <- t(as.matrix(colSums(BB)))
  MIN_intesnse <- mean( pSum)  
  Scaled_Normalized_log <- matrix(0, ncol=ncol(AA), nrow= nrow(AA))
  HH <- matrix(0, ncol=ncol(AA), nrow= 1)
  #------------------
  
  for(j in 1:ncol(G_mtrx)){
    hh <- 0
    for( i in 1:gene_cutoff){ 
      hh <- hh + log2(GG[ 1,j] - GG[ i,j] +1)   
    }
    HH[1 , j] <- hh
  }
  
  
  h_mtrx <- matrix(0, ncol = ncol(AA), nrow = nrow(AA))
  
  for(j in 1:ncol(AA)){
    for(i in 1:nrow(as.matrix(GG[ , j][GG[ , j]>0]))){ 
      h_mtrx[i ,j] <- (MIN_intesnse - pSum[1,j])*(log2(GG[ 1,j] - GG[ i,j]+1))/(HH[1,j] )
    }
  }
  
  LS1 <- AA + h_mtrx
  
  Scaled_Normalized_log <- matrix(0, ncol=ncol(AA), nrow=nrow(AA))
  
  for(j in 1:ncol(AA)){
    Scaled_Normalized_log[ as.matrix(Order_Matrix[, j]) , j] <- LS1[  , j]
  }
  

  Scaled_Normalized_log <- as.matrix(Scaled_Normalized_log)
  comm.expr <-   as.matrix(Scaled_Normalized_log[which(rowSums(Scaled_Normalized_log[ ,seq(1, ncol(AA), 1)]!=0) > 0.95*ncol(Scaled_Normalized_log) ), ] )
  MEAN_comm <- matrix(0, ncol=ncol(AA), nrow=1)
  for(j in 1:ncol(AA)){
    
    if( sum(comm.expr[j]) > 0){
      MEAN_comm[1,j] <- mean(comm.expr[j][comm.expr[j]>0]) 
    } else {
      MEAN_comm[1,j] <- 0
    }
    
  }
  
  # Mean and standard deviation of A% common genes
  Mean_Aper <- as.matrix(mean(MEAN_comm))  
  Sd_Aper <- as.matrix( sd(MEAN_comm))
  
  return((Sd_Aper/Mean_Aper)^2)
  
}




