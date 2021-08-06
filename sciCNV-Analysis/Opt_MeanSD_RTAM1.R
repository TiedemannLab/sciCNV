

##########################################################################################
######                        Optimizing Function in RTAM1                         ####### 
###### (Detects the optimal number of top-expressed genes for scaling a data-set)  #######
######   Tiedemann Lab - Princess Margaret Cancer Centre, University of Toronto    #######
######                       copyright@AliMahdipourShirayeh                        #######
##########################################################################################

# Please refer to the reference and supplemental materials described in the README for additional details.
#
# Definitions:
# Normalized_log: the matrix of Log2(.+1)-normalized data 
# Order_Matrix: the matrix of original order of gene expressions per cell
# nGene: number of expressed genes
# Min_nGene: minimum number of expressed genes across entire population
# gene_cutoff: the cutoff to calculate the summation of top (nonzero) gene expressions

Opt_MeanSD_RTAM1 <- function(Normalized_log, 
                             Order_Matrix, 
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
  
  rownames(AA)<- rownames(matrix_log)
  colnames(AA)<- colnames(general1)
  
  # ----------------------------------------------------------------
  # Step 2: RTAM1 main step, adjusting cellular gene expression
  # ----------------------------------------------------------------
  

  #Asigning R_ij matrix to L_ij or AA matrix of ordered gene expressios in each cell
  R_mtrx <- matrix(0,ncol=ncol(AA),nrow=nrow(AA))
  
  for( j in 1:ncol(AA)){
    CC <- NULL
    LLL<- nGene[1,j]
    for(i in 1:LLL){
      
      if(!(i %in% CC) ){
        
        if(i == nrow(AA)){
          R_mtrx[i,j] <- i
          
        }else if(i < nrow(AA)){
          
          if(AA[i,j]>AA[i+1,j] ){
            R_mtrx[i,j] <- i
            
          }else if(AA[i,j]==AA[i+1,j]){ 
            Equal_list <- as.matrix(which( AA[,j] == AA[i,j] ) )
            Length <- length( as.matrix(Equal_list) )
            Equal_value <- sum( seq(min(Equal_list),max(Equal_list),1) )/Length
            
            R_mtrx[ Equal_list , j ] <- Equal_value
            
            if( (i + Length-1) <= nrow(AA)){
              CC <- seq(min(Equal_list),max(Equal_list),1)
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  
  
  
  
  ############################################# 
  
  
  Y <- seq(1, gene_cutoff, 1)
  BB <- matrix(0, ncol= ncol(AA) , nrow= gene_cutoff )
  BB <- as.matrix(AA[Y, ])
  
  
  #------------------
  pSum <- t(as.matrix(colSums(BB)))
  MIN_intesnse <- mean( pSum)  
  ratio <- matrix(0, ncol=ncol(AA), nrow=1)
  
  for(j in 1:ncol(AA)){
    
    ratio[1 ,j] <- (MIN_intesnse - pSum[1,j])/( (gene_cutoff + 3/2)*log2(gene_cutoff + 3/2) 
                    - aaa/log(2)-3/2*log2(3/2))  
    
  }
  
  
  h_mtrx <- matrix(0, ncol = ncol(AA), nrow = nrow(AA))
  
  for(j in 1:ncol(R_mtrx)){
    for(i in 1:nrow(R_mtrx)){ 
      h_mtrx[i ,j] <- ratio[1 ,j]*log2( R_mtrx[i,j] + 1 )
    }
  }
  
  LS1 <- AA + h_mtrx
  Scaled_Normalized_log <- matrix(0, ncol=ncol(AA), nrow=nrow(AA))
  
  for(j in 1:ncol(AA)){
    Scaled_Normalized_log[ as.matrix(Order_Matrix[, j]) , j] <- LS[  , j]
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
  
  # Mean of A% common genes
  Mean_Aper <- as.matrix(mean(MEAN_comm))  
  Sd_Aper <- as.matrix( sd(MEAN_comm))
  
  return((Sd_Aper/Mean_Aper)^2)
}




