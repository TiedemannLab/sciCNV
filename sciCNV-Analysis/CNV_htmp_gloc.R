

#######################################################################################
######                           CNV_htmp_gloc function                         ####### 
###### Generates heatmap of sciCNV profiles, plotted by Genomic location (gloc) #######
######  Tiedemann Lab - Princess Margaret Cancer Centre, University of Toronto  #######
######                       copyright@AliMahdipourShirayeh                     #######
#######################################################################################


# Please refer to the reference and supplemental materials described in the README for additional details.

# principal variables:
# clustering: TRUE/FALSE variable specifying whether to cluster cells based on their CNV similarities. Default is "FALSE"
# clustering.type: Variable specifying the clustering method to be used in generating the heatmap. Possible options are "pearson", 
#        "euclidean", " spearman", ... "original" (retains the original cell order without unsupervised clustering).
#         Default is "pearson". Only enabled when clustering = "TRUE"
# sorting: TRUE/FALSE variable specifying whether to sort cells based on their tumor CNV score from the largest to smallest tumor scores. Default is FALSE.

# Definitions:
# CNV.mat2: copy number variation matrix 
# Gen.Loc: Genomic Location matrix with list of genes, their associated chromosome numbers, start and end locations
# CNVscore: is the tumor CNV score matrix for all cells (possibly ranked within clusters). Only used when sorting.clusters = TRUE.
# cluster.lines: is a list of values which can be used to separate cell clusters within the population; only used if multiple clusters present
# break.gloc: is a set of values each defines a vertical line that separates chromosomes
# No.test: number of test cells included in the data; can be used to delineate the populations of test and control cells in the heatmap


CNV_htmp_gloc <- function(CNV.mat2,
                          clustering = FALSE,        # TRUE or FALSE option 
                          clustering.type = c("pearson", "kendall", "spearman"),   # "pearson", "kendalln", "spearman", default: "pearson"
                          sorting = TRUE,            # TRUE or FALSE
                          CNVscore = NULL,                # Only exists when sorting = TRUE 
                          cluster.lines = NULL,   
                          break.gloc = break.gloc, 
                          No.test ){
  
  ## argument validation
  if ( missing(sorting) ){
    sorting <- FALSE
  }
  
  if  ( is.null(clustering) ){
    clustering <- FALSE
  }
  if ( clustering == "TRUE" ){
    if ( ! clustering.type %in% c("pearson", "kendall", "spearman")  ){
      stop("Please choose a proper method for unsupervised clustering.")
    } else if ( is.null(clustering.type) ){
      clustering.type <- c("pearson")
    }
    if ( sorting == "TRUE" ){
      stop("Sorting can occur when clustering = FALSE.")
    }
  } else if (clustering == "FALSE"){
    if ( (sorting == "TRUE") & Reduce("|", is.na(CNVscore)) ){
      stop("Please insert a list of CNV-scores")
    }
    if ( (sorting == "FALSE") & Reduce("|", ! is.na(CNVscore)) ){
      stop("Please change sorting status to TRUE to sort data based on CNVscore vector.")
    }
  }
  
  if  ( (is.null(No.test)) & Reduce("|", is.null(cluster.lines)) ){
    stop("Please insert the number of test cells (No.test).")
  }
  
  if ( Reduce("|", is.null(cluster.lines)) ){
    cluster.lines <- c(0, nrow(CNV.mat2) - No.test, nrow(CNV.mat2))
  }
  
  if ( Reduce("|", is.null(break.gloc)) ){
    break.gloc <- c(0, ncol(CNV.mat2))
  }
  
  
  
  ##### sorting of cells within clusters, based on tumor CNV scores, from the largest to the smallest (if applicable)
  
  if ( sorting == TRUE ){
    tst.score <- sort(CNVscore[1, 1:No.test] , decreasing=TRUE)     #MMPCs
    ctrl.score <- sort(CNVscore[1, (No.test+1):ncol(CNVscore)] , decreasing=TRUE)  #NBCs
    ranked.col <- as.matrix( c(colnames(t(as.matrix(ctrl.score))), colnames(t(as.matrix(tst.score))  )) )
    CNV.mat1 <- as.matrix( CNV.mat2[match(ranked.col, rownames(CNV.mat2)), ])
    rownames(CNV.mat1) <-  ranked.col     
    colnames(CNV.mat1) <- rownames(M_NF)
    
  } else if ( clustering == TRUE ){
    
    if ( is.na(clustering.type) ){
      CNV.mat.tst <- as.matrix(CNV.mat2[1:No.test, ])
      hclst <- hclust(as.dist(1-cor( t(CNV.mat.tst), method =  "pearson")), method = "ward.D2")
      hclst.lables <- hclst$labels[hclst$order]
      CNV.mat.clustered <- CNV.mat.tst[hclst.lables , ]
      
    } else if ( ! missing(clustering.type) ){
      CNV.mat.tst <- as.matrix(CNV.mat2[1:No.test, ])
      hclst <- hclust(as.dist(1-cor( t(CNV.mat.tst), method = clustering.type)), method = "ward.D2")
      hclst.lables <- hclst$labels[hclst$order]
      CNV.mat.clustered <- CNV.mat.tst[hclst.lables , ]
    }
    
    CNV.mat1 <- rbind( as.matrix(CNV.mat2[(No.test+1):nrow(CNV.mat2), ]) ,   as.matrix(CNV.mat.clustered)  )
    rownames(CNV.mat1) <-  c(rownames(as.matrix(CNV.mat2[(No.test+1):nrow(CNV.mat2), ])), hclst.lables)     
    colnames(CNV.mat1) <- colnames(CNV.mat2) 
    
  } else if ( (clustering == "FALSE" ) & ( sorting == "FALSE")){
    CNV.mat1 <- rbind(as.matrix(CNV.mat2[(No.test):nrow(CNV.mat2), ]) , as.matrix(CNV.mat2[1:No.test, ]) )
    
  }
  
  ROWlist <- rownames(CNV.mat1)
  COLlist <- colnames(CNV.mat1) 
  
  ## The original list of chromosomes associated with genes in iCNV-matrix
  Gen.Loc <- read.table( "../sciCNV/Dataset/10XGenomics_gen_pos_GRCh38-1.2.0.txt", sep='\t', header=TRUE)
  Specific_genes <- which( as.matrix(Gen.Loc)[, 1]   %in% colnames(CNV.mat2))
  M_sample <-  as.matrix(Gen.Loc[Specific_genes, ])
  
  ## expanding expressions towrds 0 or 1 providing their expression is less than or bigger to 0.5
  LL1 <- 0.5
  LL2 <- 1.5
  LL3 <- 2.5
  CNV.mat11 <- matrix(0, ncol = ncol(CNV.mat1), nrow = nrow(CNV.mat1))
  
  for (w in 1: ncol(CNV.mat1)){
    for(l in 1: nrow(CNV.mat1)){
      
      if( abs(CNV.mat1[ l, w]) <  LL1 ){
        CNV.mat11[ l, w] <- sign(CNV.mat1[ l, w])*( CNV.mat1[ l, w] )^2
      } else if((abs(CNV.mat1[ l, w]) >=  LL1) & (abs(CNV.mat1[ l, w]) <=  LL2)){
        CNV.mat11[ l, w] <- sign(CNV.mat1[ l, w])*sqrt( abs(CNV.mat1[ l, w]) )
      } else if(abs(CNV.mat1[ l, w]) >  LL2 ){
        CNV.mat11[ l, w] <- 2*sign(CNV.mat1[ l, w])*sqrt( abs(CNV.mat1[ l, w]) )
      }
      
    }
  }
  #-------
  LL <- 1
  TT <- 0.5
  CNV.mat1[which(CNV.mat11 >  LL)] <-  LL
  CNV.mat1[which( (CNV.mat11 <  TT) & (CNV.mat11 >  -TT) )] <- 0.0
  CNV.mat1[which(CNV.mat11 <  -LL)] <-  -LL
  
  CNV.mat11 <- as.matrix(CNV.mat11)
  #rownames(CNV.mat) <- matrix(NA, ncol=1, nrow=nrow(CNV.mat1))
  
  ##################################################
  ## Heatmap of CNV-curves against genomic locations
  ##################################################
  
  # Uploading the largest list of genes with chr number, start and end
  MM <- Gen.Loc
  MM1 <- data.frame(MM[,-1]); rownames(MM1) <- MM[,1]
  
  ############ Finding the length of each Chromosome
  M_origin <- MM
  
  ## number of segments on the genome
  No_Intrvl <- 1000
  
  ############ Finding the length of each Chromosome
  
  minn <- matrix(0, ncol=24, nrow=1)
  maxx <- matrix(0, ncol=24, nrow=1)
  
  for(i in 1: 22){ 
    minn[1,i] <- as.numeric(min(M_origin[which(as.matrix(M_origin[, 2]) == i) , 3]) )
    maxx[1,i] <- as.numeric(max(M_origin[which(as.matrix(M_origin[, 2]) == i) , 3]) )
  }
  minn[1,23] <- as.numeric(min(M_origin[which(M_origin[, 2] == "X") , 3]) )
  maxx[1,23] <- as.numeric(max(M_origin[which(as.matrix(M_origin[, 2]) == "X") , 3]) )
  
  minn[1,24] <- as.numeric(min(M_origin[which(as.matrix(M_origin[, 2]) == "Y") , 3]) )
  maxx[1,24] <- as.numeric(max(M_origin[which(as.matrix(M_origin[, 2]) == "Y") , 3]) )
  
  Total_length <- sum( (maxx[1,1:24]-minn[1,1:24]) + 1 ) 
  Seg_size <- ceiling( Total_length/No_Intrvl )
  Start_End_Chr <- rbind(minn, maxx)
  
  ###########
  
  MaxNo_intervals <- ceiling(max((maxx[1,1:24]-minn[1,1:24])+ 1)/Seg_size)
  POINTS <- matrix(0, nrow = 24, ncol = MaxNo_intervals)
  Min_chr <- rep(0,25)
  
  Min_chr[1] <- 0
  for(i in 1:24){
    StartEnd <- seq(minn[1,i] , maxx[1,i] , Seg_size)
    POINTS[i, 1:length(StartEnd) ]   <-  t(as.matrix(StartEnd  )) 
    Min_chr[i+1]  <-  maxx[1,i]              
  }
  
  require(Matrix)
  Length_POINTS <- rep(0,24)
  for(j in 1:24){
    Length_POINTS[j] <- nnzero(POINTS[j, ])
  }
  
  ALL_POINTS <- length(which(POINTS != 0 ))
  
  ###########
  
  CNV.mat4 <- matrix(0, nrow = nrow(CNV.mat11), ncol = ALL_POINTS)
  Minnn <- matrix(0, nrow=1,ncol=24)
  Maxxx <- matrix(0, nrow=1,ncol=24)
  
  for(i in 1:nrow(CNV.mat4)){
    for(z in 1:22){
      for( k in 1:(length(which(POINTS[z, ]>0))-1) ){
        X <- min(which(as.matrix(M_sample[, 2]) == z)):max(which(as.matrix(M_sample[, 2]) == z))
        ExistExpr <- which( (M_sample[ X ,3] >= POINTS[z, k]) 
                            & (M_sample[ X ,3] <=  POINTS[z, k+1]) ) 
        if( z == 1){
          if( length(ExistExpr) == 0){
            CNV.mat4[i,k] <- NA
          }else if( length(ExistExpr) == 1){
            CNV.mat4[i,k] <- CNV.mat11[i, ExistExpr]
          }else{
            CNV.mat4[i,k] <- mean( CNV.mat11[i, ExistExpr])
          }
          
        }else {
          PrePoint_No <-  sum(Length_POINTS[1:(z-1)]) 
          if( length(ExistExpr) == 0){
            CNV.mat4[i,k + PrePoint_No] <- NA
          }else if( length(ExistExpr) == 1){
            CNV.mat4[i,k + PrePoint_No] <- CNV.mat11[i, ExistExpr + min(X) - 1 ]
          }else{
            CNV.mat4[i,k + PrePoint_No] <- mean( CNV.mat11[i, ExistExpr + min(X) - 1 ])
          }
        }
      }
    }
    
    ##########
    
    for( k in 1:(length(which(POINTS[23, ]>0))-1) ){
      if( length(which(as.matrix(M_sample[, 2]) == "X")) != 0 ){
        X <- min(which(as.matrix(M_sample[, 2]) == "X")):max(which(as.matrix(M_sample[, 2]) == "X"))
        ExistExpr <- which( (M_sample[ X ,3] >= POINTS[23, k]) 
                            & (M_sample[ X ,3] <=  POINTS[23, k+1]) ) 
        PrePoint_No <- sum(Length_POINTS[1:22]) 
        if( length(ExistExpr) == 0){
          CNV.mat4[i,k + PrePoint_No] <- NA
        }else if( length(ExistExpr) == 1){
          CNV.mat4[i,k + PrePoint_No] <- CNV.mat11[i, ExistExpr + min(X) - 1 ]
        }else{
          CNV.mat4[i,k + PrePoint_No] <- mean( CNV.mat11[i, ExistExpr + min(X) - 1 ])
        }
      }
    }
    
    ##########
    
    for( k in 1:(length(which(POINTS[24, ]>0))-1) ){
      if( length(which(as.matrix(M_sample[, 2]) == "Y")) != 0 ){
        X <- min(which(as.matrix(M_sample[, 2]) == "Y")):max(which(as.matrix(M_sample[, 2]) == "Y"))
        ExistExpr <- which( (M_sample[ X ,3] >= POINTS[24, k]) 
                            & (M_sample[ X ,3] <=  POINTS[24, k+1]) )
        PrePoint_No <- sum(Length_POINTS[1:23]) 
        if( length(ExistExpr) == 0){
          CNV.mat4[i,k + PrePoint_No] <- NA
        }else if( length(ExistExpr) == 1){
          CNV.mat4[i,k + PrePoint_No] <- CNV.mat11[i, ExistExpr + min(X) - 1  ]
        }else{
          CNV.mat4[i,k + PrePoint_No] <- mean( CNV.mat11[i, ExistExpr + min(X) - 1 ])
        }
      }
    }
  }
  
  
  
  ######### Interpolations
  
  for(k in 1:ALL_POINTS){
    
    if( (k == 1) && is.na(CNV.mat4[1,k])){
      BEF <- 1  
      AFT <- 1 + min(which( CNV.mat4[ 1 , 2:ncol(CNV.mat4)] != "NA")) 
      BEF_Val <- CNV.mat4[, BEF]
      AFT_Val <- CNV.mat4[, AFT]
      #
      CNV.mat4[ , k] <- AFT_Val    #(AFT_Val + BEF_Val)/(AFT-BEF)
      #
    } else if( (k == ncol(CNV.mat4)) && is.na(CNV.mat4[1,k]) ){
      BEF <- max( which(  CNV.mat4[ 1 ,1:(k-1)] != "NA")) 
      AFT <- ncol(CNV.mat4) #min(which( MB4[i,(k+1):ncol(MB4)] != "NA" )  ) 
      BEF_Val <- CNV.mat4[, BEF]
      AFT_Val <- CNV.mat4[, AFT]
      #
      CNV.mat4[ , (k-1)+l] <- BEF_Val    #(AFT_Val + BEF_Val)/(AFT-BEF)
      #
    } else if( is.na(CNV.mat4[1,k] )){
      BEF <-  max( which( CNV.mat4[1,1:(k-1)] != "NA")) 
      AFT <-  k + min(which( CNV.mat4[1,(k+1):ncol(CNV.mat4)] != "NA"))
      BEF_Val <- CNV.mat4[, BEF]
      AFT_Val <- CNV.mat4[, AFT]
      BEF_Val[is.na(BEF_Val)] <- 0
      AFT_Val[is.na(AFT_Val)] <- 0
      #
      CNV.mat4[ , k] <- (AFT_Val + BEF_Val)/2 #(AFT_Val + BEF_Val)/(AFT-BEF)
    }
    
  }
  
  
  ###########################################
  require(robustbase)
  
  CNV.mat51 <- as.matrix(CNV.mat4)
  CNV.mat5 <- matrix(0, ncol=ncol(CNV.mat51), nrow=nrow(CNV.mat51)) # (ceiling(LLL/10)))
  Orig <- 0.5
  LLength <- c(0, cluster.lines[2:(length(cluster.lines)) ])  # c(0,205,733,834,854)
  
  for(j in 1:(length(LLength)-1)){
    NeighborNo <- min(20, floor((LLength[j+1]-LLength[j])/2))
    for(i in (LLength[j]+1):(LLength[j]+floor(NeighborNo))){
      CNV.mat5[i, ] <- CNV.mat51[i, ]*Orig + (1-Orig)*colMedians( CNV.mat51[ setdiff(seq(LLength[j]+1,(i+NeighborNo),1),i),  ])
    }
    if( (LLength[j]+NeighborNo+1) <= (LLength[j+1]-(NeighborNo))  ){
      for(i in (LLength[j]+NeighborNo+1):(LLength[j+1]-(NeighborNo) )){
        CNV.mat5[i, ] <- CNV.mat51[i, ]*Orig + (1-Orig)*colMedians( CNV.mat51[ setdiff(seq(i-NeighborNo,i+NeighborNo,1),i),  ])
      }
    }
    if( (LLength[j+1]-NeighborNo+1) <= LLength[j+1] ){
      for(i in (LLength[j+1]-NeighborNo+1):LLength[j+1]){
        CNV.mat5[i, ] <- CNV.mat51[i, ]*Orig + (1-Orig)*colMedians( CNV.mat51[setdiff(seq(i-NeighborNo,LLength[j+1], 1),i),  ])
      }
    }
    
  }
  
  ########## Defining separating lines
  labels.gloc <- t(as.matrix(c(paste("Chr", 1:22, sep = ""),"ChrX", "ChrY")))
  labels.call.gloc <- rep(NA, ncol(CNV.mat4))
  
  for(i in 1: 22){ 
    labels.call.gloc[break.gloc[i]+10 ] <- as.matrix(labels.gloc[i])
  }
  labels.call.gloc[break.gloc[23] ] <- as.matrix(labels.gloc[23])
  labels.call.gloc[break.gloc[24] ] <- as.matrix(labels.gloc[24])
  
  
  rownames(CNV.mat5) <- ROWlist 
  final.mat <- CNV.mat5
  
  ## Sketching the heatmap 
  
  COL_vec <- c( rep("steelblue",1), rep("white",2),rep("firebrick",1) )
  
  graphics.off()
  plot.new()
  par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
  heatmap.3( final.mat ,  
             main = paste("Heatmap for transcriptions relative to normal control
                          Thr 0.5 of 1, ", NeighborNo ," nearst neighbors", sep="" ), 
             xlab = "Genomic location of expressed genes", 
             ylab= "Cells",
             breaks = seq(-LL, LL, length.out =16), 
             col = colorRampPalette(COL_vec, space = "rgb")(15),  
             Rowv = FALSE, 
             Colv = FALSE, 
             trace ="none", 
             sepwidth = c(0.2,0.2),
             sepcolor = "black",
             scale = "none",  
             labRow = NA,
             labCol = labels.call.gloc,
             cexCol = 1,
             srtCol = 90,
             hieght = 50, 
             width = 400, 
             legend = TRUE, 
             margins = c(4,2),      
             key.xlab = "Transcription level",
             denscol = "grey",
             density.info = "density",
             rowsep = cluster.lines,
             add.expr = abline(v=break.gloc) 
  )
  
  
  
  
}










