

#######################################################################################
######                          CNV_htmp_glist function                         ####### 
######                   Heatmap of CNV-curves: against gene list (gnlist)      #######
######  Tiedemann Lab - Princess Margaret Cancer centre, University of Toronto  #######
######                       copyright@AliMahdipourShirayeh                     #######
#######################################################################################

# CNV.mat2: copy number variation matrix 
# Gen.Loc: Genomic Location matrix with list of genes, their associated chromosome numbers, starts and ends
# clustering: to cluster or not the cells based on their CNV similarities. Default is "FALSE"
# clustering.type: the selected method will be considered to generate heatmap. Possible options are "pearson", 
#        "euclidean", " spearman", ... "original" stands for keeping the same clusters as original without using any unsupervised clustering.
#         Defualt is "pearson" clustering method. Only when supervised = "TRUE"
# sorting: TRUE for sorting cells based on their tumor score from the largest to smallest tumor scores. Default is FASLE.
# CNVscore: is the total score matrix of all cells (possibly ranked separately for each cluster). Can be only used when sorting.clusters = TRUE.
# cluster.lines: is a list of values which separate different clusters accross entire population, onlys used if there
#        exist several clusters other than just test and control populations
# break.glist: is a set of values each defines a vertical line that separates chromosomes
# No.test: number of test cells included in the data, potentially is used for separating diverse populations of 
#       of test annd control cells in the heatmap

CNV_htmp_glist <- function(CNV.mat2,
                           Gen.Loc,
                           clustering = FALSE,            # TRUE or FALSE option 
                           clustering.type = c("pearson", "kendall", "spearman"),   # "pearson", "kendalln", " spearman", defualt: "pearson"
                           sorting = FALSE,               # TRUE or FALSE
                           CNVscore = NULL,               # Only exists when sorting = TRUE 
                           cluster.lines = NULL, 
                           break.glist = break.glist,     # Needs to be defined as separation lines for chromosomes       
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
  if ( Reduce("|", is.null(break.glist)) ){
    break.glist <- c(0, ncol(CNV.mat2))
  }
  
  ##### sorting cells within each cluster based on CNV-scores from the largest to the smallest (if applicable)
  
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
    CNV.mat1 <- rbind(as.matrix(CNV.mat2[(No.test+1):nrow(CNV.mat2), ]) , as.matrix(CNV.mat2[1:No.test, ]) )
    rownames(CNV.mat1) <-  c(rownames(as.matrix(CNV.mat2[(No.test+1):nrow(CNV.mat2), ])), rownames(as.matrix(CNV.mat2[1:No.test, ]))  )   
    colnames(CNV.mat1) <- colnames(CNV.mat2)
  }
  
  
  ROWlist <- rownames(CNV.mat1)
  COLlist <- colnames(CNV.mat1)
  
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
  CNV.mat1 <- as.matrix(CNV.mat11)
  LL <- 1
  TT <- 0.5
  CNV.mat1[which(CNV.mat1 >  LL)] <-  LL
  CNV.mat1[which( (CNV.mat1 <  TT) & (CNV.mat1 >  -TT) )] <- 0.0
  CNV.mat1[which(CNV.mat1 <  -LL)] <-  -LL
  
  CNV.mat <- as.matrix(CNV.mat1)
  rownames(CNV.mat) <- matrix(NA, ncol=1, nrow=nrow(CNV.mat1))
  
  ##############################################
  ## Heatmap of CNV-curves against list of genes
  ##############################################
  
  ## Using k-nearest neighbors (kNN)-method
  CNV.mat31 <- as.matrix(CNV.mat)
  CNV.mat3 <- matrix(0, ncol=ncol(CNV.mat31), nrow=nrow(CNV.mat31)) 
  Orig <- 0.5
  LLength <- c(0, cluster.lines[ 2:( length(cluster.lines) ) ]) 
  
  for(j in 1:( length(LLength) - 1 )){
    NeighborNo <- min(20, floor((LLength[j+1]-LLength[j])/2))
    
    for(i in (LLength[j]+1):(LLength[j]+floor(NeighborNo))){
      CNV.mat3[i, ] <- CNV.mat31[i, ]*Orig + (1-Orig)*colMedians(CNV.mat31[ setdiff(seq(LLength[j]+1,(i+NeighborNo),1),i),  ])
    }
    if( (LLength[j]+NeighborNo+1) <= (LLength[j+1]-(NeighborNo) )){
      for(i in (LLength[j]+NeighborNo+1):(LLength[j+1]-(NeighborNo) )){
        CNV.mat3[i, ] <- CNV.mat31[i, ]*Orig + (1-Orig)*colMedians( CNV.mat31[ setdiff(seq(i-NeighborNo,i+NeighborNo,1),i),  ])
      }
    }
    for(i in (LLength[j+1]-NeighborNo+1):LLength[j+1]){
      CNV.mat3[i, ] <- CNV.mat31[i, ]*Orig + (1-Orig)*colMedians(CNV.mat31[setdiff(seq(i-NeighborNo,LLength[j+1], 1),i),  ])
    }
    
  }
  
  CNV.mat3 <- as.matrix(CNV.mat3)
  
  ## Assignning location of genes on the genome
  
  label.glist <- t(as.matrix(c(paste("Chr", 1:22, sep = ""),"ChrX", "ChrY")))
  label.call.glist <- rep(NA, ncol(CNV.mat2))
  
  for(i in 1: 22){ 
    if( break.glist[i] <  ncol(CNV.mat2) - 10){
      label.call.glist[break.glist[i]+10 ] <- as.matrix(label.glist[i])
    } 
  }
  
  if( break.glist[23] <  ncol(CNV.mat2)  ){
    label.call.glist[break.glist[23] ] <- as.matrix(label.glist[23])
    label.call.glist[break.glist[24] ] <- as.matrix(label.glist[24])
  } else {
    label.call.glist[break.glist[23] ] <- as.matrix(label.glist[23])
  }
  
  
  rownames(CNV.mat3) <- ROWlist 
  final.mat <- CNV.mat3
  
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
             #treeheight_row = 0.2,
             sepwidth = c(0.2,0.2),
             sepcolor = "black",
             scale = "none",  
             labRow = NA,
             labCol = label.call.glist,
             cexCol = 1,
             srtCol = 90,
             #cutree_rows = 2, 
             hieght = 50, 
             width = 400, 
             legend = TRUE, 
             margins = c(4,2),      
             key.xlab = "Transcription level",
             denscol = "grey",
             density.info = "density",
             rowsep = cluster.lines ,
             add.expr = abline(v = break.glist) 
  )
  
  
  
}










