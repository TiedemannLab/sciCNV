
VolcanoPlot <- function(Mn, 
                         SMPL,
                         comparison,
                         gene.norm=FALSE,   # default is FALSE
                         fold_val1,
                         fold_val2,
                         pval1,
                         pval2
                         ){
  
  ##############
  
  if( is.na(gene.norm) ){
    gene.norm = TRUE
  }
  
  ##############
  COL <- c("palevioletred1","steelblue1")
  
  
  CLS1Q <- sapply(strsplit(colnames(Mn)[5], split='_', fixed=TRUE), function(x) (x[1]))
  CLSno1Q <- sapply(strsplit(colnames(Mn)[ncol(Mn)], split='_', fixed=TRUE), function(x) (x[1]))
  CLS1Q; CLSno1Q
  
  L_1Q <- length(grep(CLS1Q, colnames(Mn)))
  L_None <- length(grep(CLSno1Q, colnames(Mn)))
  L_1Q; L_None
  
  SbCl1 <- sapply(strsplit(comparison, split=' ', fixed=TRUE), function(x) (x[1]))
  SbCl2 <- sapply(strsplit(comparison, split=' ', fixed=TRUE), function(x) (x[3]))
  SbCl1;SbCl2 
  
  
  
  ###############
  
  if(gene.norm == TRUE){
  
  MMed <- matrix(0 , nrow = nrow(Mn), ncol = 1)
  for( i in 1:nrow(Mn)){
    if( sum( Mn[i, ])==0 ){
      MMed[i, 1] <- 0
    } else {
      # MMed[i, 1] <- median( Mn[i, ][ Mn[i, ] >0 ])
      MMed[i, 1] <- mean( Mn[i, ][Mn[i, ]>0] )
    }
  }
  
  # Max_MMed <- mean(MMed[ , 1][MMed[ , 1]> 0 ])
  
  Mn2 <- matrix(0 , nrow = nrow(Mn), ncol = ncol(Mn))
  for( i in 1:nrow(Mn)){
    if( MMed[i, 1] > 0 ){
      Mn2[ i, ] <- Mn[i, ]/MMed[i, 1]
    } else {
      Mn2[ i, ] <- 0
    }
    
  }
  
  Mn3 <- as.matrix(Mn2)
  colnames(Mn3) <- colnames(Mn)
  rownames(Mn3) <- rownames(Mn)
  Mn <- Mn3
  }
  
  ######################################
  
  Mn_sorted <- Mn
  Mn_subclones_1q <- Mn_sorted[, seq(1,L_1Q,1)]
  Mn_subclones_no1q <- Mn_sorted[, seq(L_1Q+1, ncol(Mn_sorted),1)]
  Mn_NPCs_Begin  <- Mn_subclones_1q
  Mn_NPCs_End  <- Mn_subclones_no1q 
  
  Mn_Beg_End <- cbind(Mn_NPCs_Begin,Mn_NPCs_End)
  #data <- Mn_Beg_End 
  #data_log2 <- log2(Mn_Beg_End + 1)
  data <- log2(Mn_Beg_End + 1)
  
  ##################
  
  Sep_Line <- L_1Q
  data1 <- data[, 1:Sep_Line]
  data2 <- data[, (Sep_Line+1):ncol(data)]
  data1.mean <- apply(data1, 1, mean)
  data2.mean <- apply(data2, 1, mean)
  
  fold = data1.mean - data2.mean
  
  ###################
  pvalue <- NULL
  tstat <- NULL
  
  for(i in 1:nrow(data)){
    x <- data1[i, ]
    y <- data2[i, ]
    t <- t.test(x,y)
    pvalue[i] <- t$p.value
    tstat[i] <- t$statistic
  }
  pvalue <- as.matrix(pvalue); fold <- as.matrix(fold);
  rownames(pvalue) <- rownames(data)
  rownames(fold) <- rownames(data)
  
 #################### Sketching

  fold_cutoff1 <- fold_val1 
  fold_cutoff2 <- -fold_val2 
  pvalue_cutoff1 <- pval1  
  pvalue_cutoff2 <- pval2 
  pvalue_cutoff <- min(pvalue_cutoff1, pvalue_cutoff2)
  fold_cutoff <- min(fold_cutoff1, abs(fold_cutoff2))
  
  filter_by_fold <-  (( fold >= fold_cutoff1  ) | ( fold <= fold_cutoff2  ))
  dim(data[filter_by_fold, ])
  
  filter_by_pvalue = pvalue <= pvalue_cutoff
  filter_combined = filter_by_fold & filter_by_pvalue
  filter_combined =  (( fold >= fold_cutoff1  ) & ( pvalue <= pvalue_cutoff1  )) | (( fold <= fold_cutoff2  ) & ( pvalue <= pvalue_cutoff2  ))
  filtered <- data[filter_combined,]
  
  filtered_ranked <- filtered[order(fold[rownames(filtered),1 ], decreasing=TRUE) , ]
  filtered <- as.matrix(filtered_ranked)
  
  ##########
  graphics.off()
  plot.new()
  
  pdf( paste(SMPL,"_",comparison,"__volcano_original_data.pdf",sep=""), width=20, height=15, paper = "a4r")
  
  par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
  par(bty="l")
  plot(fold, -log10(pvalue), 
       main = paste(SMPL,": subclone ", comparison),
       cex.lim=2,cex.lab=2, cex.axis=2, cex.main=2, cex=1.5, 
       pch=16, col.main='brown', col=alpha("black", 0.1)
       #xlim=c(1.1*min( fold ),1.1*max( fold ))
       #xlim=c(1.7,4.4), # xlim=c(1.7,2),
       #, ylim=c(0.01,50)
       # ,ylim=c(0, 1.2*max(- log10(pvalue)))
       ,xlab="log2(fold change+1)",ylab="-log10(p-value)"
       )
  points (fold[filter_combined & fold < 0],
          -log10(pvalue[filter_combined & fold < 0]),
          pch = 21, col = COL[2], bg= alpha(COL[2],0.3) ,cex=3)
  points (fold[filter_combined & fold > 0],
          -log10(pvalue[filter_combined & fold > 0]),
          pch = 21, col = COL[1], bg= alpha(COL[1],0.3), cex=3)
  
  #f39f4e
  #d32414
  
  text(fold[filter_combined & fold > fold_cutoff ],
       -log10(pvalue[filter_combined & fold > fold_cutoff]), 
       labels=rownames(data)[which(filter_combined & fold > 0)],
       pch = 16, col = "black", cex=0.7)
  
  text(fold[filter_combined & fold < -fold_cutoff ],
       -log10(pvalue[filter_combined & fold <  -fold_cutoff]), 
       labels=rownames(data)[which(filter_combined & fold <  -fold_cutoff)],
       pch = 16, col = "black", cex=0.7)
  legend("top", legend=c(paste("Fold cutoff=", fold_cutoff1,",",fold_cutoff2),
                             paste( "p-value=",pvalue_cutoff1,",",pvalue_cutoff2)), cex=1 )
  legend("bottomleft", legend=c(paste("Subclone ",SbCl1),paste("Subclone ",SbCl2)),
         col=c(alpha(COL[1],0.7), alpha(COL[2],0.7)), 
         pch=16,box.lwd = 0, cex=1.5 )
  
  dev.off()
  
  ##############
  
  graphics.off()
  plot.new()
  par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
  par(bty="l")
  plot(fold, -log10(pvalue), 
       #main = paste(SMPL,": subclone ", comparison),
       cex.lim=2,cex.lab=2, cex.axis=2, cex.main=2, cex=1.5, 
       pch=16, col.main='brown', col=alpha("black", 0.1)
       #xlim=c(1.1*min( fold ),1.1*max( fold ))
       #xlim=c(1.7,4.4), # xlim=c(1.7,2),
       #, ylim=c(0.01,50)
       # ,ylim=c(0, 1.2*max(- log10(pvalue)))
       ,xlab="log2(fold change+1)",ylab="-log10(p-value)"
  )
  points (fold[filter_combined & fold < 0],
          -log10(pvalue[filter_combined & fold < 0]),
          pch = 21, col = COL[2], bg= alpha(COL[2],1) ,cex=3)
  points (fold[filter_combined & fold > 0],
          -log10(pvalue[filter_combined & fold > 0]),
          pch = 21, col = COL[1], bg= alpha(COL[1],1), cex=3)
  

  return(filtered)
   
  
}
