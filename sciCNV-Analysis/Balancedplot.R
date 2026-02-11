
BalancedPlot <- function(Mn, 
                        SMPL,
                        nUMI,
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
  COL <- c("palevioletred1","steelblue1","black")
  
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
  Mn_NPCs_Begin <- Mn_sorted[, seq(1,L_1Q,1)]
  Mn_NPCs_End <- Mn_sorted[, seq(L_1Q+1, ncol(Mn_sorted),1)]
  
  Mn_Beg_End <- cbind(Mn_NPCs_Begin,Mn_NPCs_End)
  #data <- Mn_Beg_End 
  data <- log2(Mn_Beg_End + 1)
  
  ###############################  STEP 1
  
  Sep_Line <- L_1Q
  data1_org <- data[, 1:Sep_Line]
  data2_org <- data[, (Sep_Line+1):ncol(data)]
  
  #data_org <- data
  #data1_org <- data1
  #data2_org <- data2
  
  
  #### including zero expressions
  
  #data1.mean <- apply(data1, 1, mean)
  #data2.mean <- apply(data2, 1, mean)
  
  
  ############################## STEP 2: Balancing UMIs
 
  
  ## Non-zero expressions
  #nUMI1 <- read.csv( file.choose(), header=TRUE)
  #nUMI1[1,1:10]
  #nUMI <- nUMI1[ , -1]
  
  nUMI_1Q  <- nUMI[1,1:L_1Q]
  nUMI_No1Q <- nUMI[1,(L_1Q+1):length(nUMI)]
  
  #nUMI_1Q <- t(as.numeric(colSums(data1)))
  #colnames(nUMI_1Q ) <- colnames(data1)
  nUMI_1Q_sorted <- as.data.frame(nUMI_1Q)[ , order(nUMI_1Q, decreasing = FALSE), drop = FALSE] 
  #nUMI_1Q_sorted[1,1:10]
  #
  #nUMI_No1Q <- t(as.numeric(colSums(data2)))
  #colnames(nUMI_No1Q ) <- colnames(data2)
  nUMI_No1Q_sorted <- as.data.frame(nUMI_No1Q)[ , order(nUMI_No1Q, decreasing = FALSE), drop = FALSE] 
  #3nUMI_No1Q_sorted[1,1:10]
  
  #plot( as.matrix(nUMI_No1Q_sorted)[1,], pch=16, col=alpha("blue",0.5))
  
  
  
  ############################## STEP 3: Balancing nGenes
  
  data_orig <- log2(Mn_Beg_End+1)
  nGene <- matrix(0, ncol=ncol(data), nrow= 1)
  
  nonzero <- function(x) sum(x != 0)
  for(i in 1:ncol(data)){ 
    nGene[1 ,i] <- nonzero(data[,i])
  }
  colnames(nGene) <- colnames(Mn_Beg_End)
  
  nGene_1Q <- t(as.matrix(nGene[ , seq(1,L_1Q ,1)] ))
  colnames(nGene_1Q ) <- colnames(nGene)[seq(1,L_1Q ,1)]
  nGene_No1Q <- t(as.matrix(nGene[ , (L_1Q+1):(L_1Q+L_None)]))
  colnames(nGene_No1Q ) <- colnames(nGene)[(L_1Q+1):(L_1Q+L_None)]
  
  nGene_1Q_sorted <- as.data.frame(nGene_1Q)[ , order(nUMI_1Q, decreasing=FALSE), drop = FALSE] 
  nGene_No1Q_sorted <- as.data.frame(nGene_No1Q)[ , order(nUMI_No1Q, decreasing=FALSE), drop = FALSE]
  
  
  
  sample=c(colnames(as.matrix(nUMI_1Q_sorted)), colnames(as.matrix(nUMI_No1Q_sorted)))
  specie=c( rep(  SbCl1 , ncol(as.matrix(nUMI_1Q_sorted) )) ,  rep(SbCl2 , ncol(as.matrix(nUMI_No1Q_sorted) )) )
  gene1=c( as.matrix(nUMI_1Q_sorted) , as.matrix(nUMI_No1Q_sorted) )
  gene2=c( as.matrix(nGene_1Q_sorted) , as.matrix(nGene_No1Q_sorted)  )
  data=data.frame(sample,specie,gene1,gene2)
  
  #####################
  graphics.off()
  plot.new()
  pdf( paste(SMPL,"_",comparison,"__nUMI_nGene.pdf",sep=""), width=20, height=15, paper = "a4r")
  
  print(xyplot(gene1 ~ gene2 | specie , data=data , pch=20 , cex=3 , col=rgb(0.2,0.4,0.8,0.5), 
         cex.lim=2, cex.lab=2, cex.main=2, xlab=list(label="nUMI", cex=2), 
         main=list(label=paste(SMPL,": nUMI vs. nGenes"), cex=2, col="brown"),
         ylab=list(label="nGene", cex=2), 
         scales=list(tck=c(1,1), x=list(cex=2), y=list(cex=2)))  )
  dev.off()
  
  
  
  ######### STEP 4: Corresponding celles with same UMI and nGenes
  
  ### FIRSTLY UMIs
  
  UMI1 <- as.matrix(nUMI_1Q_sorted)
  UMI2 <- as.matrix(nUMI_No1Q_sorted)
  
  UMI1_mdf <- rep(0, max(ncol(UMI1),ncol(UMI2)))    #UMI1 Modified
  UMI2_mdf <- rep(0, max(ncol(UMI1),ncol(UMI2)))     #UMI2 Modified
  
  ##
  UMI1_mdf[1] <- UMI1[1,1]
  UMI2_mdf[1] <- UMI2[1,1]
  UMI2_mdf[2] <- UMI2[1, min(which(UMI2 > max(UMI1_mdf[1],UMI2_mdf[1]) ))]
  UMI1_mdf[2] <- UMI1[1, min(which(UMI1 > UMI2_mdf[1] ))]
  ##
  for(i in 3:max(ncol(UMI1),ncol(UMI2))){
    UMI2_mdf[i] <- UMI2[1, min(which(UMI2 > UMI1_mdf[i-1] ))]
    UMI1_mdf[i] <- UMI1[1, min(which(UMI1 > UMI2_mdf[i] ))]
  }
  
  ##
  UMI1_mdf <- t(as.matrix( UMI1_mdf[!is.na(UMI1_mdf)] ))
  UMI2_mdf <- t(as.matrix( UMI2_mdf[!is.na(UMI2_mdf)] ))
  
  L_UMI <-  min( ncol(UMI1_mdf), ncol(UMI2_mdf) )
  ##
  graphics.off()
  plot.new()
  
  pdf( paste(SMPL,"_",comparison,"__nUMI_modified.pdf",sep=""), width=20, height=15, paper = "a4r")
  
  par(mfrow=c(1,2))
  par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
  par(bty="l")
  plot(UMI1_mdf[1,], pch=21, col=COL[1], 
       bg= alpha(COL[1], 0.3), 
       cex=3, 
       main=paste(SMPL,": ",comparison, sep=""), cex.axis=2,
       cex.lab=2,cex.lim=2,cex.main=2,col.main="brown", 
       ylab="nUMI (modified)", xlab="Cells")
  par(new=TRUE); points( UMI2_mdf[1,], col=COL[2], 
                         bg= alpha(COL[2], 0.3), 
                         pch=24, axis=FALSE, cex=2.2)
  
  legend("top", legend=c(paste("Subclone ",SbCl1),paste("Subclone ",SbCl2)),
         col=alpha(c(COL[1],COL[2]),0.7), pch=16, cex=1.5 )
  
  
  plot(UMI2_mdf[1,1:L_UMI]~UMI1_mdf[1,1:L_UMI], pch=21, 
       col= COL[3], bg=alpha(COL[3],0.3), cex=3,
       main=paste(SMPL,": ",comparison), cex.axis=2,
       cex.lab=2,cex.lim=2,cex.main=2,col.main="brown", 
       ylab=paste("nUMI ",SbCl1," (modified)",sep=""), 
       xlab=paste("nUMI ",SbCl2," (modified)",sep=""))
  abline(a=0,b=1, col="brown1", lwd=2)
  
  dev.off()
  
  
  SubCell1 <- colnames(nUMI_1Q_sorted[1,which(nUMI_1Q_sorted %in% UMI1_mdf  ) ])
  SubCell2 <- colnames(nUMI_No1Q_sorted[1,which(nUMI_No1Q_sorted %in% UMI2_mdf  ) ])
  
  
  #UMI2_mdf1 <- t(as.matrix( UMI2_mdf[1, -ncol(UMI2_mdf)])
  #colnames(UMI1_mdf) <- SubCell1
  #colnames(UMI2_mdf1) <- SubCell2
  
  ####****####****####****####****####****####
  ### SECONDLY nGENEs
  
  nGene11 <- as.matrix( sort( nGene_1Q_sorted[1,SubCell1 ]) )
  colnames(nGene11) <- colnames( sort(nGene_1Q_sorted[1,SubCell1 ]) )
  #
  nGene21 <- as.matrix( sort( nGene_No1Q_sorted[1,SubCell2 ]) )
  colnames(nGene21) <- colnames( sort(nGene_No1Q_sorted[1,SubCell2 ]) )
  
  nGene1_mdf <- rep(0, max(ncol(nGene11),ncol(nGene21)))    #UMI1 Modified
  nGene2_mdf <- rep(0, max(ncol(nGene11),ncol(nGene21)))     #UMI2 Modified
  
  
  ##
  nGene1_mdf[1] <- nGene11[1,1]
  nGene2_mdf[1] <- nGene21[1,1]
  nGene2_mdf[2] <- nGene21[1, min(which(nGene21 > max(nGene1_mdf[1],nGene2_mdf[1]) ))]
  nGene1_mdf[2] <- nGene11[1, min(which(nGene11 > nGene2_mdf[1] ))]
  ##
  for(i in 3:max(ncol(nGene11),ncol(nGene21))){
    nGene2_mdf[i] <- nGene21[1, min(which(nGene21 > nGene1_mdf[i-1] ))]
    nGene1_mdf[i] <- nGene11[1, min(which(nGene11 > nGene2_mdf[i] ))]
  }
  
  ##
  nGene1_mdf <- t(as.matrix( nGene1_mdf[!is.na(nGene1_mdf)]))
  nGene2_mdf <- t(as.matrix( nGene2_mdf[!is.na(nGene2_mdf)]))
  dim(nGene1_mdf);dim(nGene2_mdf)
  Length <- min(ncol(nGene1_mdf), ncol(nGene2_mdf))
  ##
  
  graphics.off()
  plot.new()
  
  pdf( paste(SMPL,"_",comparison,"__nGene_modified_after_nUMI.pdf",sep=""), width=20, height=15, paper = "a4r")
  
  par(mfrow=c(1,2))
  par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
  par(bty="l")
  
  
  plot(nGene1_mdf[1,], pch=21, col=COL[1], 
       bg= alpha(COL[1], 0.3), 
       cex=3, 
       main=paste(SMPL,": ",comparison, sep=""), cex.axis=2,
       cex.lab=2,cex.lim=2,cex.main=2,col.main="brown", 
       ylab="nGene (modified)", xlab="Cells" )
  par(new=TRUE); points( nGene2_mdf[1,], col=COL[2], 
                         bg= alpha(COL[2], 0.3), 
                         pch=24, axis=FALSE, cex=2.2)
  
  legend("top", legend=c(paste("Subclone ",SbCl1),paste("Subclone ",SbCl2)),
                         col=c(COL[1],COL[2]), pch=16, cex=1.5 )
  
  
  plot(nGene2_mdf[1,1:Length]~nGene1_mdf[1,1:Length], pch=21, 
       col= alpha(COL[3],0.5), bg=alpha(COL[3],0.1), cex=3,
       main=paste(SMPL,": ",comparison), cex.axis=2,
       cex.lab=2,cex.lim=2,cex.main=2,col.main="brown", 
       ylab=paste("nGene ",SbCl1," (modified)",sep=""), 
       xlab=paste("nGene ",SbCl2," (modified)",sep=""))
  abline(a=0,b=1, col="brown1", lwd=2)
  
  dev.off()
  ##
  SubCell11 <- colnames( t(as.matrix(nGene11[1,which(nGene11 %in% nGene1_mdf  ) ])))
  SubCell21 <- colnames(  t(as.matrix(nGene21[1,which(nGene21 %in% nGene2_mdf  ) ])))
  
  

  ######### STEP 5: Using the selcted cells only
 
  data1 <- data1_org[, SubCell1 ]
  data2 <- data2_org[, SubCell2 ] 
  data <- cbind(data1, data2)
  
  
  data1.mean <- apply(data1, 1, mean)
  data2.mean <- apply(data2, 1, mean)
  
  ################  Distributions on NOT median normalziation
  
  ### Original
  graphics.off()
  plot.new()
  
  
  pdf( paste(SMPL,"_",comparison,"__nUMI_nGene_distribution.pdf",sep=""), width=20, height=15, paper = "a4r")
  
  par(mfrow=c(2,2))
  par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
  par(bty="l")
  plot( density(as.matrix(nUMI_1Q_sorted)[1,]), main=paste(SMPL,"-Original"), type="l", 
        lty=1, lwd=2, 
        ylim=c(min(density(as.matrix(nUMI_No1Q_sorted)[1,])$y),1.1*max(density(as.matrix(nUMI_No1Q_sorted)[1,])$y)), 
        #xlim=c(-0.1,0.1), 
        col="red", cex.axis=2,
        col.main="brown", cex.lim=2, cex.lab=2, cex.main=2,
        xlab="nUMI", ylab="Density of UMI"
  ) 
  par(new=TRUE); points(density(as.matrix(nUMI_No1Q_sorted)[1,]), type="l", lty=1, lwd=2,col="blue") 
  
  ### UMI Balanced
  par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
  par(bty="l")
  plot( density(UMI1_mdf[1,]), main=paste(SMPL,"-UMI balanced"), type="l", 
        lty=1, lwd=2, cex.axis=2,
        #ylim=c(min(density(as.matrix(nUMI_No1Q_sorted)[1,])$y),1.2*max(density(as.matrix(nUMI_No1Q_sorted)[1,])$y)), 
        #xlim=c(-0.1,0.1), 
        col="red",
        col.main="brown", cex.lim=2, cex.lab=2, cex.main=2,
        xlab="nUMI", ylab="Density of UMI"
  ) 
  par(new=TRUE); points(density(UMI2_mdf[1,]), type="l", lty=1, lwd=2,col="blue") 
  
  
  ### UMI and nGene Balanced
  par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
  par(bty="l")
  plot( density(nGene1_mdf[1,]), main=paste(SMPL,"-UMI&nGene balanced"), type="l", 
        lty=1, lwd=2, cex.axis=2,
        #ylim=c(min(density(as.matrix(nUMI_1Q_sorted)[1,])$y),max(density(as.matrix(nUMI_1Q_sorted)[1,])$y)), 
        #xlim=c(-0.1,0.1), 
        col="red",
        col.main="brown", cex.lim=2, cex.lab=2, cex.main=2,
        xlab="nUMI", ylab="Density of UMI"
  ) 
  par(new=TRUE); points(density(nGene2_mdf[1,]), type="l", lty=1, lwd=2,col="blue") 
  
  
  legend(3000, 8e-04, #2*max(density(as.matrix(nUMI_1Q_sorted)[1,])$y),
         legend=c( SbCl1, SbCl2 #, "Balancing UMIs & nGenes"
         ),
         col=c("red","blue"), fill= c("red","blue"), cex=1.5 )
  
  dev.off()
  ########################
  
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
  
  pdf( paste(SMPL,"_",comparison,"__volcano_modified_data.pdf",sep=""), width=20, height=15, paper = "a4r")
  
  par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
  par(bty="l")
  plot(fold, -log10(pvalue), 
       main = paste(SMPL,": UMI-balanced ", comparison),
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
          pch = 21, col = COL[1], bg= alpha(COL[1],0.3) , cex=3)
  
  
  
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
         col=c(COL[1],COL[2]), pch=16,box.lwd = 0, cex=1.5 )
  
  dev.off()
  
  
  ##########
  graphics.off()
  plot.new()
  par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
  par(bty="l")
  plot(fold, -log10(pvalue), 
       #main = paste(SMPL,": UMI-balanced ", comparison),
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
          pch = 21, col = COL[1], bg= alpha(COL[1],1) , cex=3)
  
  
  return(filtered)
  
 
}
  
