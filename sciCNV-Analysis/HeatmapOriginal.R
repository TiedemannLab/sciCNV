

HeatmapOriginal <- function(Mn, 
                            SMPL,
                            comparison,
                            gen.loc = FALSE,   # defaut is FASLE
                            fold_val1=fold_val1,
                            fold_val2=fold_val2,
                            pval1=pval1,
                            pval2=pval2,
                            type = "original"     # default  c("original", "balanced")
){
  
  ########
  if(is.na(type)){
    type <- "original" 
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
  
  ##################

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

RowSep <- length( which( filter_combined & (fold > fold_cutoff))  )

##################### Genomic Locations
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("EnsDb.Hsapiens.v86")

if(gen.loc == TRUE){

orig_row <- rownames(filtered)
symbols <- orig_row
mapIds(org.Hs.eg.db, symbols, "ENSEMBL",  "SYMBOL")
ensids  <- as.matrix( mapIds(org.Hs.eg.db, symbols, "ENSEMBL",  "SYMBOL") )[ ,1  ]

##
cols <- c("SYMBOL",  "MAP")
GenLoc <- select(org.Hs.eg.db, keys=ensids, columns=cols, keytype="ENSEMBL")
GenLoc[,3]
dim(as.matrix(GenLoc[,3]) )

Genomic_loc <- matrix(0,ncol=1, nrow=length(symbols))
for(i in 1:length(symbols)){
  if( symbols[i] %in% GenLoc[,2] ){
    Genomic_loc[i, 1] <- GenLoc[ min(which(GenLoc[,2]== symbols[i] )),3]
  } else 
    Genomic_loc[i, 1] <- c("")
}
dim(Genomic_loc)

rownames(filtered) <- paste(symbols,"   ",Genomic_loc[, 1])

}
####################

graphics.off()
plot.new()
COL_vec <-   c("royalblue1", rep("white",2),"brown1" )

pdf( paste(SMPL,"_DGE_volcanno__",comparison, "pvalue",fold_cutoff1,"-",pvalue_cutoff,"_",type,".pdf", sep="") )

par(mar=c(5,5,4,2)+3,mgp=c(3,1,0))
fontsize_row = 10 - nrow(filtered) / 25

pheatmap(filtered, 
         # main = paste(SMPL,": Subclonality-UMI balanced-log2+0.001", 
         main = paste(SMPL,": -log2+1",
                      comparison, ", p-value < ",pvalue_cutoff,
                      ", fold-cutoff=",fold_cutoff1),
         cluster_cols = FALSE, 
         cluster_rows = FALSE, scale = "row",
         color = colorRampPalette(COL_vec, space = "rgb")(15),
         breaks = seq(-1, 1, length.out = 16),
         show_colnames = FALSE,
         gaps_col= c(ncol(data1)), 
         gaps_row = RowSep,
         fontsize_row=fontsize_row,
         sepcolor = "black"
)

dev.off()

####################

graphics.off()
plot.new()
pdf( paste(SMPL,"_DGE_volcanno__",comparison, "pvalue",fold_cutoff1,"-",pvalue_cutoff,"_",type,".pdf", sep="") )

par(mar=c(5,5,4,2)+3,mgp=c(3,1,0))
fontsize_row = 10 - nrow(filtered) / 25
pheatmap(filtered, 
         # main = paste(SMPL,": Subclonality-UMI balanced-log2+0.001", 
         main = paste(SMPL,": -log2+1",
                      comparison, ", p-value < ",pvalue_cutoff,
                      ", fold-cutoff=",fold_cutoff1),
         cluster_cols = FALSE, 
         cluster_rows = FALSE, scale = "row",
         color = colorRampPalette(COL_vec, space = "rgb")(15),
         breaks = seq(-1, 1, length.out = 16),
         show_colnames = FALSE,
         gaps_col= c(ncol(data1)), 
         gaps_row = RowSep,
         fontsize_row=fontsize_row,
         sepcolor = "black"
)

}
