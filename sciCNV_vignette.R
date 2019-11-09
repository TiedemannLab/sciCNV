
#######################################################################################
######                             Vignette of CNV pipeline                     ####### 
###### ((  RTAM1/2- normalization, sciCNV analysis, CNV-subclonal analysis  ))  #######
######  Tiedemann Lab - Princess Margaret Cancer centre, University of Toronto  #######
######                       copyright@AliMahdipourShirayeh                     #######
#######################################################################################

.rs.restartR()
options(stringsAsFactors = FALSE)

library(base)
library(robustbase)
library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)
library(qlcMatrix)
library(devtools)
library(svd)
library(ggplot2)
library(ggridges)
library(viridis)
require(scales)
library(RColorBrewer)
library(Rtsne)
library(reticulate)
library(umap)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

##
path.code <- "./sciCNV-Analysis/"
source(file.path(path.code, "Mito_umi_gn.R"))
source(file.path(path.code, "RTAM_normalization.R"))
source(file.path(path.code, "sciCNV.R"))
source(file.path(path.code, "Scaling_CNV.R"))
source(file.path(path.code, "CNV_score.R"))
source(file.path(path.code, "sciCNV.R"))
source(file.path(path.code, "Sketch_AveCNV.R"))
source(file.path(path.code, "CNV_htmp_glist.R"))
source(file.path(path.code, "CNV_htmp_gloc.R"))
source(file.path(path.code, "Opt_MeanSD_RTAM1.R"))
source(file.path(path.code, "Opt_MeanSD_RTAM2.R"))
source(file.path(path.code, "heatmap_break_glist.R"))
source(file.path(path.code, "heatmap_break_gloc.R"))

## Reading raw data with a list of genes on the first column

raw.data1 <- read.table("./Dataset/Sample_100_CPCs__with__100_NPCs.txt", sep = '\t',header = TRUE)  
raw.data2 <- as.matrix(raw.data1[ , -1])
rownames(raw.data2) <- raw.data1[ , 1]
colnames(raw.data2) <- colnames(raw.data1[ , -1])

W <- ncol(raw.data2)
No.test <- 100                     # Number of test cells
No.control <- ncol(raw.data2)-100  # Number of control cells
Col_Sum <- t(as.numeric(colSums(raw.data2)))

##################################################
# Quality Control (QC): Eliminating damahged cells 
##################################################

## Reading UMI and mitochondrial gene expressions associated to the data
nUMI <- t(as.numeric(colSums(raw.data2)))
colnames(nUMI) <- colnames(raw.data2)

mito.genes <-  read.table("./Dataset/Sample_100_CPCs_Mitochondrial.txt", sep = '\t',header = TRUE)
mito.genes <- as.matrix(mito.genes[,-1])
percent.mito.G <- t(as.matrix(colSums(mito.genes)))/ ( Col_Sum[1:No.test] + colSums(mito.genes))

nGene1 <- matrix(0, ncol = ncol(raw.data2) , nrow = 1)
nonzero <- function(x){ sum(x != 0) }

nGene1 <- lapply( 1:ncol(raw.data2), function(i){ nonzero(raw.data2[, i])} ) 
nGene1 <- t(as.numeric(nGene1))
colnames(nGene1) <- colnames(raw.data2)

#--------------------------------------------
# Damaged Cells Removal in Entire Population 
#--------------------------------------------
MMS <- CreateSeuratObject(counts = raw.data2, project = "Sample1")
damaged_cells <- Mito_umi_gn(mat = MMS, 
                             percent.mito.G = percent.mito.G,
                             nUMI = nUMI,
                             nGene = nGene1,
                             No.test = No.test,
                             drop.mads = 1)

#----------------------------
##  Excluding Damaged Cells
#----------------------------
if( length(damaged_cells) > 0 ){
  Outliers <- damaged_cells
  raw.data <- raw.data2[, - Outliers]
  nUMI <- t(as.matrix(nUMI))[- Outliers]
} else{
  raw.data <- raw.data2
}

nUMI <- t(as.matrix(nUMI))
raw.data <- as.matrix(raw.data )
rownames(raw.data) <- raw.data1[ , 1]
colnames(nUMI) <- colnames(raw.data)

#--------------------------------------------------
##  Sorting based on UMI (from largest to smallest)
#--------------------------------------------------
raw.data <- raw.data2[, c(colnames(sort(as.data.frame(nUMI)[1:No.test], decreasing=TRUE))
                          ,colnames(sort(as.data.frame(nUMI)[(No.test+1):ncol(raw.data)], decreasing=TRUE))
                          ), drop=FALSE]
rownames(raw.data) <- rownames(raw.data2) 

#----------------------------
##  RTAM1/2 Normalziation
#----------------------------
norm.data <- RTAM_normalization(mat = raw.data,            
                                method = "RTAM2",      
                                Min_nGn  = 250,       
                                Optimizing = FALSE)
rownames(norm.data) <- rownames(raw.data2) 
colnames(norm.data) <- colnames(raw.data)

## Sketching non-zero expressions

graphics.off()
plot.new()
par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
par(xaxs="i", yaxs="i") 
par(bty="l")
plot(matrix(1,ncol=1,nrow=nrow(as.matrix(norm.data[,1][norm.data[,1]>0]))), 
     log2(as.matrix(norm.data[,1][norm.data[,1]>0]) +1  ), pch=16, cex=0.3, 
     col="darkgray" ,
     xlim=c(-ncol(norm.data)*.01, ncol(norm.data)*1.01),
     ylim=c(0,  4.2), xlab = "Cells (ranked by UMI)",
     ylab = expression("Expressions ("*Log[2]*"(.+1 ))"), 
     cex.lab = 2, cex.axis = 2, cex.main=2)

for(i in 2:ncol(norm.data)){
  par(new=TRUE)
  points(  matrix(i,ncol=1,nrow=nrow(as.matrix(norm.data[,i][norm.data[,i]>0]))), 
           log2(as.matrix(norm.data[,i][norm.data[,i]>0]) +1  ), pch=16, cex=0.3, 
           col="darkgray")
}
title( paste("Sample1, RTAM2-normalization, cutoff ", 250," nGene ",250,sep=""), 
       col.main = "brown", cex.main = 2)

#---- 95% commonly expressed genes
Sqnce <- seq(1,ncol(norm.data),1)
Common.mat <-   as.matrix(norm.data[which(rowSums(norm.data[,Sqnce ] != 0) > 0.95*ncol(norm.data) ), ] )

COLMED <- log2(colMeans(as.matrix(Common.mat[Common.mat>0])) +1)
COLMED <- t(as.matrix(COLMED) )
for(j in 1:ncol(norm.data)){
  par( new=TRUE)
  points(matrix(j,ncol=1,nrow=1),
         log2(mean(as.matrix(Common.mat[,j][Common.mat[,j]>0])) + 1 ) ,
         axis = FALSE , col="red" , pch=15, cex =0.5      
  )
} 
legend(0,0.75,bty="n",pch=16,col=c("red",NA), cex=1.5, legend=paste("Mean of 95% commonly expressed genes"))

#---- Average expression of housekeeping genes
Houskeeping_gene_list <- read.table( "./Dataset/HouseKeepingGenes.txt", sep = '\t',header = TRUE)
HK_mat  <- norm.data[which(raw.data1[ , 1]%in% t(as.matrix(Houskeeping_gene_list))), ]
HK_mat <- as.matrix(HK_mat)
colnames(HK_mat) <- colnames(norm.data)

Mean_HK_mat <- matrix(0, ncol = ncol(norm.data) , nrow = 1)
for(k in 1: (ncol(norm.data))){
  Mean_HK_mat[1,k] <- as.numeric(mean(HK_mat[,k][HK_mat[,k]>0]) )
}
par( new=TRUE)
points(log2(Mean_HK_mat[1,] +1 ), col="blue" , pch=15, cex =0.5 )
legend(0,0.5,bty="n",pch=16,col=c("blue",NA), cex=1.5, legend=paste("Mean of Houskeeping gene expressions"))

#---------------------------------
##  Excluding non-expressed genes
#---------------------------------
train1 <- as.matrix(norm.data)
colnames(train1) <- colnames(norm.data)
train <- as.matrix(train1[which(rowSums(train1 != 0) >= 1  ), ] )

#---------------------------------
##  clustering the normalized data
#---------------------------------
## Finding matrix of single cells (MSC) as a seurat object for clustering
MSC <- CreateSeuratObject(counts = train, project = "sample1")
MSC <- FindVariableFeatures(object = MSC)
all.genes <- rownames(MSC)
MSC <- ScaleData(MSC, features = all.genes)
MSC <- RunPCA(object = MSC, do.print = TRUE, pcs.print = 1:10, genes.print = 10)
MSC <- FindNeighbors(object = MSC)
MSC <- FindClusters(object = MSC)
MSC <- RunTSNE(object = MSC, reduction.type = "pca",   #"cca.aligned",
               dims.use = 1:5, do.fast = TRUE)
DimPlot(object = MSC,reduction = "tsne", pt.size = 3, label = TRUE, label.size = 4)

## Running UMAP
MSC <- RunUMAP(MSC, dims = 1:10)
DimPlot(MSC, reduction = "umap", pt.size = 3, label = TRUE, label.size = 4)

#---------------------------------------------------------
##  Running sciCNV function to derive CNV-curves per cell
#---------------------------------------------------------
tst.index  <- seq(1, No.test , 1)                      # No of test cells
ctrl.index <- seq(No.test+1, ncol(norm.data), 1)      # No of controcl cells

## generating infered-CNV data for (test and/or control) cells 
CNV.data <- sciCNV(norm.mat = norm.data, 
                   No.test = No.test, 
                   sharpness  = 1, 
                   baseline_adj  = FALSE,  
                   baseline = 0)

#--------------------------------------------------------------
##  Scaling and Filtering noise of the iCNV curves
#--------------------------------------------------------------

## Scaling CNV-curves to adjust one copy number gain/losses to height +1/-1 if applicable
CNV.data.scaled <- Scaling_CNV(V7Alt = CNV.data, 
                               n.TestCells = No.test, 
                               scaling.factor = 1)

## In case one needs to re-scale data can change the scaling.factor and run  
## the Scaling.CNV function again
CNV.data.scaled <- Scaling_CNV(CNV.data, 
                               n.TestCells = No.test, 
                               scaling.factor = 0.4)

## Defining M_NF as Noise-Free Matrix of test cells (and control cells)
M_NF <- CNV.data.scaled

#######  Noise Filteration after scaling
noise.thr = 0.4   # Noise threshold
for(w in 1:ncol(M_NF) ){
  for(j in 1:nrow(M_NF)){
    if( (M_NF[j,w] > -noise.thr) && (M_NF[j,w] < noise.thr)  ){
      M_NF[j,w] <- 0
    }
  }
}
M_NF <- as.matrix(M_NF)

#### Taking Square Rate
for(w in 1:ncol(M_NF)){
  for(j in 1:nrow(M_NF)){
    if (M_NF[j,w] > 0){
      M_NF[j,w]<- sqrt( as.numeric(M_NF[j,w]))
    } else if (M_NF[j,w] < 0){
      M_NF[j,w]<- -sqrt( -as.numeric(M_NF[j,w]))
    }
  }
}

rownames(M_NF) <- rownames(CNV.data.scaled)
colnames(M_NF) <- c(colnames(CNV.data.scaled)[-length(colnames(CNV.data.scaled))], "AveTest")

## Assigning chromosome number to each gene sorted based on chromosme number, 
## starts and ends to sketch the average iCNV curve of test cells

Gen.loc <- read.table("./Dataset/10XGenomics_gen_pos_GRCh38-1.2.0.txt", sep = '\t', header=TRUE)
Specific_genes <- which( Gen.loc[, 1] %in% rownames(CNV.data.scaled))
Assoc.Chr <-  Gen.loc[Specific_genes, 2]
#Assoc.Chr <-  apply(Assoc.Chr, 2, as.numeric)

#### Finalizing the iCNV-matrix by attaching gene-name and chromosome number lists
M_NF1 <- cbind(as.matrix(Gen.loc[Specific_genes, 1]), 
               as.matrix(Gen.loc[Specific_genes, 2]), 
               M_NF)
colnames(M_NF1) <- c("Genes", "Chromosomes", colnames(norm.data), "Ave test")
M_NF2 <- as.matrix(M_NF1)
M_NF3 <- as.matrix(M_NF2[ , ncol(M_NF2)] )
rownames(M_NF3) <- as.matrix(M_NF2[ , 1] )
#pdf( paste("AveiCNVcurve_testVScontrol_scaling factor",scaling.factor,"_Noise threshold:",noise.thr,".pdf", sep=""),
#     width = 6, height = 4, paper = 'special')

Sketch_AveCNV( Ave.mat = M_NF[, ncol(M_NF)] )

#---------------------------------------------------------
#--------------------------------------------------------------
##  Calculating tumor CNV score for test and control cells
#--------------------------------------------------------------
TotScore <- CNV_score( M_nf = M_NF )

## sketching tumor scores for all cells showing segregation of test/control tumor scores
TotScoreSort0 <- sort(TotScore[1,1:No.test])
TotScoreSortNPC <- sort(TotScore[1,(No.test+1):ncol(TotScore)])

colors <- c( "royalblue1","brown1")
Labels <- c("Control","Test")
names <- as.matrix(c(rep(Labels[1], length(TotScoreSortNPC)), rep(Labels[2],length(TotScoreSort0) )) )
value <- c(as.matrix(TotScoreSortNPC), as.matrix(TotScoreSort0))
data=data.frame(names,value)
data$factor <- factor(data$names, levels=Labels)
##
graphics.off()
plot.new()
par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
par(bty="l")
boxplot(data$value ~ data$factor , col=alpha(colors,0.6),
        ylim=c(min(TotScore)-1,max(TotScore)+1), 
        cex.lab = 2, cex.axis = 2, cex.main=2,
        xlab = "Type of individuals"
        , ylab = "CNV-score",
        bty='l', boxcol="gray" ,
        outpch=16, outcex=1)
title("CNV-score of individuals", col.main = "brown", cex.main = 2.5)
##
mylevels<-levels(data$factor)
levelProportions<-summary(data$factor)/nrow(data)
for(i in 1:length(mylevels)){
  thislevel<-mylevels[i]
  thisvalues<-data[data$factor==thislevel, "value"]
  myjitter<-jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(myjitter, thisvalues, pch=20, cex=2, xaxt = "n",yaxt = "n",col=alpha(colors[i], 0.6)) 
}

#--------------------------------------------------------------
##  Heatmap of CNV-curves and detecting rare subclones
#--------------------------------------------------------------
CNV.matrix <- t( M_NF[, -ncol(M_NF)])   
rownames(CNV.matrix) <- colnames(M_NF[, -ncol(M_NF)])
colnames(CNV.matrix) <- rownames(CNV.data)

###### generating heatmap
## In order to use genomic locations in our heatmap we read Gen.Loc matrix 
## with list of genes, chromosome numbers, starts and ends:

## Heatmap against list of genes

break.glist <- rep(0, 24)
break.glist <- heatmap_break_glist(CNV.mat2 = CNV.matrix )

CNV_htmp_glist( CNV.mat2 = CNV.matrix,
                Gen.Loc = Gen.Loc,
                clustering = FALSE,        
                sorting = TRUE,        
                CNVscore = TotScore,
                break.glist = break.glist,
                No.test = No.test )

## Heatmap against genomic location
break.gloc <- rep(0, 24)
break.gloc <- heatmap_break_gloc()

CNV_htmp_gloc( CNV.mat2 = CNV.matrix,
               Gen.Loc = Gen.Loc,
               clustering = FALSE,
               sorting = TRUE,        
               CNVscore = TotScore,
               break.gloc = break.gloc,
               No.test = No.test )


# CNV.mat2 = CNV.matrix; clustering = FALSE; sorting = TRUE; CNVscore = TotScore; cluster.lines = NULL










