# sciCNV

# A Short Explanation to sciCNV Pipeline
## Introduction
An important challenge in copy number variation calling using single-cell RNA-seqs is 
to provide a clean data which can reflect real biological interface and actual functionality 
of genes across entire genome and in comparison to a negative control. This can be 
even more challenging if we merge two datasets (e.g. test samples to control samples).
Here, we introduce a new bioinformatics tool which is prepared to overcome this problem and to provide an accurate 
inferred-copy-number-alteration analysis to study the evolutionary pathogenesis of diseases.
Our analysis can be widely applied to any malignancy/abnormality and in any context of intra-tumoral/malignant/infectious 
heterogeneous system.

Single-cell-inferred Copy Number Variation (sciCNV) pipeline is a novel 
strategy to likely answer all these challenges. Our pipe is consist of the 
following steps:

* [Reading data and quality control](#Reding_and_QC)
* [RTAM1/RTAM2 Normalization](#RTAM1/RTAM2-Normalization)
    * [Cleaning data](##Cleaning_data)
    * [RTAM1/2 normalization](##RTAM1/2)
    * [Justification of normalized data](##justification-of-normalized-data)
* [Clustering to cell-types](#clustering-to-celltypes)
* [iCNV Analysis from RNA-seq data](#infered-CNV-analysis)
    * [generating infered-CNV curves for (test and/or control) cells ](##sciCNV-on-normalized-data)
    * [Scaling and Filtering noise of the iCNV curves ](##Scaling_Noise_Filtering)
    * [Sketching the average MMPCs iCNV-curve after correction](##Sketching_ave_iCNV)
* [Malignancy CNV-score](#malignancy_score)
* [Heatmap of CNV-curves and detecting rare subclones](#heatmap)
    * [Generating heatmap](##generating_heatmap)
    * [Detecting subclones](##deteccting_subclones)




***
# Reading data and quality control

```
.rs.restartR()
options(stringsAsFactors = FALSE)

library(grid)
library(BiocGenerics)
library(Rtsne)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(heatmap3)
library(colorspace)
library(GetoptLong)
library(MASS)
library(Matrix)
library(mgcv)
library(RColorBrewer)
library(viridis)
library(gplots)
library(devtools)
library(Seurat)
library(GMD)
library(cluster)
library(robustbase)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

```

## Reading raw data with a list of genes on the first column

Starting with raw data in which gene symbols are as features and cell identities as column names, 
we firstly try to exclude damaged cells when the mitochondrial gene expressions of those cells are relatively higher 
than a specific threshold than can be different in different cell-types.

```
raw.data1 <- read.csv(file.choose(), header=TRUE)  
raw.data2 <- raw.data1[ , -1]
rownames(raw.data2) <- raw.data1[ , 1]
W <- ncol(raw.data2)
Col_Sum <- t(as.numeric(colSums(raw.data2)))
```

***
## Quality Control (QC): Eliminating damahged cells 

### Reading UMI and mitochondrial gene expressions associated to the data


```
nUMI <-  read.csv(file.choose(), header=TRUE)
nUMI <- as.matrix(nUMI[,-1])

mito.genes <-  read.csv(file.choose(), header=TRUE)
mito.genes <- as.matrix(mito.genes[,-1])
percent.mito.G <- t(as.matrix(colSums(mito.genes)))/Col_Sum

nGene1 <- matrix(0, ncol = ncol(raw.data2) , nrow = 1)
nonzero <- function(x) sum(x != 0)

nGene1 <- lapply( 1:ncol(raw.data2), nonzero(raw.data2[, i]) ) 
colnames(nGene1) <- colnames(raw.data2)
```

***
## Removing damaged cells accross entire population 

### Reading the 10x matrix

```
seurat.data <- Read10X("~/Desktop/Sample1/filtered_feature_bc_matrix/")

dense.size <- object.size(as.matrix(seurat.data))
dense.size
sparse.size <- object.size(seurat.data)
sparse.size
dense.size/sparse.size
seurat.obj <- CreateSeuratObject(seurat.data ,project = "sample1", min.cells = 0,
                           min.genes = 0, is.expr = 0, normalization.method = NULL,
                           scale.factor = 10000, do.scale = FALSE, do.center = FALSE,
                           names.field = 1, names.delim = "_", meta.data = NULL)
```

### Defining threshold to remove damaged cells

```
drop_mitoMads <- 3
threshold <- max(0.05, mean(percent.mito.G) + (drop_mitoMads)*(mad(percent.mito.G)) ) 

#--- Mitocchondrial percentage vs nUMI

layout(matrix(c(2,1,0,3,5,4,0,6),2,4) ,c(4.5,1,4.5,1),c(1,5), respect = TRUE) 
par(mar=c(3,3,0,0),mgp=2:0)
plot(percent.mito.G ~ nUMI, col=alpha("black",0.2),  pch=16, cex=1.2,
     xlab="nUMI", ylab="Mitochondrial expression (%)" , ylim=c(0,0.2)   # ylim=c(0,1)
)
with(seurat.obj, abline(h = threshold, lwd = 2, lty = 2, col = alpha("red", 0.8)))
legend("topright", bty = "n", lty = 2, col = alpha("red", 0.8), pt.bg = alpha("red", 0.8),
       legend=paste(drop_mitoMads, "MADs above mean :", threshold))

par(mar=c(0,3,1,0))
HST <- hist(nUMI,breaks=100,col="grey",main=NULL,xaxt="n")
MTPLR <- HST$counts / HST$density
Dnsty <- density(nUMI)
Dnsty$y <- Dnsty$y * MTPLR[1]
lines(Dnsty,col=alpha("blue",0.7))

par(mar=c(3,0,0,1))
Dnsty <-  density(as.matrix(percent.mito.G)) 
HST <- hist(percent.mito.G ,breaks=100,plot=F)
BAR <- barplot(HST$density,horiz=T,space=0,col="grey",main=NULL,xlab="Density")
SLOPE <- (max(BAR) - min(BAR)) / (max(HST$mids) - min(HST$mids))
lines(y=Dnsty$x * SLOPE + (min(BAR) - min(HST$mids) * SLOPE),
      x=Dnsty$y,lwd=2,col=alpha("blue",0.7))


#--- Selecting damged cells
damaged_cells <- as.matrix(which(percent.mito.G > threshold) )

par(mar=c(3,3,0,0),mgp=2:0)
plot(nGene1 ~ nUMI, 
     col=alpha("black",0.2),  
     pch=16, cex=1.2, xlab="nUMI", ylab="nGene")
points(nGene1[1, damaged_cells] ~ nUMI[1, damaged_cells] ,
       pch=21,cex=1.2,col=alpha("red",0.5),bg=alpha("red",0.3))
legend("topleft",bty="n",pch=21,col=alpha("red",0.8),pt.bg=alpha("red",0.8),
       legend="Damaged cells")


par(mar=c(0,3,1,0))
HST <- hist(nUMI,breaks=100,col="grey",main=NULL,xaxt="n")
MTPLR <- HST$counts / HST$density
Dnsty <- density(nUMI)
Dnsty$y <- Dnsty$y * MTPLR[1]
lines(Dnsty,col=alpha("blue",0.7))


par(mar=c(3,0,0,1))
Dnsty <-  density(as.matrix(nGene1)) 
HST <- hist(nGene1 ,breaks=100,plot=F)
BAR <- barplot(HST$density,horiz=T,space=0,col="grey",main=NULL,xlab="Density")
SLOPE <- (max(BAR) - min(BAR)) / (max(HST$mids) - min(HST$mids))
lines(y=Dnsty$x * SLOPE + (min(BAR) - min(HST$mids) * SLOPE),
      x=Dnsty$y,lwd=2,col=alpha("blue",0.7))


title("Detecting damaged cells based on mitochonrial expression level",
      outer = TRUE, line = -2,
      cex = 2, col ="blue")

```

***
##  Excluding Damaged Cells


```
Outliers <- damaged_cells
raw.data <- raw.data2[, - Outliers]
rownames(raw.data) <- raw.data1[ , 1]
dim(raw.data2)
```

***
#RTAM1/RTAM2 Normalization


```
norm.data <- RTAM_normalization(raw.data, method = "RTAM2", Min_nGn = 250, Optimizing = TRUE) 
```

## Sketching non-zero expressions


```
graphic.off()
plot.new()
par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
par(xaxs="i", yaxs="i") 
par(bty="l")
plot(matrix(1,ncol=1,nrow=nrow(as.matrix(norm.data[,1][norm.data[,1]>0]))), 
     log2(as.matrix(norm.data[,1][norm.data[,1]>0]) +1  ), pch=16, cex=0.3, 
     col="darkgray" ,
     xlim=c(-100, ncol(norm.data)+100)
     ,ylim=c(0,  4.2), xlab = "Cells (ranked by UMI)",
     ylab = expression("Expressions ("*Log[2]*"(.+1 ))"), 
     cex.lab = 2, cex.axis = 2, cex.main=2)


for(i in 2:ncol(norm.data)){
  par(new=TRUE)
  points(  matrix(i,ncol=1,nrow=nrow(as.matrix(norm.data[,i][norm.data[,i]>0]))), 
           log2(as.matrix(norm.data[,i][norm.data[,i]>0]) +1  ), pch=16, cex=0.3, 
           col="darkgray")
}
title( paste("Sample1, ",method, "-normalization, cutoff",Min_nGn," - nGene",Nt,sep=""), 
       col.main = "brown", cex.main = 2)
       
```     

### Checking the balance of 95% commonly expressed genes

```
W3 <- ncol(norm.data)
XX1 <- seq(1,W3,1)
Z3 <-   as.matrix(norm.data[which(rowSums(norm.data[,XX1] != 0) > 0.95*W3 ), ] )

library(robustbase)
COLMED <- log2(colMeans(as.matrix(Z3[Z3>0])) +1)
COLMED <- t(as.matrix(COLMED) )

for(j in 1:(W3)){
  par( new=TRUE)
  points(matrix(j,ncol=1,nrow=1),
  log2(mean(as.matrix(Z3[,j][Z3[,j]>0])) + 1 ) ,axis = FALSE , col="red" , pch=15, cex =0.5      
  )
} 
```

### Checking the balancce of average expression of housekeeping genes

```
Houskeeping_gene_list <- read.csv(file.choose(), header=TRUE)
W3 <- ncol(norm.data)

BB12  <- norm.data[which(rownames(norm.data)%in% t(as.matrix(Houskeeping_gene_list))), ]
BB12 <- as.matrix(BB12)
colnames(BB12) <- colnames(norm.data)
dim(BB12)

Mean_BB12 <- matrix(0, ncol= W3 , nrow= 1)
for(k in 1: (W3)){
Mean_BB12[1,k] <- as.numeric(mean(BB12[,k][BB12[,k]>0]) )
}
par( new=TRUE)
points(log2(Mean_BB12[1,] +1 ) ,axis = FALSE , col="blue" , pch=15, cex =0.5 )

legend("bottomright",bty="n",pch=16,col=c("blue",NA), cex=1.5,
       legend=paste("Mean of Houskeeping gene expressions"))
```

***
# Clustering to cell-types


***
##  Excluding non-expressed genes

```
train1 <- as.matrix(norm.data)
rownames(train1) <- rownames(Scaled_Normalized_log)
colnames(train1) <- colnames(Scaled_Normalized_log)

train <- as.matrix(train1[which(rowSums(train1 != 0) >= 1  ), ] )
dim(train)
```

***
##  clustering the normalized data

```
library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(svd)
library(ggplot2)
library(ggridges)
library(viridis)
require(scales)
```

## Finding matrix of single cells (MSC) as a seurat object for clustering

```
MSC <- CreateSeuratObject( train ,project = "RTBM241", min.cells = 0,
                           min.genes = 0, is.expr = 0, normalization.method = NULL,
                           scale.factor = 10000, do.scale = FALSE, do.center = FALSE,
                           names.field = 1, names.delim = "_", meta.data = NULL)
## OR 
MSC  <- CreateSeuratObject(counts = train, project = "sample1", normalization.method = NULL)

MSC<- FindVariableGenes(object = MSC, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0, x.high.cutoff = 20, y.cutoff = -5)

length(x = MMS@var.genes)


MSC <- ScaleData(object = MSC)


MSC <- RunPCA(object = MSC, pc.genes = MSC@var.genes, 
              do.print = TRUE, pcs.print = 1:10, genes.print = 10)
MSC <- ProjectPCA(object = MSC, do.print = FALSE)

MSC <- FindClusters(object = MSC, 
                    reduction.type = "pca",
                    dims.use = 1:10, 
                    resolution = 0.6, print.output = 0, save.SNN = TRUE,
                    force.recalc = TRUE)

PrintFindClustersParams(object = MSC)

MSC <- RunTSNE(object =  MSC , reduction.type = "pca"   #  "cca.aligned"
               , dims.use = 1:10, do.fast = TRUE )
```

***
# iCNV Analysis of RNA-seq data

Now we run sciCNV function on a dataset comprising of test and control cells to derive iCNV-curves per cell across 
entire genome. So we firstly read the normalized matrix of test and control cells which is called _norm.tst.ctrl_ matrix,

```
norm.tst.ctrl1 <- read.csv( file.choose(), header = TRUE) 
tst.labels <- sapply(strsplit(colnames(Mn)[1], split='_', fixed=TRUE), function(x) (x[1]))
L_tst <- length(grep(tst.labels, colnames(Mn)))

ctrl.index <- seq( L_tst + 1, ncol(norm.tst.ctrl), 1)   # No of controcl cells
tst.index <- seq(1, L_tst , 1) # No of test cells
```

## generating infered-CNV curves for (test and/or control) cells 

Then we perform our scciCNV function to generate iCNV curves for both tests and control cells.

```
CNV.data <- sciCNV (norm.tst.ctrl, n.TestCells = L_tst,  sharpness = 50, baseline_adj = FALSE)
```

***
## Scaling and Filtering noise of the iCNV curves 

In here, we scaling iCNV-curves to adjust one copy number gain/losses to height +1/-1 or so if applicable.

```
CNV.data.scaled <- Scaling.CNV(CNV.data, n.TestCells = L_tst, scaling.factor = 1 )

CNV.scaled.output <- Scaling.CNV(CNV.data, n.TestCells = L_tst, scaling.factor = 0.64)
CNV.scaled <- CNV.scaled.output[1]
Ave_MMPCs1 <- CNV.scaled.output[2]
```

### defining M_NF as Noise-Free Matrix of test cells (and control cells)

```
M_NF <- matrix(0, ncol=ncol(CNV.scaled)+1, nrow=nrow(CNV.scaled))
M_NF[, seq(1,ncol( CNV.scaled))] <- CNV.scaled
M_NF[, ncol(M_NF)] <- Ave_MMPCs1
```

### Noise Filteration after scaling

Based on the average bulk iCNVs of test cells, one may want to remove redundant signals.

```
noise.thr = 0.3   # Noise filteration threshold for filtering noise after sclainng
for(w in 1:ncol(M_NF) ){
  for(j in 1:nrow(M_NF)){
    if( (M_NF[j,w] > -noise.thr) && (M_NF[j,w] < noise.thr)  ){
      M_NF[j,w] <- 0
    }
  }
}
MM_NF <- as.matrix(MM_NF)
```

### Taking Square Root of iCNVs

At this stage we take square root of all iCNV values to converge to 1 those values which are 
closer values to 1 and make smaller values than 

```
for(w in 1:ncol(M_NF)){
  for(j in 1:nrow(M_NF)){
    if (M_NF[j,w] > 0){
      M_NF[j,w]<- sqrt( as.numeric(M_NF[j,w]))
    } else if (M_NF[j,w] < 0){
      M_NF[j,w]<- -sqrt( -as.numeric(M_NF[j,w]))
    }
  }
}
```

Then we assigning chromosome number to each gene sorted based on chromosme number, starts and ends 
to sketch the average iCNV curve of test cells.

```
Gen.Loc <- read.csv( "~/Desktop/10X_Project/NEW_10XGENOMICS ANALYSIS/10XGenomics_gen_pos_GRCh38-1.2.0.csv", header=TRUE)
Assoc.Chr <-  as.matrix(Gen.Loc[rownames(norm.data), 2])
```

### Finalizing the iCNV-matrix by attaching gene-name and chromosome number lists

```
M_NF1 <- cbind(Gen.Loc[rownames(norm.data), c(1,2)], M_NF)
colnames(MM_NF1) <- c("Genes", "Chromosomes", colnames(norm.tst.ctrl), "Ave MMPCs")
```

## Sketching the average MMPCs iCNV-curve after correction

```
M_NF <- as.matrix(M_NF1)
M_NF2 <- t(as.matrix(M_NF[ , ncol(M_NF)] ))
rownames(MM_NF2) <- as.matrix(M_NF[ , 1] )

pdf( paste("AveiCNVcurve_testVScontrol_scaling factor",scaling.factor,"_Noise threshold:",noise.thr,".pdf", sep=""),
     width = 6, height = 4, paper = 'special')

Sketch.AveCNV( Ave.mat = M_NF2, Assoc.Chr )
dev.off()
```

***
# Malignancy CNV-score

Perfoming _CNV.score_  function we calculate malignancy CNV-score for test and control cells, which shows the likeness of 
test and control cells to the average iCNV-curve of bulk test cells.

```
TotScore <- CNV.score( M_NF )
```

## sketching tumor scores for all cells showing segregation of test/control tumor scores

```
TotScoreSort0 <- sort(TotScore[1,1:(ncol(TotScore)-205)])
TotScoreSortNPC <- sort(TotScore[1,(ncol(TotScore)-204):ncol(TotScore)])

colors <- c( "green3","brown1")
Labels <- c("NPCs","CLS0")


plot.new()
par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
par(bty="l")
plot(1:205, TotScoreSortNPC ,col=alpha(colors[1],0.4) ,
     xlab="Cells", ylab="Tumor Score", 
     xlim=c(-5,ncol(TotScore)+10), 
     ylim = c(min(TotScore),max(TotScore)+10),
     pch = 16, 
     cex.lab = 2, cex.axis = 2, cex.main=2, cex=2)
par(new=TRUE)
points(206:ncol(TotScore),TotScoreSort0,col=alpha(colors[2],0.4) ,axes=FALSE, pch = 16, cex=2  )

title('Tumor karyotype likness',col.main="brown", cex.main=2.5)

legend(1000,-20,bty="n",pch=16,col=c(colors[1],NA),legend=paste(Labels[1]), 
       cex=1.7)
legend(1000,0,bty="n",pch=16,col=c(colors[2],NA),legend=paste(Labels[2]), 
       cex=1.7)
```

***
# Heatmap of CNV-curves and detecting rare subclones

To see the final result of iCNV-analysis we sketch all CNV-curves together and try to 
segregate diverse subclones based on their CNV-similarities.

```
CNV.mat <- t( M_NF)   #  /quantile(as.matrix(M1),0.99) )  #max(M1))
L_ctrl <- nrow(CNV.mat) - L_tst

tst.score <- sort(TotScore[1,1:(ncol(TotScore)-L_ctrl)] , decreasing=TRUE)     #MMPCs
ctrl.score <- sort(TotScore[1, (ncol(TotScore)-L_ctrl+1):ncol(TotScore)] , decreasing=TRUE)  #NBCs
ranked.score <- as.matrix( c(colnames(ctrl.score), colnames(tst.score)) )

CNV.mat1 <- as.matrix(CNV.mat[ t(as.matrix(ranked.score)) , ])
rownames(CNV.mat1) <-  ranked.score     
colnames(CNV.mat1) <- rownames(M_NF)

CNV.mat1_saved <- CNV.mat1
ROWlist <- rownames(CNV.mat1)
COLlist <- colnames(CNV.mat1) 
```
## Generating heatmap

Here is the function to generate heatmap using _heatmap3_ funcction against either 
list of genes or genomic locations.

```
CNV.heatmap(CNV.mat2,   sorting= TRUE,  TotScore,   L_tst,  heatmap.type =  "gnloc")
```

## Detecting rare subclones

CNV-similarities may separate cells in a different way in comparison to clustering based 
on gene expressions. 










