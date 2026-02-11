options(stringsAsFactors = FALSE)

library(devtools)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(grid)
library(heatmap3)
library(colorspace)
library(GetoptLong)
library(RColorBrewer)
library(viridis)
library(gplots)
library(scales)
library(bitops)
library(RCurl)
library(scales)
library(lattice)

library(pheatmap)
library(grid)
library(RColorBrewer)
library(viridis)


library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library('biomaRt')

#M <- read.csv(file.choose(), header=TRUE)      # The V7alt curve matrix
#M1 <- train #as.matrix(MMS@meta.data)

#Mnorm <- read.csv(file.choose(), header=TRUE)  # The normalized version to calculate the BC/PC-score
Mnorm <- read.table(file.choose(), sep="\t", header=TRUE)  # The normalized version to calculate the BC/PC-score
dim(Mnorm) 

SMPL <- "sample_name"
comparison <- paste("+1Q vs No-1Q")

Mn <- as.matrix(Mnorm[,-1])
rownames(Mn) <- Mnorm[,1]

##### 1Q vs No-1Q

fold_val11=0.65
fold_val21=0.25
pval11=1E-5
pval21=1E-2

Original <- VolcanoPlot( Mn=Mn, 
             SMPL=SMPL, 
             comparison=comparison,
             gene.norm=FALSE,   # default is TRUE
             fold_val1=fold_val11,
             fold_val2=fold_val21,
             pval1=pval11,
             pval2=pval21
)


######
nUMI <- as.matrix(t(as.numeric(colSums(Mn))))
colnames(nUMI) <- colnames(Mn)

fold_val12=0.4
fold_val22=0.4
pval12=5*1E-3
pval22=5*1E-3

Balanced <- BalancedPlot( Mn=Mn, 
              SMPL=SMPL,
              nUMI=nUMI,
              comparison=comparison,
              gene.norm=FALSE,   # default is TRUE
              fold_val1=fold_val12,
              fold_val2=fold_val22,
              pval1=pval12,
              pval2=pval22
)

############################### HEATMAPS
## (A) Pre-balancing
HeatmapVolcano( Mn=Mn, 
                SMPL=SMPL, 
                comparison=comparison,
                gen.loc = FALSE,   
                fold_val1=fold_val11,
                fold_val2=fold_val21,
                pval1=pval11,
                pval2=pval21,
                gene.norm =FALSE,
                type = "original"  
)

## (B) Post-balancing
HeatmapVolcano(Mn=Mn, 
               SMPL=SMPL, 
               comparison=comparison,
               nUMI=nUMI,
               gen.loc = FALSE,  
               fold_val1=fold_val12,
               fold_val2=fold_val22,
               pval1=pval12,
               pval2=pval22,
               gene.norm =FALSE,
               type = "balanced"     
)
