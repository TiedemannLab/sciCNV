

#######################################################################################
######                             Mito_umi_gn function                         ####### 
###### Function to sketch mitochondrial transcript content, nUMI and nGene data #######
######      at cellular level to detect damaged cells within the population     #######
######  Tiedemann Lab - Princess Margaret Cancer Centre, University of Toronto  #######
######                       copyright@AliMahdipourShirayeh                     #######
#######################################################################################

# Please refer to the reference and supplemental materials described in the README for additional details.
#
# Definitions:
# mat: raw data matrix 
# percent.mito.G: percentage of mitochondrial transcripts relative to total number of reads per cell
# nUMI: number of unified molecular identification (UMI)
# nGene: number of expressed genes
# No.test: number of test cells included in the data; can be used to delineate the populations of test and control cells in the heatmap
# drop.mads: number of the median absolute deviation (MAD) to drop damaged cells using the percentage of mitochondrial expression level

Mito_umi_gn <- function( mat = MMS, 
                         percent.mito.G = percent.mito.G,
                         nUMI = nUMI,
                         nGene = nGene1,
                         No.test,
                         drop.mads = 3
){
  
  if( missing(percent.mito.G)){
    stop("Please insert proper percentage of mitochondrial matrix (percent.mito.G).")
  }
  if( missing(nUMI)){
    stop("Please insert proper number of UMI per cell matrix (nUMI).")
  }
  if( missing(nGene)){
    stop("Please insert proper number of gene per cell matrix (nGene).")
  }
  if( is.null(drop.mads)){
    drop.mads <- 3
  }
  
  threshold <- max(0.05, mean(percent.mito.G) + (drop.mads)*(mad(percent.mito.G)) ) 
  
  layout(matrix(c(2,1,0,3,5,4,0,6),2,4) ,c(4.5,1,4.5,1),c(1,5), respect = TRUE) 
  par(mar=c(3,3,0,0),mgp=2:0)
  plot(percent.mito.G ~ t(as.matrix(nUMI[1:No.test])), col=alpha("black",0.2),  pch=16, cex=1.2,
       xlab="nUMI", ylab="Mitochondrial expression (%)",
       cex.lab=1.5,cex.lim=1.5,cex.axis=1.5
  )
  with(mat, abline(h = threshold, lwd = 2, lty = 2, col = alpha("red", 0.8)))
  legend("topright", bty = "n", lty = 2, col = alpha("red", 0.8), pt.bg = alpha("red", 0.8),
         legend=paste(drop.mads, "MADs above mean :", threshold))
  
  par(mar=c(0,3,1,0))
  HST <- hist(t(as.matrix(nUMI[1:No.test])) ,breaks=100,col="grey",main=NULL,xaxt="n")
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
  
  
  #----- Selecting damged cells
  damaged_cells <- NULL
  damaged_cells <- as.matrix(which(percent.mito.G > threshold) )
  
  par(mar=c(3,3,0,0),mgp=2:0)
  plot( t(as.matrix(nGene1))[1:No.test] ~ t(as.matrix(nUMI))[1:No.test], 
        col=alpha("black",0.2),  
        pch=16, cex=1.2, xlab="nUMI", ylab="nGene",
        cex.lab=1.5,cex.lim=1.5,cex.axis=1.5)
  points( t(as.matrix(nGene1))[damaged_cells] ~ t(as.matrix(nUMI))[damaged_cells] ,
          pch=21,cex=1.2,col=alpha("red",0.5),bg=alpha("red",0.3))
  legend("topleft",bty="n",pch=21,col=alpha("red",0.8),pt.bg=alpha("red",0.8),
         legend="Damaged cells")
  
  
  par(mar=c(0,3,1,0))
  HST <- hist(t(as.matrix(nUMI))[1:No.test],breaks=100,col="grey",main=NULL,xaxt="n")
  MTPLR <- HST$counts / HST$density
  Dnsty <- density(nUMI)
  Dnsty$y <- Dnsty$y * MTPLR[1]
  lines(Dnsty,col=alpha("blue",0.7))
  
  
  par(mar=c(3,0,0,1))
  Dnsty <-  density(as.matrix(nGene1[1:No.test])) 
  HST <- hist( nGene1 ,breaks=100,plot=F)
  BAR <- barplot(HST$density,horiz=T,space=0,col="grey",main=NULL,xlab="Density")
  SLOPE <- (max(BAR) - min(BAR)) / (max(HST$mids) - min(HST$mids))
  lines(y=Dnsty$x * SLOPE + (min(BAR) - min(HST$mids) * SLOPE),
        x=Dnsty$y,lwd=2,col=alpha("blue",0.7))
  
  
  title("Detecting damaged cells based on mitochonrial expression level",
        outer = TRUE, line = -2,
        cex.main = 2, col.main ="brown")
  
  
  return(damaged_cells)
  
  
}






