

#######################################################################################
######                             Sketch_AveCNV function                       ####### 
######         Sketching the average test CNV-curve agaisnt genomic location    #######
######  Tiedemann Lab - Princess Margaret Cancer centre, University of Toronto  #######
######                        copyright@AliMahdipourShirayeh                    #######
#######################################################################################

Sketch_AveCNV <- function(Ave.mat){
  
  Gen.loc <- read.table("./sciCNV/Dataset/10XGenomics_gen_pos_GRCh38-1.2.0.csv", sep = ',', header=TRUE)
  Specific_genes <- which( as.matrix(Gen.loc)[, 1]   %in% rownames(as.matrix(Ave.mat)))
  Assoc.Chr <-  as.matrix(Gen.loc[Specific_genes, 2])

  plot.new()
  par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
  plot( Ave.mat, 
        col="brown1", 
        type="l", 
        lty=1, 
        lwd=3, 
        xlim=c(-nrow(as.matrix(Ave.mat))*.02,nrow(as.matrix(Ave.mat))*1.02),
        ylim=c(-2,2),
        cex=1, 
        cex.lab=2, 
        cex.axis=2, 
        xlab="Genomic location", 
        ylab="CNV-expression")
  abline(h=0, col="black")
  
  M <- as.matrix(Assoc.Chr) # Chromosome numbers
  Break <- matrix(0, ncol = 24, nrow = 1)
  for(i in 1: 22){ 
    Break[1,i] <- apply(as.matrix(Assoc.Chr) == i, 2, which.max)
  }
  if( length(which( as.matrix(Assoc.Chr) == "X") > 0) ){
    Break[1,23] <- apply(as.matrix(Assoc.Chr) == "X", 2, which.max)
    if( length(which( as.matrix(Assoc.Chr) == "Y") > 0) ){
      Break[1,24] <- apply(as.matrix(Assoc.Chr) == "Y", 2, which.max)
    } else{
      Break[1,24] <- nrow(as.matrix(Ave.mat))
    }
  } else{
    Break[1,23] <- nrow(as.matrix(Ave.mat))-1
    Break[1,24] <- nrow(as.matrix(Ave.mat))
  }
  
  
  Break_lines <- as.matrix(c(Break, nrow(Assoc.Chr)))
  abline(v=Break_lines, col="gray65", lty=2)
   par(new=TRUE); 
  points( Ave.mat, 
          col = "brown1", 
          type = "l", 
          lty = 1, 
          lwd = 3  
          )
   points(Break , matrix(c(1.8, 1.7 ), ncol=24, nrow=1), pch=16, col="royalblue1", cex=4)
   text(Break , matrix(c(1.8, 1.7 ), ncol = 24, nrow = 1), c(seq(1, 22, 1),"X","Y"), 
        col = "white", cex = 1.2)
   
  
  title("Ave iCNV-curve of test bulk relative to normal control", cex.main=2, col.main="brown")
  
}
