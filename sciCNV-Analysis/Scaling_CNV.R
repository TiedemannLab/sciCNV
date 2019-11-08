#######################################################################################
######                             Scaling.CNV function                         ####### 
######                 To scale CNV-curves derived for ssRNA-seq data           #######
######                       copyright@AliMahdipourShirayeh                     #######
######                 Tiedemann Lab - Princess Margaret Cancer centre          #######
#######################################################################################

Scaling_CNV <-  function(V7Alt, n.TestCells, scaling.factor ){
  
  # argument validation
  if ( Reduce("|",is.na(V7Alt)) ){
    stop( "The CNV-matrix has not been properly prepared.")
  }
  
  if ( is.na( n.TestCells) ){
    stop( "Number of tumor cells needs to be inserted.")
  }
  
  if ( is.na(scaling.factor) ){
    scaling.factor <- 1.0
  }
  
  Ave_control <- matrix(0, ncol=1, nrow= nrow(V7Alt ))
  Ave_test <- matrix(0, ncol=1, nrow= nrow(V7Alt ))
  Z2 <- seq( n.TestCells + 1, ncol(V7Alt), 1)   #NPCs
  Z3 <- seq(1, n.TestCells , 1)     #MMPCs
  
  ###### Sample 3 -choosing ONLY PCs and NBCs/NPCs (205 cells)
  for(i in 1:nrow(V7Alt)){
    if( !sum(V7Alt[ i, Z2]) == 0 ){
      Ave_control[i, 1] <- mean(as.matrix(V7Alt[ i, Z2][which(!V7Alt[ i, Z2] == 0 ) ]))
    } else {
      Ave_control[i, 1] <- 0
    }
    Ave_test[i, 1] <- mean(as.matrix(V7Alt[ i, Z3]))
  }
  
  ## Scaling data if needed
    V7Alt <- V7Alt*(1/scaling.factor)
    Ave_control <- Ave_control*(1/scaling.factor)
    Ave_test <- Ave_test*(1/scaling.factor)
    
    Ave_control[is.na(Ave_control)] <- 0
    Ave_test[is.na(Ave_test)] <- 0
  
    
  
  ## To scale data we used some genes that are correlated to specific number of CNV
  plot.new()
  par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
  plot(Ave_control[, 1], 
       col= "blue", 
       pch=16, 
       cex=0.5,
       ylim=c(-2,2), 
       cex.lab=2,
       cex.axis=2, 
       xlab="Genomic location",
       ylab="CNV-expression"
  )
  points(Ave_test[, 1], col= "red", pch=16, cex=0.5 )
  abline(h=0, col="black")

  legend("bottomleft",
         legend=c( "Ave control","Ave test"), 
         inset=.02, col=c("blue","red"), 
         fill=c("blue", "red"), 
         horiz=TRUE, 
         cex=1.2)
  
  title(paste("Average iCNV curves of test/control cells - Threshold=",
              round(scaling.factor, 
              digits=2)), 
              col.main="brown", 
              cex.main=2)
  
  Final_Mat <- cbind(V7Alt, Ave_test)
  rownames(Final_Mat) <- rownames(V7Alt)
  colnames(Final_Mat) <- c(colnames(V7Alt),c("AveTest"))
  #### this function returnes the scaled CNV-matrix
  return( Final_Mat )
  
  
  
}

