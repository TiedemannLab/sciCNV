




#######################################################################################
######                             CNV_infer function                           ####### 
######          for Single Cell Inferred Copy Number Variation (sciCNV)         #######
######  Tiedemann Lab - Princess Margaret Cancer centre, University of Toronto  #######
######                       copyright@AliMahdipourShirayeh                     #######
#######################################################################################

# Please see the reference and supplmental materials described in the README file for additional information.
#
# Definitions:
# ss.expr: scRNA-seq based expression for each test/control cell
# mean.ctrl: The average gene expression of control cells
# chr.n: List of chromosome numbers associated with the list of genes
# Resolution: Adjusts the resolution, nrow(MSC)/(50*sharpness), used for the sciCNV-curve calculation. Default sharpness is =1.0 (best sharpnesses range between 0.6-1.4).
# baseline: An optional correction to adjust the CNV zero setpoint (copy number gain =0) which is otherwise the median CNV of all genes. 
#           Default = 0 and typical range is between -0.5 and 0.5.
# baseline_adj: The baseline adjustment is only applied to test cells if it is TRUE. Default is FALSE.
# P12: The variable which is associated to resolution, is automatically calculated and is used for CNV calling.
# mat.fab: is a merged matrix of FF, AW, BD vectors per cell. These vectors are used to calculate the relative expression of test to control. 

CNV_infer  <- function( ss.expr , 
                        mean.ctrl,
                        gen.Loc , 
                        resolution, 
                        baseline_adj = FALSE,
                        baseline = 0,
                        chr.n,
                        P12,
                        mat.fab
){ 
  
  ## argument value/validation
  if (is.na(baseline)){
    baseline <- 0    
  }
  
  if (is.na(baseline_adj)){
    baseline_adj = FALSE    
  }
  
  if (is.na(sharpness)){
    sharpness <-  1.0   
  } else if( (sharpness > 1.4) || (sharpness < 0.6)){
    stop("Sharpness is not in the range [0.6, 1.4].")
  }
  
  # argument correction for Baseline Correction (BLC)
  if ( baseline_adj == "TRUE" ){
    BLC <- baseline          # Baseline Correction value
  } else {
    BLC <- 0
  }
  
  ## defining parameters        
  P31 <- P12/2;      E22 <- 1
  E23 <- 1;          O51 <- 0.985
  O52 <- 1 - O51;    P31 <- P12/2
  E22 <- 1;          E23 <- 1
  O25 <- 0.78;       O26 <- 0.2
  C297 <- resolution/10  
  C298 <- 2;         J311 <- 1
  m <- 1.0001
  t <-  0.0002
  Lambda <- 0.00001
  row_n <- nrow(ss.expr)
  
  ## Recalling FF, AW, BD vectors and defining G matrix to start with average relative expression of test to control
  G <- as.matrix(seq(1,row_n,1))
  FF <- mat.fab[,1]
  AW <- mat.fab[,2]
  BD <- mat.fab[,3]
  
  #---------------  T: moving average of gene expression in control cells & U: moving average of expression in test cell
  TT <- rep(0, row_n); U <- rep(0, row_n)
  max_1 <- pmax( as.matrix(AW - P12), matrix(-floor(P12/C297), row_n) )
  min_1 <- pmin( as.matrix(P12 - BD), matrix(floor(P12/C297), row_n) )
  
  for( i in 1:row_n ){
    bb <- seq(i + max_1[i], i + min_1[i] , 1)
    bb <- bb[! is.na(ss.expr[bb])]
    TT[i] <- mean( ss.expr[bb] )
    U[i] <- mean( mean.ctrl[bb] )
    
  }
  
  
  #---------------------- V: ma-Wratio (weighted disparity score comparing gene expression in test and control cells within the chromosome segment defined by the moving averages)
  V <- (TT - U)/(TT + U + Lambda )
  
  P100 <- quantile( as.matrix(V)[seq(1,row_n,1)], 0.5 - BLC )
  WW <- V - P100
  
  
  #---------------  X: moving average of test cell, before weighting. & Y: moving average of control cells, before weighting 
  X <- rep(0, row_n); Y <- rep(0, row_n)
  max_2 <- pmax( as.matrix(AW - P12), matrix(-floor(P12/C298), row_n) )
  min_2 <- pmin( as.matrix(P12 - BD), matrix(floor(P12/C298), row_n) )
  
  for(i in 1:row_n){
    bb <- seq(i + max_2[i], i + min_2[i], 1)
    bb <- bb[! is.na(ss.expr[bb])]
    X[i] <- mean( ss.expr[bb] )
    Y[i] <- mean( mean.ctrl[bb] )
  }
  
  #---------------------- Z: ma-Wratio (weighted expression disparity score comparing gene expression in test and control cells within the chromosome segment defined by the moving averages) 
  #-----------------------AA: median normalized Z. AB: Linearized solution of AA from quadratic equation (aka W-lin)
  Z <- (X - Y)/(X + Y + Lambda)
  
  P101 <- median( Z[seq(1,row_n,1) ]) 
  AA <- Z - P101
  AB <- 3.3222*(AA)^4 + 5.6399*(AA)^3 + 4.2189*(AA)^2 + 3.8956*(AA)
  
  #---------------------- ACC: ma-DMS ("Demoncratic moving score" [1 CNV vote per gene]- calculated on post moving average post median normalized data)
  
  ACC <- rep(0, row_n)
  
  ACC[1] <- 0
  for(i in 2:row_n){
    if( WW[i] > 0){ 
      ACC[i] <-  ACC[i-1] + 1
    } else if(WW[i] < 0){ 
      ACC[i] <-  ACC[i-1] - 1
    } else
      ACC[i] <- ACC[i-1]
  }
  
  #---------------  NEW: AC: ma-DMS (on post MovAve median)  &  DA: ma-Wratio  non-Moving Ave. +/-  
  AC <- ACC + ( BLC * G )
  DA <- (ss.expr - mean.ctrl) %/% (ss.expr + mean.ctrl + Lambda)
  
  
  #---------------------- AC: non-ma DMS
  V56 <- length( DA[ DA > 0 ] )
  V57 <- length( DA[ DA < 0 ] )
  V60 <- (V56-V57)/row_n
  DD <-  rep(0, row_n)
  
  DD[1] <- 0
  for(i in 2:row_n){
    if( DA[i] > 0){ 
      DD[i] <-  DD[i-1] + 1 - V60
    } else if(DA[i] < 0){ 
      DD[i] <-  DD[i-1] - 1 - V60
    } else
      DD[i] <- DD[i-1] - V60
  }
  
  
  SCALER <- (max(DD)-min(DD))/(max(ACC)-min(ACC))      # Scaling non-ma to ma-DMS
  DE <- DD/SCALER
  DG3 <- baseline
  DG <- DE + ( baseline * G)
  DH <- (AC*0.66) + (DG*0.34)
  
  #---------------------- AD: V7alt[ma- and non-ma] Centrograde : Slope of AC (DMS)
  DI <- rep(0, row_n)
  
  
  for(i in 1:row_n){
    if( (BD[i] < P31) && (AW[i] < P31) ){
      aaa <- seq(i - P31 ,i + P31,1)
      aaa <- aaa[! is.na(DH[aaa])]
      Ya <- as.matrix(DH[aaa])
      summary(lmResult <- lm(Ya~aaa))
      DI[i] <- as.matrix(coef(lmResult)["aaa"][[1]])
      
    } else if( (BD[i] < P31) && (AW[i] >= P31) ){
      if(  (AW[i] - P12) < P31 ){
        aaa <- seq(i + AW[i] - P12, i + P31 , 1)
        aaa <- aaa[! is.na(DH[aaa])]
        Ya <- as.matrix(DH[aaa])
        summary(lmResult <- lm(Ya~aaa))
        DI[i] <- as.matrix(coef(lmResult)["aaa"][[1]])
        
      } else
        aaa <- seq(i + AW[i] - P12, i + P31 , -1)
      aaa <- aaa[! is.na(DH[aaa])]
      Ya <- as.matrix(DH[aaa])
      summary(lmResult <- lm(Ya~aaa))
      DI[i] <- as.matrix(coef(lmResult)["aaa"][[1]])
      
    } else if( (BD[i] >= P31) && (AW[i] < P31)  ){
      if( ( - P31) < ( P12 - BD[i]) ){
        aaa <- seq(i - P31, i + P12 - BD[i] , 1)
        aaa <- aaa[! is.na(DH[aaa])]
        Ya <- as.matrix(DH[aaa])
        summary(lmResult <- lm(Ya~aaa))
        DI[i] <- as.matrix(coef(lmResult)["aaa"][[1]])
        
      } else
        aaa <- seq(i - P31, i + P12 - BD[i] , -1)
      aaa <- aaa[! is.na(DH[aaa])]
      Ya <- as.matrix(DH[aaa])
      summary(lmResult <- lm(Ya~aaa))
      DI[i] <- as.matrix(coef(lmResult)["aaa"][[1]])
      
    } else if( (BD[i] >= P31) && (AW[i] >= P31)  ){
      if( (P12 - BD[i]) < (AW[i] -P12) ){
        aaa <- seq(i + P12 - BD[i], i + AW[i] -P12 , 1)
        aaa <- aaa[! is.na(DH[aaa])]
        Ya <- as.matrix(DH[aaa])
        summary(lmResult <- lm(Ya~aaa))
        DI[i] <- as.matrix(coef(lmResult)["aaa"][[1]])
        
      } else if( (P12 - BD[i]) == (AW[i] - P12) ){
        DI[i] <- 0
      } else{
        aaa <- seq(i + P12 - BD[i], i + AW[i] - P12 , -1)
        aaa <- aaa[! is.na(DH[aaa])]
        Ya <- as.matrix(DH[aaa])
        summary(lmResult <- lm(Ya~aaa))
        DI[i] <- as.matrix(coef(lmResult)["aaa"][[1]])
      }
    }
  }
  
  
  
  #---------------------- AF: maV7-alt x maW-lin curve & AG: maV7-alt x maW-lin curve with Manual Scale
  AF <- rep(0, row_n)
  
  for(i in 1:row_n){
    if( (! is.na(DI[i]*AB[i])) && (DI[i]*AB[i] > 0) ){
      AF[i] <- DI[i]*abs(AB[i])
    } else
      AF[i] <- 0
  }
  
  AG <- AF/J311
  
  
  #### returns the square of the single cell inferred CNV profile matrix ('squared' estimate of CNV; requires square root)
  
  return(AG)
  
  
  
  
}


