
#######################################################################################
######                             CNV_infer function                           ####### 
######                 Infered Copy Numver Vatriation for each single cell      #######
######                       copyright@AliMahdipourShirayeh                     #######
######                 Tiedemann Lab - Princess Margaret Cancer centre          #######
#######################################################################################

# ss.expr: RNA-seq for each test/control cell
# mean.cntrl: The average expressions of control cells
# chr.n: List of chromosome numbers associated with the list of genes
# sharpness: Adjusts the resolution used for the sciCNV-curve calculation. Default =1.0 (range approximately 0.6-1.4).
# baseline: An optional correction to adjust the baseline CNV zero setpoint (copy number gain =0). 
#           Default = 0 and the given value for baseline ccorrecction is supposed to be
#           between -0.5 and 0.5.
# baseline_adj: The baseline adjustment is only applied to test cells if it is TRUE. Default is FALSE.

CNV_infer  <- function( ss.expr , 
                        mean.cntrl,  
                        sharpness, 
                        baseline_adj,
                        baseline 
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
  
  #---------------  T: mov av Control, before weighting & U: mov av test, before weighting
  TT <- rep(0, nrow(ss.expr)); U <- rep(0, nrow(ss.expr))
  max_1 <- pmax( as.matrix(AW - P12), matrix(-floor(P12/C297), nrow(ss.expr)) )
  min_1 <- pmin( as.matrix(P12 - BD), matrix(floor(P12/C297), nrow(ss.expr)) )
  
  for( i in 1:nrow(ss.expr) ){
    bb <- seq(i + max_1[i], i + min_1[i] , 1)
    bb <- bb[! is.na(ss.expr[bb])]
    TT[i] <- mean( ss.expr[bb] )
    U[i] <- mean( mean.cntrl[bb] )
    
  }
  
  
  #---------------------- V: ma-Wratio (from Moving Average)
  V <- (TT - U)/(TT + U + Lambda )
  
  P100 <- quantile( as.matrix(V)[seq(1,nrow(ss.expr),1)], 0.5 - BLC )
  WW <- V - P100
  
  
  #---------------  X: mov av test, before weighting. & Y: mov av control, before weighting 
  X <- rep(0, nrow(ss.expr)); Y <- rep(0, nrow(ss.expr))
  max_2 <- pmax( as.matrix(AW - P12), matrix(-floor(P12/C298), nrow(ss.expr)) )
  min_2 <- pmin( as.matrix(P12 - BD), matrix(floor(P12/C298), nrow(ss.expr)) )
  
  for(i in 1:nrow(ss.expr)){
    bb <- seq(i + max_2[i], i + min_2[i], 1)
    bb <- bb[! is.na(ss.expr[bb])]
    X[i] <- mean( ss.expr[bb] )
    Y[i] <- mean( mean.cntrl[bb] )
  }
  
  #---------------------- Z ma-Wratio  Moving Ave. +/-  &  AA Median Norm & AB: ma Weight Linirization
  Z <- (X - Y)/(X + Y + Lambda)
  
  P101 <- median( Z[seq(1,nrow(ss.expr),1) ]) 
  AA <- Z - P101
  AB <- 3.3222*(AA)^4 + 5.6399*(AA)^3 + 4.2189*(AA)^2 + 3.8956*(AA)
  
  #---------------------- ACC: ma-DMS (on post MovAv post median normalized data)  [AN in new excel file]
  
  ACC <- rep(0, nrow(ss.expr))
  
  ACC[1] <- 0
  for(i in 2:nrow(ss.expr)){
    if( WW[i] > 0){ 
      ACC[i] <-  ACC[i-1] + 1
    } else if(WW[i] < 0){ 
      ACC[i] <-  ACC[i-1] - 1
    } else
      ACC[i] <- ACC[i-1]
  }
  
  #---------------  NEW: AC: ma-DMS (on post MovAve median) [AO in our excel file]  & DA ma-Wratio  non-Moving Ave. +/-  
  AC <- ACC + ( BLC * G )
  DA <- (ss.expr - mean.cntrl) %/% (ss.expr + mean.cntrl + Lambda)
  
  
  #---------------------- AC: non-ma DMS   [DD in new excel file V28J4]
  V56 <- length( DA[ DA > 0 ] )
  V57 <- length( DA[ DA < 0 ] )
  V60 <- (V56-V57)/nrow(ss.expr)
  DD <-  rep(0, nrow(ss.expr))
  
  DD[1] <- 0
  for(i in 2:nrow(ss.expr)){
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
  
  #---------------------- AD: V7alt[ma- and non-ma] Centrograde : Slope of AC
  DI <- rep(0, nrow(ss.expr))
  
  
  for(i in 1:nrow(ss.expr)){
    if( (BD[i] < P31) & (AW[i] < P31) ){
      aaa <- seq(i - P31 ,i + P31,1)
      aaa <- aaa[! is.na(DH[aaa])]
      Ya <- as.matrix(DH[aaa])
      summary(lmResult <- lm(Ya~aaa))
      DI[i] <- as.matrix(coef(lmResult)["aaa"][[1]])
      
    } else if( (BD[i] < P31) & (AW[i] >= P31) ){
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
      
    } else if( (BD[i] >= P31) & (AW[i] < P31)  ){
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
      
    } else if( (BD[i] >= P31) & (AW[i] >= P31)  ){
      if( (P12 - BD[i]) < (AW[i] -P12) ){ 
        aaa <- seq(i + P12 - BD[i], i + AW[i] -P12 , 1)
        aaa <- aaa[! is.na(DH[aaa])]
        Ya <- as.matrix(DH[aaa])
        summary(lmResult <- lm(Ya~aaa))
        DI[i] <- as.matrix(coef(lmResult)["aaa"][[1]])
        
      } else if( (P12 - BD[i])==(AW[i] - P12) ){ 
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
  AF <- rep(0, nrow(ss.expr))
  
  for(i in 1:nrow(ss.expr)){
    if( (! is.na(DI[i]*AB[i])) & (DI[i]*AB[i] > 0) ){
      AF[i] <- DI[i]*abs(AB[i])
    } else
      AF[i] <- 0
  }
  
  AG <- AF/J311
  
  
  #### returning the infered CNV-curve matrix
  
  return(AG)
  
  
  
  
}