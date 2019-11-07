

#######################################################################################
######                              CNV-Score function                          ####### 
######             ((  Assigning tumor/test likness score to single cells  ))   #######
######                     Tiedemann Lab, Princess Margaret Cancer centre       #######
######                     copyright @AliMahdipourShirayeh                      #######
#######################################################################################

# M_NF: is the matrix of icNV curves among with an extra column of average of 
#          iCNV-curve of all test cells

CNV_score <-  function( M_nf ){
  
  # argument validation
  if ( is.matrix(M_nf) == FALSE ){
    stop( "Please insert a proper matrix of CNV-curves")
  }
  
  #-------
  LL1 <- 0.5
  LL2 <- 1.5
  LL3 <- 2.5
  M_nf.scaled <- matrix(0, ncol=ncol(M_nf), nrow=nrow(M_nf))
  
  for (w in 1:(ncol(M_nf)-1)  ){
    for(l in 1: nrow(M_nf)){
      if( abs( M_nf[ l, w] ) <  LL1 ){
        M_nf.scaled[ l, w] <- sign(M_nf[ l, w] )*( M_nf[ l, w] )^2
      } else if( (abs(M_nf[ l, w]) >=  LL1) & ( abs(M_nf[ l, w] ) <=  LL2) ){
        M_nf.scaled[ l, w] <- sign(M_nf[ l, w] )*sqrt( abs(M_nf[ l, w] ) )
      } else if( M_nf[ l, w]  >  LL2 ){
        M_nf.scaled[ l, w] <- 2*sign( M_nf[ l, w] )*sqrt( M_nf[ l, w]  )
      }
      
    }
  }
  M_nf.scaled[ , ncol(M_nf)] <- t(as.matrix(M_nf[ , ncol(M_nf)]))

  ##################
  Score <- matrix(0, nrow = nrow(M_nf.scaled), ncol = ncol(M_nf.scaled) ) 
  TotScore <- matrix(0, nrow = 1, ncol = (ncol(M_nf.scaled)-1) ) 
  Wr <- ncol(M_nf.scaled)-1
  
  for(w in 1:(ncol(M_nf.scaled)-1)){
    for(l in 1:nrow(M_nf.scaled)){
      if(M_nf.scaled[l,ncol(M_nf.scaled)] != 0){
        Score[l,w] <- sign(M_nf.scaled[l,ncol(M_nf.scaled)])*sign(M_nf.scaled[l,w])*sqrt(abs(M_nf.scaled[l,w]))
      } else {
        Score[l,w] <- -sqrt( abs(M_nf.scaled[l,w]) )
      }
    }
  }
  
  #--------------------------------
  Score[1:nrow(M_nf.scaled), ncol(M_nf.scaled)] <- as.matrix(sign(M_nf.scaled[1:nrow(M_nf.scaled),
                                                   ncol(M_nf.scaled)])*sqrt( abs(M_nf.scaled[,ncol(M_nf.scaled)])))
  TotScore[1, ] <- 100*colSums(Score[, 1:(ncol(M_nf.scaled)-1)])/sum(Score[,ncol(M_nf.scaled)])
  colnames(TotScore) <- colnames( M_nf[, -ncol(M_nf)]) 
  
  plot(TotScore[1, ],col="blue")
  return( TotScore )
  
}
