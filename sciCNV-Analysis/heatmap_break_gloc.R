

#################################################################################################
######                             heatmap_break_gloc function                           ####### 
###### Finding separating lines among diverse chromosomes to be used in sketching heatmap #######
######  Tiedemann Lab - Princess Margaret Cancer centre, University of Toronto            #######
######                            copyright@AliMahdipourShirayeh                          #######
#################################################################################################

heatmap_break_gloc <- function( ){

M_origin <- read.table( "./sciCNV/Dataset/10XGenomics_gen_pos_GRCh38-1.2.0.csv", sep = ',', header=TRUE)
  
## number of segments on the genome
No_Intrvl <- 1000

############
minn <- matrix(0, ncol=24, nrow=1)
maxx <- matrix(0, ncol=24, nrow=1)

for(i in 1: 22){ 
  minn[1,i] <- as.numeric(min(M_origin[which(as.matrix(M_origin[, 2]) == i) , 3]) )
  maxx[1,i] <- as.numeric(max(M_origin[which(as.matrix(M_origin[, 2]) == i) , 3]) )
}
minn[1,23] <- as.numeric(min(M_origin[which(as.matrix(M_origin[, 2]) == "X") , 3]) )
maxx[1,23] <- as.numeric(max(M_origin[which(as.matrix(M_origin[, 2]) == "X") , 3]) )

minn[1,24] <- as.numeric(min(M_origin[which(as.matrix(M_origin[, 2]) == "Y") , 3]) )
maxx[1,24] <- as.numeric(max(M_origin[which(as.matrix(M_origin[, 2]) == "Y") , 3]) )

Total_length <- sum( (maxx[1,1:24]-minn[1,1:24]) + 1 ) 
Seg_size <- ceiling( Total_length/No_Intrvl )

########

MaxNo_intervals <- ceiling(max((maxx[1,1:24]-minn[1,1:24])+ 1)/Seg_size)
Intrvls <- matrix(0, nrow = 24, ncol = MaxNo_intervals)
Min_chr <- rep(0,25)

Min_chr[1] <- 0
for(i in 1:24){
  StartEnd <- seq(minn[1,i] , maxx[1,i] , Seg_size)
  Intrvls[i, 1:length(StartEnd) ]  <-  t(as.matrix( StartEnd )) + Min_chr[i]
  Min_chr[i+1]  <-  maxx[1,i]              
}



break.gloc <- rep(0,24)
break.gloc[1] <- 1

for(i in 2: 24){ 
  break.gloc[i] <- length( Intrvls[1:(i-1),][ Intrvls[1:(i-1),] != 0])
}

return(break.gloc)

}
