

#######################################################################################
######                             sciCNV Method                                ####### 
######   Single Cell Inferred Copy Number Variations, detected from scRNA-seq   #######
######  Tiedemann Lab - Princess Margaret Cancer centre, University of Toronto  #######
######                       copyright@AliMahdipourShirayeh                     #######
#######################################################################################

# Definitions
# norm.mat: Matrix of normalized single-cell gene expression (e.g. normalized scRNA-seq). Formatted as follows: 
#       1st column = gene name list. 
#       1st Row = cell identifier.  
#       Cells are arranged with test cells in leftward columns, followed by control cells (normal diploid cells of matching lineage) in rightward columns.
# No.test: number of test cells
# sharpness: a variable that adjusts the resolution used for the sciCNV-curve calculation (by defining a moving window size over which gene expression values are averaged). 
#       Default =1.0 (range approximately 0.6-1.4). 
#       A lower sharpness can be used to offset data sparsity and provide more reliable detection of large CNVs.
#       A higher sharpness provides more resolution for the detection of smaller CNVs but is more susceptible to noise from data sparsity and thus requires greater data density.
# baseline_adj: The baseline adjustment is only applied to test cells if TRUE is specified. Default is FALSE. 
# baseline: An optional correction to adjust the CNV zero setpoint (copy number gain =0) and improve CNV detection when CN gains and losses are substantially unbalanced.
#       Consider using for markedly hyperdiploid or hyodiploid cells with multiple trisomies or monosomies where the chromosome number deviates substantially from 46.
#       Default = 0. Ideal setting: a fraction representing the net fractional genomic change from diploidy (using a positive fraction for gain and a negative fraction for net genomic loss).
#       If the CNV profile is unknown, run the CNV analysis with default zero setting, review the preliminary CNV profile, and consider re-running with baseline ccorrection.


sciCNV <- function(norm.mat, 
                   No.test, 
                   sharpness, 
                   baseline_adj = FALSE,   # TRUE or FALSE
                   baseline = 0    # default is 0. Typical range -0.5 to +0.5.
){
  
  source(file.path(path.code, "CNV_infer.R"))
  
  ## argument validation
  if (  is.na(No.test) ){
    stop( "Number of tumor cells needs to be inserted.")
  }
  
  if (is.na(sharpness)){
    sharpness <- 1    # typical range 0.6-1.4
  }
  
  if (is.na(baseline_adj)){
    baseline_adj = FALSE    
  }
  
  if ( (baseline_adj == TRUE) & (is.na(baseline)) ){
    baseline <- 0
  }
  
  if ( (baseline_adj == FALSE) & ( baseline != 0) ){
    stop( "The baseline value is only accepted for baseline_adj = TRUE.")
  }else if ( (baseline_adj == FALSE)  ){
    baseline <- 0
  }
  
  if( (baseline > 0.5) || (baseline < -0.5) ){
    stop( "The baseline value inserted is not supported (baseline is 
          assumed to be in [-0.5, 0.5] interval).")
  }
  
  
  ## Assigning chromosome number to each gene sorted based on chromosme number, starts and ends 
  gen.Loc <- read.table( "../sciCNV/Dataset/10XGenomics_gen_pos_GRCh38-1.2.0.txt", sep = '\t', header=TRUE)
  Spec.genes <- which(as.matrix(gen.Loc[,1]) %in% as.matrix(rownames(norm.mat)))
  gen.Loc <-  gen.Loc[Spec.genes, c(1,2)]
  
  
  ## Attaching the average gene expression of the control cells 
  
  ctrl.index <- seq(No.test+1, ncol(norm.mat), 1)   # index for control cells
  
  ##
  Ave_NCs <- matrix(0, ncol=1, nrow=nrow(norm.mat))     #the average gene expression of normal control cells
  for(k in 1:nrow(norm.mat)){
    Ave_NCs[k, 1] <- mean(as.numeric(norm.mat[k, ctrl.index]) )
  }
  
  
 
  ######
  norm.mat1 <- cbind(gen.Loc, as.matrix(norm.mat), as.matrix(Ave_NCs) )
  rownames(norm.mat1) <- rownames(norm.mat)
  
  ## Focusing the sciCNV analysis on genes with expression >2% of cells (open to user specification/optimization).
  ## This enriches the analysis for informative genes.
  ## The variable excludes infrequently detected genes, that are detected in less than a minimum % of cells (assessed across both test and control cells).
  
  percent <- 0.02
  MSC1 <- norm.mat1[which(rowSums(norm.mat != 0) >= ncol(norm.mat)*percent ), ] 
  colnames(MSC1) <- c("genes","chromosomes", colnames(norm.mat),"AveNCs") 
  #chr.n <- sort(MSC1[ , 2], decreasing=FALSE)
  MSC <- MSC1[ , -c(1,2)]
  
  #####################################
  ## single cell inferred CNV function
  #####################################
  V7Alt <- matrix(0, ncol = (ncol(MSC)-1), nrow = nrow(MSC))
  
  # Average gene expression of control cells
  mean.cntrl <- as.matrix(as.numeric(MSC[, ncol(MSC)]) ) 
  
  
  
  ## sciCNV for test cells  
  V7Alt1 <- apply( as.matrix(MSC[ ,1:(No.test)]), 2, 
            function(Vector){CNV_infer(ss.expr = as.matrix(Vector), 
                                       mean.cntrl = mean.cntrl,  
                                       sharpness = sharpness, 
                                       baseline_adj = baseline_adj,
                                       baseline = baseline,
                                       new.genes = rownames(MSC1))}   )
  
  
  ## sciCNV for control cells
  V7Alt2 <- apply( as.matrix(MSC[ ,(No.test+1):(ncol(MSC)-1)]), 2, 
            function(Vector){CNV_infer(ss.expr = as.matrix(Vector), 
                                       mean.cntrl = mean.cntrl,
                                       sharpness = sharpness, 
                                       baseline_adj = baseline_adj,
                                       baseline = 0,
                                       new.genes = rownames(MSC1))}   )
  

  V7Alt <- cbind(V7Alt1, V7Alt2)
  colnames(V7Alt) <- colnames(MSC[,-ncol(MSC)])
  rownames(V7Alt) <- MSC1[ , 1]


  #### 
  
  return(V7Alt)
  
  
}


