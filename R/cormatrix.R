#############function pca
library(Rcpp)
sourceCpp("./cormatrix.cpp")

##########proCorMatrix by Rcpp
proCorMatrix_c <- function(pMat)
{
  pMatS=matrix(0, nrow(pMat), ncol(pMat))

  for(i in 1:ncol(pMat))
  {
    pMatS[,i] = scale(pMat[,i])
  }
  
  pmats0=pMatS
  pmats0[which(is.na(pmats0))]=0

  exCor=CorMatrix(pmats0)
  return(exCor)
}

