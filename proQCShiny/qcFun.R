#############function pca
library(Rcpp)
sourceCpp("cormatrix.cpp")

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

#############function pca
pEvecCat <- function(pEg, Cat, PC=c(1,2))
{
  TN=names(table(Cat))
  layout(matrix(1:length(TN), floor(sqrt(length(TN))), ceiling(sqrt(length(TN)))))
  for(i in 1:length(TN))
  {
    idx=which(Cat == TN[i])
    plot(pch=16, cex=0.7, main=TN[i], bty="l", xlab=paste0("Eigenvector ", PC[1]), ylab=paste0("Eigenvector ", PC[2]), pEg$vectors[idx,PC[1]], pEg$vectors[idx,PC[2]], xlim=range(pEg$vectors[,PC[1]])*1.1, ylim=range(pEg$vectors[,PC[2]])*1.1)
    abline(v=0, h=0, col="grey", lty=2)
    points(mean(pEg$vectors[idx,PC[1]]), mean(pEg$vectors[idx,PC[2]]), cex=1.5, lwd=2, col="green")
  }
}

pEvecCatCov <- function(pEg, Cat, Cov, m, PC=c(1,2))
{

  TN=names(table(Cat))
  layout(matrix(1:length(TN), floor(sqrt(length(TN))), ceiling(sqrt(length(TN)))))
  for(i in 1:length(TN))
  {
    idx=which(Cat == TN[i])
    colM=max(Cov/m)
    colS=min(Cov/m)
    plot(pch=16, cex=0.7, main=TN[i], bty="l", xlab=paste0("Eigenvector ", PC[1]), ylab=paste0("Eigenvector ", PC[2]), pEg$vectors[idx,PC[1]], pEg$vectors[idx,PC[2]], xlim=range(pEg$vectors[,PC[1]])*1.1, ylim=range(pEg$vectors[,PC[2]])*1.1, col=rgb((Cov/length(Cov)-colS)/colM, 1-(Cov/m-colS)/colM, 1))

    abline(v=0, h=0, col="grey", lty=2)
    points(mean(pEg$vectors[idx,PC[1]]), mean(pEg$vectors[idx,PC[2]]), cex=1.5, lwd=2, col="black")
  }
}

pEvecPlot <- function(pEg, PC=c(1,2), COL="grey", ma=0.3)
{
  pc=length(PC)
  par(mai=c(ma, ma, ma, ma))
  nt=pc*(pc-1)/2
  mat=matrix(nt+pc+1, pc, pc)
  mat[col(mat) < row(mat)] = 1:nt
  diag(mat)=(nt+1):(nt+pc)
  layout(mat)
  for(i in 1:(pc-1))
  {
    for(j in (i+1):pc)
    {
      plot(pEg$vectors[,PC[i]], pEg$vectors[, PC[j]], cex=0.5, pch=16, col=COL, xlim=1.2*range(pEg$vectors[,PC[i]]), ylim=1.2*range(pEg$vectors[,PC[j]]), bty="l", axes = F)
      abline(h=0, v=0, col="grey70", lty=2)
    }
  }
  for(i in 1:pc)
  {
    if (i== 1)
    {
      plot(x=NULL, y=NULL, xlim=1.2*range(pEg$vectors[,PC[i]]), ylim=1.2*range(pEg$vectors[,PC[i]]), bty="l")
    } else if (i == pc)
    {
      plot(x=NULL, y=NULL, xlim=1.2*range(pEg$vectors[,PC[i]]), ylim=1.2*range(pEg$vectors[,PC[i]]), bty="l")
    } else {
      plot(x=NULL, y=NULL, xlim=1.2*range(pEg$vectors[,PC[i]]), ylim=1.2*range(pEg$vectors[,PC[i]]), bty="l")
    }
    text(mean(1.1*range(pEg$vectors[,PC[i]])), mean(1.1*range(pEg$vectors[,PC[i]])), labels = paste("EV", PC[i]))
  }
}

##########proCorMatrix
proCorMatrix <- function(pMat)
{
  pMatS=matrix(0, nrow(pMat), ncol(pMat))

  for(i in 1:ncol(pMat))
  {
    pMatS[,i] = scale(pMat[,i])
  }

  exCor=matrix(0, nrow(pMat), nrow(pMat))
  datM=matrix(0, nrow=nrow(exCor) * (nrow(exCor)+1)/2, 4)
  cnt=1
  for(i in 1:nrow(pMatS))
  {
    for(j in 1:i)
    {
      idxD=which(!is.na(pMatS[i,]) & !is.na(pMatS[j,]))
      exCor[i,j]=exCor[j,i]=mean(pMatS[i,idxD] * pMatS[j,idxD])
      datM[cnt,3] = length(idxD)
      cnt=cnt+1
    }
  }

  return(exCor)
}

##pCorMatHistgram
pCorMatHist <- function(pCorMat)
{
  hist(pCorMat[row(pCorMat) < col(pCorMat)], xlab="Individual relatedness", main = "Proteomic relatedness")
}

##ProteinMissing
proMissing <- function(pMat)
{
  indMiss=array(0, ncol(pMat))
  for(i in 1:ncol(pMat))
  {
    indMiss[i]=length(which(is.na(pMat[,i])))
  }
  return(indMiss)
}

indMissing <- function(pMat)
{
  proMiss=array(0, nrow(pMat))
  for(i in 1:nrow(pMat))
  {
    proMiss[i]=length(which(is.na(pMat[i,])))
  }
  return(proMiss)

}

##
proDensityPlot <- function(pMat, TAG)
{
  layout(matrix(c(1,1,4,1,1,5,2,3,6), 3, 3, byrow = T))
  for(i in 1:nrow(pMat)) {
    dn = density(log10(as.numeric(pMat[i,!is.na(pMat[i,])])))
    if(i == 1) {
      plot(dn, main="All individuals", pch=16, cex=0.2, bty='l', xlab="")
    } else {
      lines(dn$x,dn$y, lwd=0.2)
    }
  }

  TS = names(table(TAG))
  for (i in 1:length(TS)) {
    idxC=which(TAG == TS[i])
    for (j in 1:length(idxC)) {
      dn = density(log10(as.numeric(pMat[idxC[j],!is.na(pMat[idxC[j],])])))
      if(j == 1) {
        plot(dn, main=TS[i], pch=16, cex=0.2, bty='l', xlab="")
      } else {
        lines(dn$x,dn$y, lwd=0.2)
      }
    }
  }
}

##pvca
pVCA <- function(pEg, pc=5, Cat)
{
  layout(matrix(1:2, 1, 2))
  TN=names(table(Cat))

  Wt=pEg$values[1:pc]/sum(pEg$values)
  wRsq=array(0, 3)
  idxT=which(Cat %in% TN)
  x=Cat[idxT]

  Rsq=array(0,pc)
  names(Rsq)=seq(1, pc)
  for (i in 1:pc)
  {
    an=aov(pEg$vectors[idxT,i]~x)
    Rsq[i]=summary(an)[[1]][1,2]/sum(summary(an)[[1]][,2])
  }

  barplot(main="", Rsq, ylim=c(0, max(Rsq)*1.2), ylab=expression(R^2), xlab="PC", border = F)
  barplot(Rsq*Wt, ylim=c(0, max(Rsq*Wt)*1.2), ylab=expression(paste("Weighted ", R^2, " (PVCA)")), xlab="PC", border = F)
}

###
pcaCatPlot <- function(pEg, eve = 5, Cat)
{
  TB = names(table(Cat))
  vA = matrix(0, length(TB), eve)
  rownames(vA) = TB
  for(i in 1:eve)
  {
    for(j in 1:length(TB))
    {
      vA[j,i] = var(pEg$vectors[which(Cat==TB[j]),i])
    }
  }
  barplot(t(vA), beside = T, border = F, ylab="Variation")
}
