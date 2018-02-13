#############function pca
pEvecCat <- function(pEvec, cat, PC=c(1,2))
{
  tName=names(table(cat))
  layout(matrix(1:length(tName), floor(sqrt(length(tName))), ceiling(sqrt(length(tName)))))
  for(i in 1:length(tName))
  {
    idx=which(cat == tName[i])
    plot(pch=16, cex=0.7, main=tName[i], bty="l", xlab=paste0("Eigenvector ", PC[1]), ylab=paste0("Eigenvector ", PC[2]), pEvec[idx,PC[1]], pEvec[idx,PC[2]], xlim=range(pEvec[,PC[1]])*1.1, ylim=range(pEvec[,PC[2]])*1.1)
    abline(v=0, h=0, col="grey", lty=2)
    points(mean(pEvec[idx,PC[1]]), mean(pEvec[idx,PC[2]]), cex=1.5, lwd=2, col="green")
  }
}

pEvecCatCov <- function(pEvec, cat, cov, m,PC=c(1,2))
{

  tName=names(table(cat))
  layout(matrix(1:length(tName), floor(sqrt(length(tName))), ceiling(sqrt(length(tName)))))
  for(i in 1:length(tName))
  {
    idx=which(cat == tName[i])
    colM=max(MisP/m)
    colS=min(MisP/m)
    plot(pch=16, cex=0.7, main=tName[i], bty="l", xlab=paste0("Eigenvector ", PC[1]), ylab=paste0("Eigenvector ", PC[2]), pEvec[idx,PC[1]], pEvec[idx,PC[2]], xlim=range(pEvec[,PC[1]])*1.1, ylim=range(pEvec[,PC[2]])*1.1, col=rgb((MisP/ncol(exDat)-colS)/colM, 1-(MisP/m-colS)/colM, 1))# (MisP/m-colS)/colM))

    abline(v=0, h=0, col="grey", lty=2)
    points(mean(pEvec[idx,PC[1]]), mean(pEvec[idx,PC[2]]), cex=1.5, lwd=2, col="black")
  }
}

pEvecPlot <- function(pEvec, PC=c(1,2), COL="grey", ma=0.3)
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
      plot(pEvec[,PC[i]+2],pEvec[, PC[j]+2], cex=0.5, pch=16, col=COL, xlim=1.2*range(pEvec[,PC[i]+2]), ylim=1.2*range(pEvec[,PC[j]+2]), bty="l", axes = F)
      abline(h=0, v=0, col="grey70", lty=2)
    }
  }
  for(i in 1:pc)
  {
    if (i== 1)
    {
      plot(x=NULL, y=NULL, xlim=1.2*range(pEvec[,PC[i]+2]), ylim=1.2*range(pEvec[,PC[i]+2]), bty="l")
    } else if (i == pc)
    {
      plot(x=NULL, y=NULL, xlim=1.2*range(pEvec[,PC[i]+2]), ylim=1.2*range(pEvec[,PC[i]+2]), bty="l")
    } else {
      plot(x=NULL, y=NULL, xlim=1.2*range(pEvec[,PC[i]+2]), ylim=1.2*range(pEvec[,PC[i]+2]), bty="l")
    }
    text(mean(1.2*range(pEvec[,PC[i]+2])), mean(1.2*range(pEvec[,PC[i]+2])),labels = paste("EV", PC[i]))
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
