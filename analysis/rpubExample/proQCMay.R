##read data
proFile="2018-05-21pppb_prot.txt"
prot=read.table(proFile, as.is = T, header = T)
rownames(prot)=prot[,1]
prot=prot[,-c(1,2)]
PPP0=read.csv("PPPB1_lineup.csv", as.is = T, header = T)
PPP1=PPP0
PPP1$PPPA_ID=tolower(PPP1$PPPA_ID)

datmp=t(prot)

###line up protein matrix and datmp
exDat=matrix(0, nrow(datmp), ncol(datmp))
rID=matrix(0, nrow(datmp), 1)
for (i in 1:nrow(PPP1)) {
  idx=which(rownames(datmp) == PPP1$PPPA_ID[i])
  rID[i,1]=idx
  exDat[i,] = datmp[idx,]
}
rownames(exDat) = PPP1$PPPA_ID[rID[,1]]
colnames(exDat) = rownames(prot)

#evaluate missing value
misD=matrix(1, nrow(exDat), 4)
for(i in 1:nrow(misD))
{
  naIdx=which(is.na(exDat[i,]))
  if (length(naIdx) < ncol(exDat) )
  {
    misD[i, 1] = length(naIdx)/ncol(exDat)
    misD[i, 2] = min(datmp[i,-naIdx])
    misD[i, 3] = median(exDat[i, -naIdx])
    misD[i, 4] = mean(exDat[i, -naIdx])
  }
}
blankInd = which(misD[,1]> 0.8)

hist(misD[,1], breaks = 50, main="Missing per individual")
if(length(blankInd)> 0)
{
  exDat=exDat[-blankInd,]
  PPP1=PPP1[-blankInd,]
}

###quantile normalization
###https://en.wikipedia.org/wiki/Quantile_normalization
library(preprocessCore)

#exDat1=exDat
#qdat=t(normalize.quantiles(t(exDat1)))

qdat=exDat
rownames(qdat)=rownames(exDat)
colnames(qdat)=colnames(exDat)

for(i in 1:ncol(datmp))
{
  idx=which(rownames(qdat) == rownames(datmp)[i])
  if(length(idx) > 0)
  {
    datmp[i,] = qdat[idx,]
  }
}

write.table(t(datmp), "Apr_PPP1.txt", row.names = T, col.names = T, quote=F)


###proQC package
qdat=t(read.table("Apr_PPP1.txt", as.is = T, header = T)) #row for individual, column for protein

PPP1=read.csv("PPPB1_lineup.csv", as.is = T, header = T)
indM=array(0, nrow(qdat))
for(i in 1:nrow(qdat)) {
  indM[i]=length(which(is.na(qdat[i,])))
}
idx=which(indM/ncol(qdat) > 0.8)

qdat=qdat[-idx,]
PPP1=PPP1[-idx,]
indM1=indM[-idx]

TAG=c("A", "C", "N", "T")

layout(matrix(c(1,1,2,1,1,3,4,5,6),3,3, byrow=F))
hist(indM1/ncol(qdat), main="All samples", xlab="Missing Rate")
for(i in 1:length(TAG)) {
  hist(indM1[PPP1$Tissue==TAG[i]]/ncol(qdat), main=TAG[i], xlim=c(0, 1), xlab="Missing Rate")
}

indMPro=list(qdat[PPP1$Tissue==TAG[1],], qdat[PPP1$Tissue==TAG[2],],
             qdat[PPP1$Tissue==TAG[3],],qdat[PPP1$Tissue==TAG[4],])
layout(matrix(1:4, 1,4))
for(k in 1:4) {
  refAmean=colMeans(indMPro[[k]], na.rm = T)
  refAmiss=array(1, dim=length(refAmean))
  for(i in 1:length(refAmean)) {
    if(length(which(!is.na(indMPro[[k]][,i])))>0) {
      refAmiss[i] = length(which(!is.na(indMPro[[k]][,i])))
    }
  }
  plot(main=TAG[k], xlab="Missing rate", ylab="Expression", xlim=c(0, 1), ylim=c(0, 23), 1-refAmiss/nrow(indMPro[[k]]), refAmean,cex=0.5,pch=16, col="grey")
  xx=seq(0, 1, length.out = ncol(indMPro[[k]]))
  idNA=which(is.na(refAmean))

  ss=smooth.spline(1-refAmiss[-idNA]/nrow(indMPro[[k]]), refAmean[-idNA], df=5)
  if(k==1) {
    ssA=smooth.spline(1-refAmiss[-idNA]/nrow(indMPro[[k]]), refAmean[-idNA], df=5)
  }
  lines(ssA, col="blue", lty=2, lwd=2)
  lines(ss, col=ifelse(k!=1, "red", "blue"), lty=1, lwd=2)
}

pmat=matrix(NA, 6, ncol(indMPro[[1]]))
pIdx=1
layout(matrix(1:6, 2, 3))
for(i in 1:3) {
  dn=indMPro[[i]]
  for(j in (i+1):4) {
    dm=indMPro[[j]]
    for(k in 1:ncol(dn)) {
      if( (length(which(!is.na(dn[,k]))) > 1) & (length(which(!is.na(dm[,k]))) >1)) {
        ft= var.test(dn[,k], dm[,k], na.rm = T)
        pmat[pIdx,k] = ft$p.value
      }
    }
    pIdx=pIdx+1
  }
}

pIdx=1
layout(matrix(c(1,2,3,8,4,5,7,7,6), 3, 3))
gc=array(0,6)
for(i in 1:3) {
  for(j in (i+1):4) {
    hist(main=paste0(TAG[j], "/", TAG[i]), xlab="p-value (F test)", 
         pmat[pIdx,which(!is.na(pmat[pIdx,]))], 
         breaks = 20, col=pIdx)
    abline(h=length(which(!is.na(pmat[pIdx,])))/20, lty=2, col="grey")
    gc[pIdx]=qchisq(median(pmat[pIdx,], na.rm = T), df=1, lower.tail = F)/0.455
    legend("topright", bty = 'n', legend = paste("GC=", format(gc[pIdx], digits = 3)))
    pIdx=pIdx+1
  }
}

pIdx=1
for(i in 1:3) {
  for(j in (i+1):4) {
    if(pIdx==1) {
      qqplot(col=pIdx, ylim=c(0,70), xlab="Expected -log10(p)", ylab="Observed -log10(p)",
             -log10(runif(length(which(pmat[pIdx,]<1)))),-log10(pmat[pIdx,]), pch=16, cex=0.5)
      abline(a=0, b=1, col="grey", lty=2)
    } else {
      points(sort(-log10(runif(length(which(pmat[pIdx,]<1))))),sort(-log10(pmat[pIdx,])), col=pIdx, pch=16,cex=0.5)
    }
    pIdx=pIdx+1
  }
}
legend("topleft", bty='n', legend = paste0("GC=", format(gc, digits=3)), pch=15, col=1:6)
#
idxA=which(PPP1$Tissue == "A")
avA=aov(indM1[idxA]/ncol(qdat) ~ PPP1$MS_ID[idxA])
idxC=which(PPP1$Tissue == "C")
avC=aov(indM1[idxC]/ncol(qdat) ~ PPP1$MS_ID[idxC])

idxTN=which(PPP1$Tissue == "T" | PPP1$Tissue == "N")
avTN=aov(indM1[idxTN]/ncol(qdat) ~ PPP1$MS_ID[idxTN] + PPP1$Tissue[idxTN])

pmat=matrix(NA, nrow = ncol(qdat), 6)
qdatN=qdat[PPP1$Tissue == "N",]
qdatT=qdat[PPP1$Tissue == "T",]
pmat[,1] = proMissing(qdatN)/nrow(qdatN)
pmat[,2] = proMissing(qdatT)/nrow(qdatT)
plot(pmat[,1], pmat[,2], bty="l", xlab="Norm missing rate", ylab="Tumor missing rate", pch=16, cex=0.5)
text(0.2,0.8, labels = format(cor(pmat[,1], pmat[,2]), digits = 4))
abline(a=0, b=1, col="red", lty=2, lwd=3)
for(i in 1:ncol(qdat)) {
  if(pmat[i,1]< 0.95 && pmat[i,2] < 0.95) {
    pmat[i,3] = t.test(qdatN[,i], qdatT[,i])$p.value
    pmat[i,4] = t.test(qdatN[,i], qdatT[,i])$statistic
    pmat[i,5] = median(qdatN[,i], na.rm = T)
    pmat[i,6] = median(qdatT[,i], na.rm = T)
  }
}

#chiT=(pmat[,1]-pmat[,2])^2/(pmat[,1]*(1-pmat[,1])+pmat[,2]*pmat[,2])
layout(matrix(1:2, 1, 2))
Cat1=which(pmat[,4]<0)
plot(ylim=c(0, 150), main=paste0("Tumor High ", length(Cat1)), pmat[Cat1,1]-pmat[Cat1,2], -log10(pmat[Cat1,3]), xlab="Missing(N)-Missing(T)", ylab="p-value for expression differentiation", col=ifelse(-log10(pmat[Cat1,3]) > -log10(0.05/length(Cat1)), "red", "grey"), pch=16, cex=0.4, bty='l')
abline(h=c(40,60,80), col="blue", lty=2)

plot(ylim=c(0, 150), main=paste0("Tumor Low ", ncol(qdat)-length(Cat1)), pmat[-Cat1,1]-pmat[-Cat1,2], -log10(pmat[-Cat1,3]), xlab="Missing(N)-Missing(T)", ylab="p-value for expression differentiation", col=ifelse(-log10(pmat[-Cat1,3]) > -log10(0.05/(ncol(qdat)-length(Cat1))), "red", "grey"), pch=16, cex=0.4, bty='l')
abline(h=c(40,60,80), col="blue", lty=2)

layout(matrix(1:2, 1, 2))
plot(pmat[,5], pmat[,1], pch=16, cex=0.3)
plot(pmat[,6], pmat[,2], pch=16, cex=0.3)

#generate correlation matrix
pCorM=proCorMatrix(qdat)

######density plot
proDensityPlot(qdat, PPP1$Tissue)

misInd=proMissing(qdat)
#######
pEg=eigen(pCorM)
barplot(pEg$values[1:20], border = F)

pEvecPlot(pEg, c(1,2,3,4), COL=as.numeric(as.factor(PPP1$Tissue)))
pEvecCat(pEg, PPP1$Tissue, c(1,2)) #tissue
pEvecCat(pEg, PPP1$MS_ID, c(1,2)) #tissue

pcaCatPlot(pEg, 5, PPP1$Tissue)
pcaCatPlot(pEg, 5, PPP1$MS_ID)

pVCA(pEg, 10, PPP1$Tissue)
pVCA(pEg, 10, PPP1$MS_ID)


indMissing <- function(pMat)
{
  proMiss=array(0, nrow(pMat))
  for(i in 1:nrow(pMat))
  {
    proMiss[i]=length(which(is.na(pMat[i,])))
  }
  return(proMiss)
  
}