source("../qcFun.R")
##read data
proFile="../sTable2_proteinMatrix.txt"

prot=read.table(proFile, as.is = T, header = T)

PPP0=read.csv("../PPPB1_lineup.csv", as.is = T, header = T)
PPP1=PPP0
PPP1$PPPA_ID=tolower(PPP1$PPPA_ID)

#remove missing
dat=t(prot[,c(2:ncol(prot))])
ms=array(0, dim=nrow(dat))
for(i in 1:length(ms)) {
  ms[i]=length(which(is.na(dat[i,])))
}
idxMS=which(ms>2000)

PPP1=PPP1[-idxMS,]
dat=dat[-idxMS,]
prot=prot[,-(idxMS+1)]

###quantile normalization
###https://en.wikipedia.org/wiki/Quantile_normalization
library(preprocessCore)

qdat=t(normalize.quantiles(t(dat)))
dat=qdat

colnames(dat)=prot$protein_group
rownames(dat)=rownames(t(prot[,c(2:ncol(prot))]))

#line up data
exDat=matrix(0, nrow(dat), ncol(dat))
rID=matrix(0, nrow(dat), 1)
for (i in 1:nrow(PPP1))
{
  idx=which(rownames(dat) == PPP1$PPPA_ID[i])
  rID[i,1]=idx
  exDat[i,] = dat[idx,]
}
rownames(exDat) = PPP1$PPPA_ID[rID[,1]]
colnames(exDat)=colnames(dat)

####pca analysis
outRoot="prot"
PCA(exDat, outRoot, PPP1$PPPA_ID)

####common
TSnames=c("AQUA", "CTRL", "Normal", "Tumor")
Aidx=which(PPP1$Tissue == "A")
Cidx=which(PPP1$Tissue == "C")
Nidx=which(PPP1$Tissue == "N")
Tidx=which(PPP1$Tissue == "T")
IDX=c(length(Aidx), length(Cidx), length(Nidx), length(Tidx))

TSidx=list(Aidx)
TSidx[2]=list(Cidx)
TSidx[3]=list(Nidx)
TSidx[4]=list(Tidx)

eve=read.table(paste0(outRoot, ".eigenvec"), as.is = T)

#meanplot
pdf("EigenMean.pdf")
library(gplots)
layout(matrix(1:4,2,2))
for(i in 1:4)
{
  plotmeans(eve[,2+i]~PPP1$Tissue, xlab="Group", ylab=paste("Eigenvector", i), main="Mean plot\n with 95%CI", bty='l')
  fit=aov(eve[,2+i]~PPP1$Tissue)
  TukeyHSD(fit)
}
dev.off()

#######Eigenvalue plot
grm=read.table(gzfile(paste0(outRoot,".grm.gz")), as.is = T)
eva=read.table(paste0(outRoot, ".eigenval"), as.is = T)

#jpeg("EigenVal.jpeg", width = 800, height = 600)
pdf("EigenVal.pdf")
layout(matrix(1:2, 1, 2))
hist(grm[,4], breaks=50, main = "Protein relationship matrix", xlab="Score")
barplot(eva[1:20,1], main=paste(nrow(eva), "individuals"), xlab="Eigenvalue", col=ifelse(row(eva) < 6, "red", "grey"), border = "NA")
dev.off()

######################pcA for tissue
#jpeg("PCA_tissue1.jpeg", width = 800, height = 800)
pdf("PCA_tissue1.pdf")
mat=matrix(1:4, 2, 2)

layout(mat)

for(i in 1:length(TSidx))
{
  plot(eve[TSidx[[i]],3], eve[TSidx[[i]],4], xlab="PC 1", ylab="PC 2", bty='n', axes = T, xlim=range(eve[,3], na.rm = TRUE)*1.1, ylim=1.1*range(eve[,4], na.rm = TRUE), col=i, pch=16, cex=0.5)
  points(mean(eve[TSidx[[i]],3]), mean(eve[TSidx[[i]],4]), pch=1, cex=1, lwd=3, col="grey")
  legend("bottomleft", legend = TSnames[i], bty = 'n')
}
dev.off()
