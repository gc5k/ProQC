####pca analysis
exDat=read.table("./demo/QN_PPP1.txt", as.is = T, header = T)
PPP1=read.table("./demo/ID_PPP1.txt", as.is = T, header = T)
pCorM=proCorMatrix(exDat)
pCorMatHist(pCorM)

misPro=proMissing(exDat) #missing per protein (column)
misInd=indMissing(exDat) #missing per individual (row)

######density plot
proDensityPlot(exDat, PPP1$Tissue)

#######
pEg=eigen(pCorM)

pcaCatPlot(pEg, 5, PPP1$Tissue)
pcaCatPlot(pEg, 5, PPP1$MS_ID)

pVCA(pEg, 10, PPP1$Tissue)

pEvecPlot(pEg$vectors, c(1,2,3,4), COL=as.numeric(as.factor(PPP1$Tissue)))
pEvecCat(pEg$vectors, PPP1$Tissue, c(1,2)) #tissue
pEvecCat(pEg$vectors, PPP1$MS_ID, c(1,2)) #tissue

pEvecCatCov(pEg$vectors, PPP1$Tissue, PPP1$MS_ID, ncol(exDat), c(1,2)) #tissue

library(gplots)
layout(matrix(1:4,2,2))
for(i in 1:4)
{
  plotmeans(pEg$vectors[,i]~PPP1$Tissue, xlab="Group", ylab=paste("Eigenvector", i), main="Mean plot\n with 95%CI", bty='l')
  fit=aov(pEg$vectors[,i]~PPP1$Tissue)
  TukeyHSD(fit)
}
