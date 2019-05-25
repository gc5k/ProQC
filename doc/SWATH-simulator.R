Ucs=2.5
Ucl=2.3
Ukeep=2.3

n_cs=200
n_cl=200

wHouse=0.6
csDriver=0.8
M=30000
m=2000

ProExp=array(0, dim=c(3, M))

for(i in 1:M) {
  isHouseKeep=runif(1)
  if (isHouseKeep < wHouse) {
    ProExp[1, i] = ProExp[2, i] = qchisq(runif(1, 0, 1), Ukeep)
    ProExp[3, i] = 3
  } else {
    Q=runif(1, 0, 1)
    ep1=qchisq(Q, Ucs)
    ep2=qchisq(Q, Ucl)
    if (runif(1, 0, 1) < csDriver) {
      ProExp[1, i] = ep1
      ProExp[2, i] = ep2
      ProExp[3, i] = 1
    } else {
      ProExp[1, i] = ep2
      ProExp[2, i] = ep1
      ProExp[3, i] = 2
    }
  }
}

ExpCS=matrix(0, n_cs, M)
ExpCL=matrix(0, n_cl, M)
for(i in 1:n_cs) {
  ExpCS[i,] = rnorm(M, ProExp[1,], 1)
  idx=which(ExpCS[i,] < 10)
  ExpCS[i,idx] = 0
}

for(i in 1:n_cl) {
  ExpCL[i,] = rnorm(M, ProExp[2,], 1)
  idx=which(ExpCL[i,] < 10)
  ExpCL[i,idx] = 0
}

Pmat=rbind(ExpCS, ExpCL)
idx=which(colMeans(Pmat) != 0)
PmatFinal=Pmat[,idx]
od=order(colMeans(PmatFinal))
PmatFinal=PmatFinal[,od]

ProSum=array(0, dim=c(4, ncol(PmatFinal)))
for(i in 1:ncol(ProSum)) {
  ProSum[1,i] = mean(PmatFinal[which(PmatFinal[1:n_cs, i]!=0), i])
  ProSum[2,i] = length(which(PmatFinal[1:n_cs,i]!=0))
  ProSum[3,i] = mean(PmatFinal[n_cs+which(PmatFinal[(1+n_cs):(n_cs+n_cl),i]!=0), i])
  ProSum[4,i] = length(which(PmatFinal[(1+n_cs):(n_cs+n_cl),i]!=0))
}
layout(matrix(1:2, 2, 1))
plot(ProSum[1,], ProSum[2,], pch=16, xlab="Expression", ylab="Cnt", bty="L")
plot(ProSum[3,], ProSum[4,], pch=16, xlab="Expression", ylab="Cnt", bty="L")
