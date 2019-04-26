Ucl=4.5
Ucs=4.7
Mcl=sort(rchisq(25000, Ucl))
Mcs=sort(rchisq(25000,Ucs))

layout(matrix(1:4, 2, 2))
hist(Mcl, breaks = 25, main="Control", col="red")
abline(v=mean(Mcl), lty=2, lwd=2)
plot(density(Mcl), main=paste0("Control (", length(which(Mcl>10)), ")"),
     col="red", bty="l", xlab="Proteome Expression (ctrl)")
abline(v=10)

hist(Mcs, breaks = 25, main="Case", col="blue")
abline(v=mean(Mcs), lty=2, lwd=2)
plot(density(Mcs), main=paste0("Case (", length(which(Mcs>10)), ")"),
    col="blue", bty="l", xlab="Proteome Expression (cs)")
abline(v=10)

REP=1000
PT=array(REP)
PTR=array(REP)
n=115
sg=0.75
for(i in 1:REP) {
  s1=rnorm(n, 14, sg)
  s2=rnorm(n, 14.5, sg)
  PT[i]=t.test(s1, s2)$p.value
  x=c(rep(1, length(s1)), rep(0, length(s2)))
  mod=lm(c(s1,s2)~x)
  PTR[i]=summary(mod)$coefficient[2,3]
}
length(which(PT < 0.05/2000))/REP
length(which(PTR < 0.05/2000))/REP

idx=which(Mcs>10)
Rpow=array(0,length(idx))
for(i in 1:length(idx)) {
  ncs=pchisq(Mcs[idx[i]], Ucs) * n
  ncl=pchisq(Mcl[idx[i]], Ucl) * n
#  ncs=n
#  ncl=n
  NCP_m=(Mcs[idx[i]]-Mcl[idx[i]])^2/sigma^2*ncs*ncl/(ncs+ncl)
  qC=rchisq(1, 1, NCP_m)
  Rpow[i]=pchisq(qC, 1)
}
length(which(1-Rpow < 0.05/1000))
