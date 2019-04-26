rep=1000
pvec=array(0, dim=c(3,rep))
ncs=100
ncl=125
m=1
mu1=13
mu2=13.5
sigma=0.75
for(i in 1:rep) {
  z1=rnorm(ncs, mu1, sigma)
  z2=rnorm(ncl, mu2, sigma)
  tt=t.test(z1, z2)
  pvec[1,i]=tt$p.value

  Z=c(z1, z2)
  B=c(rep(1, ncs), rep(0, ncl))
  mod=lm(Z~B)
  pvec[2,i]=summary(mod)$coefficient[2,4]
  pvec[3,i]=summary(mod)$coefficient[2,3]^2
}
length(which(pvec[1,] < 0.025/m))
length(which(pvec[2,] < 0.025/m))
length(which(pvec[3,] > qchisq(0.025/m, 1, lower.tail = F)))

NCP=(mu1-mu2)^2/sigma^2*ncs*ncl/(ncs+ncl)
qqplot(pvec[3,], rchisq(rep, 1, ncp=NCP), pch=16)
abline(a=0, b=1)


