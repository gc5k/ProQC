m=2000 #number of proteins observed
kappa=1 # even sample size

alpha=0.05/m
beta=0.20
SP=1000

Pvec=array(0, dim=m)
for(i in 1:m) {
  mu=runif(1, 10, 20)
  mu0=mu+rnorm(1, -0.3)
  sd=0.75

  ncl=(1+1/kappa)*(sd*(qnorm(1-alpha/2)+qnorm(1-beta))/(mu-mu0))^2
  ceiling(ncl)# 32

  mRate=(mu-10)/10
  NCP=mRate*(mu-mu0)^2/sd^2*(ncl*ncl*kappa)/(ncl+kappa*ncl)
  Power=length(which(pchisq(rchisq(SP,1,ncp=NCP),1,lower.tail=F) < alpha))
  Pvec[i]=Power
}
print(mean(Pvec/SP))
plot(sort(Pvec/SP))
