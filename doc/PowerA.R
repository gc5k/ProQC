rep=1000
pvec=array(0, rep)
for(i in 1:rep) {
  z1=rnorm(125, 13, 0.75)
  z2=rnorm(125, 13.5, 0.75)  
  tt=t.test(z1, z2)
  pvec[i]=tt$p.value
}
length(which(pvec < 0.025/2000))
