library(data.table)
library(dplyr)
library(doParallel)
library(foreach)
library(parallel)
library(ggplot2)
library(viridis)
setwd(dir="C:/Users/Duchenne/Documents/cheating/initial")


foreach(jj=1:500)%dopar%{
library(truncnorm)
library(bipartite)
library(plot3D)
library(Rmpfr)
library(circular)
library(CircStats)
library(doBy)
library(phangorn)
library(gridExtra)
library(deSolve)
library(rootSolve)
library(expm)
library(igraph)
library(reshape2)
library(ggplot2)
library(cowplot)
library(dplyr)
library(truncnorm)
#OTHER PARAMETERS
nbsp_a=20 #round(runif(1,10,30))
nbsp_p=20 #round(runif(1,10,30))

##### MORPHOLOGICAL TRAIT:
genmu_a=runif(nbsp_a,-1.5,1.5)
genmu_p=runif(nbsp_p,-1.5,1.5)
gensd_a=exp(runif(nbsp_a,log(0.1),log(1)))
gensd_p=exp(runif(nbsp_p,log(0.1),log(1)))

#### INTERACTION MATRIX
APOLL=matrix(0,nbsp_p,nbsp_a)
for(i in 1:nbsp_a){
for(i2 in 1:nbsp_p){
func <- function(x) {
 f1 <- dnorm(x,genmu_a[i],gensd_a[i])
 f2 <- dnorm(x,genmu_p[i2],gensd_p[i2])
  return(list(inte=pmin(f1,f2),poll=f1,flow=f2))
}
inte=func(seq(-10,10,0.01))
APOLL[i2,i]=sum(inte$inte)*0.01#*((1-sqrt(gensd_p[i]*gensd_f[i2]))^0.3)
}
}
IPOLL=APOLL
boxplot(IPOLL)
### COMPETITION MATRICES
CSp=matrix(runif(nbsp_p*nbsp_p,0,0),nbsp_p,nbsp_p)
diag(CSp)=1
CSa=matrix(runif(nbsp_a*nbsp_a,0,0),nbsp_a,nbsp_a)
diag(CSa)=1

#growth rates:
shape_a=runif(1,log(0.7),log(3))
shape_p=runif(1,log(0.7),log(3))
#r=c(runif(nbsp_h,-1,0),runif(nbsp_p,0,0.5))
r=c(-1*rbeta(nbsp_a,1,exp(shape_a))/2,-1*rbeta(nbsp_p,1,exp(shape_p))/2)

#save initial state
setwd(dir="C:/Users/Duchenne/Documents/cheating/initial")
save(IPOLL,APOLL,r,shape_a,shape_p,nbsp_a,nbsp_p,genmu_a,genmu_p,gensd_a,gensd_p,CSp,CSa,
file=paste0("replicate_",jj,".RData"),version = 2)
}


