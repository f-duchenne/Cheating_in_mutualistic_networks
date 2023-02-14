library(data.table)
library(deSolve)
library(rootSolve)
library(expm)
library(igraph)
library(dplyr)
library(bipartite)
setwd(dir="C:/Users/Duchenne/Documents/cheating/initial_empir/")

seuil=1e-5
conv=1e-14
handling=1
efficience=1.5
interfp=1
interff=interfp
maxiter=6000
scenarios=c("specialists","generalists")
nb_replicates=length(unique(scenarios))



# Collect command arguments
args <- commandArgs(trailingOnly = TRUE)
args_contents <- strsplit(args, ' ')
# Get first argument
jj <- as.numeric(args_contents[[1]])

set.seed(jj)
print(jj)




sites=c("AMIG","BOQU")
cost_vec=seq(0,0)
essais=1:5
simulations=c("tout","without_cheat")
tab=expand.grid(sites,cost_vec,essais,simulations)

for(jj in 1:nrow(tab)){

netcar=NULL
time1=Sys.time()
site=tab$Var1[jj]
cost=tab$Var2[jj]
essai=tab$Var3[jj]
simulation=tab$Var4[jj]

#fonctions
derivs <-function(t, y,parms){
dy=rep(0,1)
for(i in 1:length(Nini)){
if(i<=nbsp_a){
I2=T[,i]*T*y[(1+nbsp_a):(nbsp_a+nbsp_p)]
I2=matrix(I2,ncol=nbsp_a,nrow=nbsp_p)
vec=((apply(I2,2,sum))/sum(y[(1+nbsp_a):(nbsp_a+nbsp_p)]*T[,i]))
vec[which(is.na(vec))]=0
vec[which(is.infinite(vec))]=0
eq1=(r[i]+
(efficience*sum(T[,i]*y[(1+nbsp_a):(nbsp_a+nbsp_p)])-cost*sum(M[,i]*y[(1+nbsp_a):(nbsp_a+nbsp_p)]))/(1+handling*sum(T[,i]*y[(1+nbsp_a):(nbsp_a+nbsp_p)])+
interfp*sum(vec*y[1:nbsp_a]))-
y[i])*y[i]

myenv=new.env(i)
environment(eq1)=myenv
myenv$i=i
myenv$Nini=Nini
myenv$I=I
myenv$I2=I2
myenv$vec=vec
myenv$r=r
myenv$nbsp_a=nbsp_a
myenv$nbsp_p=nbsp_p
myenv$interfp=interfp
}else{
I2=t(M)[,(i-nbsp_a)]*t(M)*y[1:nbsp_a]
I2=matrix(I2,ncol=nbsp_p,nrow=nbsp_a)
vec=(apply(I2,2,sum))/sum(y[1:nbsp_a]*M[(i-nbsp_a),])
vec[which(is.na(vec))]=0
vec[which(is.infinite(vec))]=0
if(cost>0){
I2=t(C)[,(i-nbsp_a)]*t(C)*y[1:nbsp_a]
I2=matrix(I2,ncol=nbsp_p,nrow=nbsp_a)
vec2=(apply(I2,2,sum))/sum(y[1:nbsp_a]*C[(i-nbsp_a),])
vec2[which(is.na(vec2))]=0
vec2[which(is.infinite(vec2))]=0
eq1=(r[i]+
(efficience-cost)*sum(M[(i-nbsp_a),]*y[1:nbsp_a])/(1+handling*sum(M[(i-nbsp_a),]*y[1:nbsp_a])+
interff*sum(vec*y[(1+nbsp_a):(nbsp_a+nbsp_p)]))-
cost*sum(C[(i-nbsp_a),]*y[1:nbsp_a])/(1+handling*sum(C[(i-nbsp_a),]*y[1:nbsp_a])+
interff*sum(vec2*y[(1+nbsp_a):(nbsp_a+nbsp_p)]))-y[i])*y[i]
}else{
eq1=(r[i]+
(efficience-cost)*sum(M[(i-nbsp_a),]*y[1:nbsp_a])/(1+handling*sum(M[(i-nbsp_a),]*y[1:nbsp_a])+
interff*sum(vec*y[(1+nbsp_a):(nbsp_a+nbsp_p)]))-y[i])*y[i]
}

myenv=new.env(i)
environment(eq1)=myenv
myenv$i=i
myenv$I=I
myenv$Nini=Nini
myenv$I2=I2
myenv$vec=vec
myenv$r=r
myenv$nbsp_a=nbsp_a
myenv$nbsp_p=nbsp_p
myenv$interff=interff
}
if(y[i]>=seuil){dy[i]=eq1}else{dy[i]=0}}
return(list(dy))}

load(paste0(site,"_",essai,".RData"))

nbsp_a_dep=length(list_hum)
nbsp_p_dep=length(list_plant)
nbsp_a=length(list_hum)
nbsp_p=length(list_plant)

M=mut_mat
M[M>0]=1
if(simulation=="without_cheat"){C=matrix(0,nrow(M),ncol(M))}else{C=cheat_mat}
C[C>0]=1
prop_mat=cheat_mat/(cheat_mat+mut_mat)
prop_mat[is.na(prop_mat)]=0
M=M*(1-prop_mat)
C=C*prop_mat
T=M+C
Tini=T

Nini=rep(1,nbsp_a+nbsp_p)
initial=Nini

#Solve system until convergence
a=1
b=0
while(a>conv & b<maxiter){
if(a==1){out <- ode(y=Nini, times=seq(0,20,0.05), func = derivs,parms = NULL,method="lsoda")}else{
out <- ode(y=Nini, times=seq(0,20,1), func = derivs,parms = NULL,method="lsoda")}
colnames(out)[1]="Time"
Nini=t(out)[2:(1+nbsp_a+nbsp_p),nrow(out)]
Nini[which(Nini<seuil)]=0
if(is.na(mean(Nini))){a=-2}else{
if(b==0){dyna=out}else{
out[,"Time"]=out[,"Time"]+b
dyna=rbind(dyna,out[-1,])}
if(nrow(dyna)>10){a=max(apply(t(dyna[(nrow(dyna)-10):nrow(dyna),2:(1+nbsp_a+nbsp_p)]),1,var))}
if(max(Nini,na.rm=TRUE)>1e5){a=-1}
b=b+20
}
}
#matplot(dyna[,1],dyna[,-1],type="l")

#Store basic data at species level
Nini[is.na(Nini)]=0
popf_p=Nini[(nbsp_a+1):(nbsp_a+nbsp_p)]
popf_a=Nini[1:nbsp_a]
nbsp_p_per=length(popf_p[popf_p>seuil])
nbsp_a_per=length(popf_a[popf_a>seuil])

if(nbsp_p_per>0 & nbsp_a_per>0){
C=C[which(popf_p>=seuil),which(popf_a>=seuil)]
M=M[which(popf_p>=seuil),which(popf_a>=seuil)]
T=T[which(popf_p>=seuil),which(popf_a>=seuil)]
r=r[which(Nini>=seuil)]
popf_a=popf_a[which(popf_a>=seuil)]
popf_p=popf_p[which(popf_p>=seuil)]
nbsp_p=length(popf_p)
nbsp_a=length(popf_a)
Nini=c(popf_a,popf_p)

T=matrix(T,ncol=nbsp_a_per,nrow=nbsp_p_per)
M=matrix(M,ncol=nbsp_a_per,nrow=nbsp_p_per)
C=matrix(C,ncol=nbsp_a_per,nrow=nbsp_p_per)

jacob=rootSolve::jacobian.full(Nini,derivs,pert = 1e-6)
eig=eigen(jacob,only.values=TRUE)$values

netcar=rbind(netcar,data.frame(interf=interff,essai=essai,simulation=simulation,site=site,cost=cost,nbsp_p_dep=nbsp_p_dep,nbsp_a_dep=nbsp_a_dep,nbsp_a=nbsp_a_per,nbsp_p=nbsp_p_per,
pers_tot=(nbsp_a_per+nbsp_p_per)/(nbsp_a_dep+nbsp_p_dep),
variance=a,valprop=max(Re(eig))))
}else{
netcar=rbind(netcar,data.frame(interf=interff,essai=essai,simulation=simulation,site=site,cost=cost,nbsp_p_dep=nbsp_p_dep,nbsp_a_dep=nbsp_a_dep,nbsp_a=nbsp_a_per,nbsp_p=nbsp_p_per,
pers_tot=(nbsp_a_per+nbsp_p_per)/(nbsp_a_dep+nbsp_p_dep),
variance=a,valprop=NA))
}



fwrite(netcar,paste0("C:/Users/Duchenne/Documents/cheating/eq_empir/netcar_",jj,".txt"),row.names=FALSE,sep="\t")
}
