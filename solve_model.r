library(data.table)
library(deSolve)
library(rootSolve)
library(expm)
library(igraph)
library(dplyr)
library(bipartite)
setwd(dir="/home/duchenne/cheating/initial/")

seuil=1e-5
conv=1e-13
handling=1
efficience=1.4
interfp=1
interff=interfp
maxiter=3000
scenarios=c("specialists","generalists")
nb_replicates=length(unique(scenarios))

# Collect command arguments
args <- commandArgs(trailingOnly = TRUE)
args_contents <- strsplit(args, ' ')
# Get first argument
jj <- as.numeric(args_contents[[1]])

set.seed(jj)

netcar=NULL
prop_cheaters_vec=seq(0.1,0.6,0.1)
prop_cheating_vec=seq(0.1,1,0.1)
prop_innovative_vec=seq(0,1,0.2)
cost_vec=seq(0,0.6,0.3)
connectance_vec=seq(0.2,0.4,0.2)
tab=rbind(data.frame(Var1=0,Var2=0,Var3=0,Var4=0,Var5=unique(connectance_vec)),expand.grid(prop_cheaters_vec,prop_cheating_vec,prop_innovative_vec,cost_vec,connectance_vec))
for(jjj in 1:nrow(tab)){
time1=Sys.time()
prop_cheaters=tab$Var1[jjj]
prop_cheating=tab$Var2[jjj]
prop_innovative=tab$Var3[jjj]
cost=tab$Var4[jjj]
connectance=tab$Var5[jjj]
if(prop_cheaters>0){nbrep=nb_replicates}else{nbrep=1}

for(random in 1:nbrep){

scenario=scenarios[random]

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
efficience*sum(T[,i]*y[(1+nbsp_a):(nbsp_a+nbsp_p)])/(1+handling*sum(T[,i]*y[(1+nbsp_a):(nbsp_a+nbsp_p)])+
interfp*sum(vec*y[1:nbsp_a]))-
sum(CSa[,i]*y[1:nbsp_a]))*y[i]

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
efficience*sum(M[(i-nbsp_a),]*y[1:nbsp_a])/(1+handling*sum(M[(i-nbsp_a),]*y[1:nbsp_a])+
interff*sum(vec*y[(1+nbsp_a):(nbsp_a+nbsp_p)]))-
cost*sum(C[(i-nbsp_a),]*y[1:nbsp_a])/(1+handling*sum(C[(i-nbsp_a),]*y[1:nbsp_a])+
interff*sum(vec2*y[(1+nbsp_a):(nbsp_a+nbsp_p)]))-
sum(CSp[,(i-nbsp_a)]*y[(1+nbsp_a):(nbsp_a+nbsp_p)]))*
y[i]
}else{
eq1=(r[i]+
efficience*sum(M[(i-nbsp_a),]*y[1:nbsp_a])/(1+handling*sum(M[(i-nbsp_a),]*y[1:nbsp_a])+
interff*sum(vec*y[(1+nbsp_a):(nbsp_a+nbsp_p)]))-
sum(CSp[,(i-nbsp_a)]*y[(1+nbsp_a):(nbsp_a+nbsp_p)]))*
y[i]
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

load(paste0("replicate_",jj,".RData"))
IPOLL[IPOLL<quantile(c(IPOLL),probs=(1-connectance))]=0
IPOLL[IPOLL>0]=1
nbsp_a_dep=nbsp_a
nbsp_p_dep=nbsp_p
Iini=IPOLL
IM=IPOLL
IO=1-IM

r[r>(-0.001)]=-0.001

vec_sha=apply(IPOLL,2,vegan::diversity)
cheaters=rep(0,nbsp_a)
if(scenario=="generalists"){cheaters[which(rank(vec_sha,ties.method="first")>(nbsp_a-round(prop_cheaters*nbsp_a)))]=1}
if(scenario=="specialists"){cheaters[which(rank(vec_sha,ties.method="first")<=round(prop_cheaters*nbsp_a))]=1}
cheaters_ini=cheaters

cheating_mat=do.call(rbind, replicate(nbsp_p,cheaters*prop_cheating,simplify=F))
M=IM*(1-cheating_mat)
Mini=M
C=cheating_mat*((1-prop_innovative)*IM+prop_innovative*IO)
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
IM=IM[which(popf_p>=seuil),which(popf_a>=seuil)]
IO=IO[which(popf_p>=seuil),which(popf_a>=seuil)]
CSp=CSp[which(popf_p>=seuil),which(popf_p>=seuil)]
CSa=CSa[which(popf_a>=seuil),which(popf_a>=seuil)]
r=r[which(Nini>=seuil)]
cheaters=cheaters[which(popf_a>=seuil)]
popf_a=popf_a[which(popf_a>=seuil)]
popf_p=popf_p[which(popf_p>=seuil)]
mati=sqrt(popf_p%*%t(popf_a))
resf=mati*IM
nbsp_p=length(popf_p)
nbsp_a=length(popf_a)
Nini=c(popf_a,popf_p)

T=matrix(T,ncol=nbsp_a_per,nrow=nbsp_p_per)
IM=matrix(IM,ncol=nbsp_a_per,nrow=nbsp_p_per)
IO=matrix(IO,ncol=nbsp_a_per,nrow=nbsp_p_per)
T=matrix(T,ncol=nbsp_a_per,nrow=nbsp_p_per)
M=matrix(M,ncol=nbsp_a_per,nrow=nbsp_p_per)
C=matrix(C,ncol=nbsp_a_per,nrow=nbsp_p_per)
CSa=matrix(CSa,ncol=nbsp_a_per,nrow=nbsp_a_per)
CSp=matrix(CSp,ncol=nbsp_p_per,nrow=nbsp_p_per)

jacob=rootSolve::jacobian.full(Nini,derivs,pert = 1e-6)
eig=eigen(jacob,only.values=T)$values

inv=-1*solve(jacob)
jacob2=jacob
inv2=inv
for(zz in 1:nrow(jacob)){
for(zzz in 1:ncol(jacob)){
inv2[zz,zzz]=inv[zz,zzz]/(inv[zz,zz]*inv[zzz,zzz]-inv[zz,zzz]*inv[zzz,zz])
}}
diag(jacob2)=NA
diag(inv2)=NA

ai_direct_aa=mean(jacob2[(1:nbsp_a),(1:nbsp_a)],na.rm=TRUE)
ai_direct_pp=mean(jacob2[(nbsp_a+1):(nbsp_a+nbsp_p),(nbsp_a+1):(nbsp_a+nbsp_p)],na.rm=TRUE)
ai_direct_pa=mean(jacob2[(1:nbsp_a),(nbsp_a+1):(nbsp_a+nbsp_p)],na.rm=TRUE)
ai_direct_ap=mean(jacob2[(nbsp_a+1):(nbsp_a+nbsp_p),(1:nbsp_a)],na.rm=TRUE)
ai_net_ap=mean(inv2[(nbsp_a+1):(nbsp_a+nbsp_p),(1:nbsp_a)],na.rm=TRUE)
ai_net_aa=mean(inv2[(1:nbsp_a),(1:nbsp_a)],na.rm=TRUE)
ai_net_pp=mean(inv2[(nbsp_a+1):(nbsp_a+nbsp_p),(nbsp_a+1):(nbsp_a+nbsp_p)],na.rm=TRUE)
ai_net_pa=mean(inv2[(1:nbsp_a),(nbsp_a+1):(nbsp_a+nbsp_p)],na.rm=TRUE)
ai_indirect_pa=mean(inv2[(1:nbsp_a),(nbsp_a+1):(nbsp_a+nbsp_p)]-jacob2[(1:nbsp_a),(nbsp_a+1):(nbsp_a+nbsp_p)],na.rm=TRUE)
ai_indirect_ap=mean(inv2[(nbsp_a+1):(nbsp_a+nbsp_p),(1:nbsp_a)]-jacob2[(nbsp_a+1):(nbsp_a+nbsp_p),(1:nbsp_a)],na.rm=TRUE)
ai_indirect_aa=mean(inv2[(1:nbsp_a),(1:nbsp_a)]-jacob2[(1:nbsp_a),(1:nbsp_a)],na.rm=TRUE)
ai_indirect_pp=mean(inv2[(nbsp_a+1):(nbsp_a+nbsp_p),(nbsp_a+1):(nbsp_a+nbsp_p)]-jacob2[(nbsp_a+1):(nbsp_a+nbsp_p),(nbsp_a+1):(nbsp_a+nbsp_p)],na.rm=TRUE)
ai_contrib_pa=mean((inv2[(1:nbsp_a),(nbsp_a+1):(nbsp_a+nbsp_p)]-jacob2[(1:nbsp_a),(nbsp_a+1):(nbsp_a+nbsp_p)])/(abs(inv2[(1:nbsp_a),(nbsp_a+1):(nbsp_a+nbsp_p)]-jacob2[(1:nbsp_a),(nbsp_a+1):(nbsp_a+nbsp_p)])+abs(jacob2[(1:nbsp_a),(nbsp_a+1):(nbsp_a+nbsp_p)])),na.rm=TRUE)
ai_contrib_ap=mean((inv2[(nbsp_a+1):(nbsp_a+nbsp_p),(1:nbsp_a)]-jacob2[(nbsp_a+1):(nbsp_a+nbsp_p),(1:nbsp_a)])/(abs(inv2[(nbsp_a+1):(nbsp_a+nbsp_p),(1:nbsp_a)]-jacob2[(nbsp_a+1):(nbsp_a+nbsp_p),(1:nbsp_a)])+abs(jacob2[(nbsp_a+1):(nbsp_a+nbsp_p),(1:nbsp_a)])),na.rm=TRUE)
ai_contrib_aa=mean((inv2[(1:nbsp_a),(1:nbsp_a)]-jacob2[(1:nbsp_a),(1:nbsp_a)])/(abs(inv2[(1:nbsp_a),(1:nbsp_a)]-jacob2[(1:nbsp_a),(1:nbsp_a)])+abs(jacob2[(1:nbsp_a),(1:nbsp_a)])),na.rm=TRUE)
ai_contrib_pp=mean((inv2[(nbsp_a+1):(nbsp_a+nbsp_p),(nbsp_a+1):(nbsp_a+nbsp_p)]-jacob2[(nbsp_a+1):(nbsp_a+nbsp_p),(nbsp_a+1):(nbsp_a+nbsp_p)])/(abs(inv2[(nbsp_a+1):(nbsp_a+nbsp_p),(nbsp_a+1):(nbsp_a+nbsp_p)]-jacob2[(nbsp_a+1):(nbsp_a+nbsp_p),(nbsp_a+1):(nbsp_a+nbsp_p)])+abs(jacob2[(nbsp_a+1):(nbsp_a+nbsp_p),(nbsp_a+1):(nbsp_a+nbsp_p)])),na.rm=TRUE)

if(length(cheaters[cheaters==1])>0){
cheaters_to_plant=mean(inv2[(nbsp_a+1):(nbsp_a+nbsp_p),cheaters==1])
cheaters_to_poll=mean(inv2[1:nbsp_a,cheaters==1],na.rm=TRUE)
}else{
cheaters_to_plant=NA
cheaters_to_poll=NA
}

if(length(cheaters[cheaters==0])>0){
no_cheaters_to_plant=mean(inv2[(nbsp_a+1):(nbsp_a+nbsp_p),cheaters==0])
no_cheaters_to_poll=mean(inv2[1:nbsp_a,cheaters==0],na.rm=TRUE)
}else{
no_cheaters_to_plant=NA
no_cheaters_to_poll=NA
}



netcar=rbind(netcar,data.frame(interf=interff,essai=jj,scenario=scenario,prop_cheaters=prop_cheaters,prop_cheating=prop_cheating,prop_innovative=prop_innovative,cost=cost,connectance,
nbsp_p_dep=nbsp_p_dep,nbsp_a_dep=nbsp_a_dep,nbsp_a=nbsp_a_per,nbsp_p=nbsp_p_per,
pers_tot=(nbsp_a_per+nbsp_p_per)/(nbsp_a_dep+nbsp_p_dep),generalism_cheaters=if(length(cheaters_ini[cheaters_ini==1])>0){mean(apply(Iini[,cheaters_ini==1],2,mean))}else{NA},generalism=mean(apply(Iini,2,mean)),
variance=a,valprop=max(Re(eig)),
ai_direct_aa=ai_direct_aa,ai_direct_pp=ai_direct_pp,ai_direct_pa=ai_direct_pa,
ai_direct_ap=ai_direct_ap,ai_net_ap=ai_net_ap,ai_net_aa=ai_net_aa,ai_net_pp=ai_net_pp,ai_net_pa=ai_net_pa,
ai_indirect_aa=ai_indirect_aa,ai_indirect_pp=ai_indirect_pp,ai_indirect_pa=ai_indirect_pa,ai_indirect_ap=ai_indirect_ap,
ai_contrib_ap=ai_contrib_ap,ai_contrib_aa=ai_contrib_aa,ai_contrib_pp=ai_contrib_pp,ai_contrib_pa=ai_contrib_pa,shape_a=exp(shape_a),shape_p=exp(shape_p),nb_cheaters_dep=length(cheaters_ini[cheaters_ini==1]),nb_cheaters=length(cheaters[cheaters==1]),
cheaters_to_plant=cheaters_to_plant,cheaters_to_poll=cheaters_to_poll,no_cheaters_to_plant=no_cheaters_to_plant,no_cheaters_to_poll=no_cheaters_to_poll,
NODF_Mini=networklevel(Mini,index="weighted NODF")[[1]],NODF_M=if(nbsp_a>1 & nbsp_p>1){networklevel(M,index="weighted NODF")[[1]]}else{NA},NODF_Tini=networklevel(Tini,index="weighted NODF")[[1]],
NODF_T=if(nbsp_a>1 & nbsp_p>1){networklevel(T,index="weighted NODF")[[1]]}else{NA},
NODF_Iini=networklevel(Iini,index="weighted nestedness")[[1]],NODF_IM=if(nbsp_a>1 & nbsp_p>1){networklevel(IM,index="weighted nestedness")[[1]]}else{NA}))
}else{
netcar=rbind(netcar,data.frame(interf=interff,essai=jj,scenario=scenario,prop_cheaters=prop_cheaters,prop_cheating=prop_cheating,prop_innovative=prop_innovative,cost=cost,connectance,
nbsp_p_dep=nbsp_p_dep,nbsp_a_dep=nbsp_a_dep,nbsp_a=nbsp_a_per,nbsp_p=nbsp_p_per,
pers_tot=(nbsp_a_per+nbsp_p_per)/(nbsp_a_dep+nbsp_p_dep),generalism_cheaters=if(length(cheaters_ini[cheaters_ini==1])>0){mean(apply(Iini[,cheaters_ini==1],2,mean))}else{NA},generalism=mean(apply(Iini,2,mean)),
variance=a,valprop=NA,ai_direct_aa=NA,ai_direct_pp=NA,ai_direct_pa=NA,
ai_direct_ap=NA,ai_net_ap=NA,ai_net_aa=NA,ai_net_pp=NA,ai_net_pa=NA,
ai_indirect_aa=NA,ai_indirect_pp=NA,ai_indirect_pa=NA,ai_indirect_ap=NA,
ai_contrib_ap=NA,ai_contrib_aa=NA,ai_contrib_pp=NA,ai_contrib_pa=NA,shape_a=exp(shape_a),shape_p=exp(shape_p),nb_cheaters_dep=length(cheaters_ini[cheaters_ini==1]),nb_cheaters=0,
cheaters_to_plant=NA,cheaters_to_poll=NA,no_cheaters_to_plant=NA,no_cheaters_to_poll=NA,NODF_Mini=NA,NODF_M=NA,NODF_Tini=NA,
NODF_T=NA,NODF_Iini=NA,NODF_IM=NA))
}

}
print(time1-Sys.time())


}


fwrite(netcar,paste0("/home/duchenne/cheating/eq/netcar_",jj,".txt"),row.names=FALSE,sep="\t")

