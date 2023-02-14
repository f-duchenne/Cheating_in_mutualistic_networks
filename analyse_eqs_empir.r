library(data.table)
library(dplyr)
setwd(dir="C:/Users/Duchenne/Documents/cheating/eq_empir/")

lili=list.files()
res=NULL
for(i in lili){
#bidon=fread(paste0("netcar_",i,".txt"))
bidon=fread(i)
res=rbind(res,bidon)
}

setwd(dir="C:/Users/Duchenne/Documents/cheating/")
fwrite(unique(res),"data_for_analyse_empir.txt")

library(bipartite)
library(Rmpfr)
library(circular)
library(CircStats)
library(plot3D)
library(ggradar)
library(ggplot2)
library(gridExtra)
library(data.table)
library(dplyr)
library(rgbif)
library(ggforce)
library(ggeffects)
library(ggExtra)
library(viridis)
library(lme4)
library(cowplot)
library(scales)
library(car)
library(DHARMa)
library(glmmTMB)
library(qgraph)
library(igraph)
library(piecewiseSEM)
library(randomForest)
library(FactoMineR)
require(factoextra)
library("ggdendro")
library(dendextend)
library(ape)
library(ggtree)
library(ggnewscale)
library(geiger)
library(diversityForest)
library(caper)
library(phytools)
library(ggthemes)
library(grid)
library(tidyverse)
source("C:/Users/Duchenne/Documents/cheating/script/theme_border.r")
setwd(dir="C:/Users/Duchenne/Documents/cheating/")

#LOAD SITE METADATARAIT
site_cr=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_2022-11-07/Site_metadata_Costa-Rica.txt")
site_ec=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/Ecuador_2022-11-07/Site_metadata_Ecuador.txt")
sites=rbind(site_cr,site_ec)
resi=fread("resi.txt")
sites=merge(sites,resi,by="site")

res=fread("data_for_analyse_empir.txt")
res=merge(res,sites,by="site")
res$resilience=-1*res$valprop

res$simulation[res$simulation=="tout"]="Complete network"
res$simulation[res$simulation=="without_cheat"]="Without cheating"
res$pers_p=res$nbsp_p/res$nbsp_p_dep
res$pers_a=res$nbsp_a/res$nbsp_a_dep


#### SELECT SIMULATIONS:
res=subset(res, efficience<=2 & cheater_cost=="Escape")

b2=merge(subset(res,simulation=="Complete network"),subset(res,simulation=="Without cheating",select=c("site","cheater_cost","interf","cost","pers_tot","efficience","essai")),
by=c("site","cheater_cost","cost","efficience","essai","interf"))
b2$pers=b2$pers_tot.x-b2$pers_tot.y

b=b2 %>% group_by(interf,cost,efficience,Country,site) %>% summarise(pers_moy=mean(pers),sde=sd(pers)/sqrt(length(pers)))
b$Country=factor(b$Country,levels=c("Ecuador","Costa-Rica"))
b$site=factor(b$site,levels=c("AMIG","BOQU","CUSI","GIGA","LONG","MILL","NIMB","NUBE","QUEB","RIOM","SANM","TOLO","Alaspungo","Alaspungo_disturbed","Verdecocha","Yanacocha","Yanacocha_disturbed"))

b3=b2 %>% group_by(interf,cost,efficience) %>% summarise(pers_moy=mean(pers),sde=sd(pers)/sqrt(length(pers)))

dodge=0.5

pl1=ggplot(data=b,aes(x=as.factor(cost),y=pers_moy))+
geom_hline(yintercept=0,color="red")+
geom_errorbar(aes(ymin=pers_moy-1.96*sde,ymax=pers_moy+1.96*sde,col=as.factor(Country),group=site),width=0,position=position_dodge(width=dodge),alpha=0.5)+
geom_point(aes(fill=as.factor(Country),group=site),position=position_dodge(width=dodge),shape=21,col="white",alpha=0.5)+
geom_pointrange(data=b3,aes(x=as.factor(cost),y=pers_moy,ymin=pers_moy-1.96*sde,ymax=pers_moy+1.96*sde),col="black",shape=4)+
theme_classic()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")))+
facet_grid(cols=vars(efficience),rows=vars(interf),labeller = label_bquote(cols= italic(alpha)== .(efficience),rows= italic(c)== .(interf)),scales="free")+
scale_color_manual(values=c("#0B4F6C","#CBB9A8"))+scale_fill_manual(values=c("#0B4F6C","#CBB9A8"))+
ylab("Effect of cheating on network persistence")+xlab(expression(paste("Cost associated with mutualism (",Lambda,")")))+
theme(strip.placement = "outside")+labs(col="Country",fill="Country")

setwd(dir="C:/Users/Duchenne/Documents/cheating")
pdf("figure_5.pdf",width=10,height=6)
pl1
dev.off();














ggplot(data=subset(b2,efficience %in% c(1.4,2)))+
#geom_jitter(data=res,aes(x=as.factor(site),y=pers_tot,col=simulation),alpha=0.5,size=0.5,position=position_jitterdodge(jitter.width=0.1,dodge.width=0.2))+
geom_point(aes(x=resi,y=pers,col=Country))+geom_hline(yintercept=0,color="red")+
theme_classic()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")))+
facet_grid(cols=vars(cost),rows=vars(efficience),labeller = label_bquote(cols=Lambda == .(cost)),scales="free")+
scale_color_manual(values=c("#0B4F6C","#CBB9A8"))+ylab("Effect of cheating on network persistence")+scale_fill_manual(values=c("#0B4F6C","#CBB9A8"))+xlab("Overall level of cheating (residuals)")+
theme(strip.placement = "outside")+stat_smooth(method="lm",aes(x=resi,y=pers,col=Country,fill=Country),alpha=0.1)


ggplot()+
#geom_jitter(data=res,aes(x=as.factor(site),y=pers_tot,col=simulation),alpha=0.5,size=0.5,position=position_jitterdodge(jitter.width=0.1,dodge.width=0.2))+
geom_point(data=b2,aes(x=min_transect_elev,y=pers,col=Country))+geom_hline(yintercept=0,color="red")+
theme_classic()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")))+
facet_grid(cols=vars(cost),rows=vars(cheater_cost),labeller = label_bquote(cols=Lambda == .(cost)),scales="free")+
scale_color_manual(values=c("#0B4F6C","#CBB9A8"))+ylab("Effect of cheating on network persistence")+scale_fill_manual(values=c("#0B4F6C","#CBB9A8"))+xlab("Overall level of cheating (residuals)")+
theme(strip.placement = "outside")+stat_smooth(method="lm",data=b2,aes(x=min_transect_elev,y=pers,col=Country,fill=Country),alpha=0.1)



ggplot()+
#geom_jitter(data=res,aes(x=as.factor(site),y=pers_tot,col=simulation),alpha=0.5,size=0.5,position=position_jitterdodge(jitter.width=0.1,dodge.width=0.2))+
geom_point(data=b2,aes(x=resi,y=feas,col=Country))+geom_hline(yintercept=0,color="red")+
theme_classic()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")))+
facet_grid(cols=vars(cost),rows=vars(cheater_cost),labeller = label_bquote(cols=Lambda == .(cost)),scales="free")+
scale_color_manual(values=c("#0B4F6C","#CBB9A8"))+ylab("Effect of cheating on network persistence")+scale_fill_manual(values=c("#0B4F6C","#CBB9A8"))+xlab("Overall level of cheating (residuals)")+
theme(strip.placement = "outside")+stat_smooth(method="lm",data=b2,aes(x=resi,y=feas,col=Country,fill=Country),alpha=0.1)




ggplot()+
geom_point(data=subset(b,self_cros=="yes"),aes(x=as.factor(site),y=feas,fill=simulation),col="white",shape=21,size=2,position=position_dodge(width=0.2))+
theme_classic()+coord_flip()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")))+
facet_grid(cols=vars(cost),rows=vars(interf),labeller = label_bquote(cols=Lambda == .(cost)), space ="free",scales="free")+
scale_color_manual(values=c("#16E0BD","#98838F"))+ylab("Network persistence")+scale_fill_manual(values=c("#16E0BD","#98838F"))+ylab("Network persistence")+
xlab("Sites")+
theme(strip.placement = "outside")



b=res %>% group_by(site,simulation,cost,cheater_cost,efficience,min_transect_elev,interf) %>% summarise(pers_moy=mean(pers_tot),sde=sd(pers_tot)/sqrt(length(pers_tot)),
feas=length(pers_tot[pers_tot==1])/length(pers_tot),res=mean(resilience,na.rm=T))

ggplot(data=subset(b,cheater_cost=="Reduce" & simulation=="Complete network"),aes(x=efficience,y=cost,fill=feas))+
geom_tile()+
theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_wrap(~site+interf,ncol=4)+
scale_color_manual(values=c("#0B4F6C","#CBB9A8"))+ylab("Effect of cheating on network persistence")+scale_fill_viridis()+
theme(strip.placement = "outside")+coord_fixed(expand=F)



b=res %>% group_by(Country,site,simulation,cost,cheater_cost,min_transect_elev,interf,efficience,resi) %>% summarise(pers_moy=mean(pers_tot),sde=sd(pers_tot)/sqrt(length(pers_tot)),
feas=length(pers_tot[pers_tot==1])/length(pers_tot),res=mean(resilience,na.rm=T),pers_p=mean(pers_p),pers_a=mean(pers_a))
b2=merge(subset(b,simulation=="Complete network"),subset(b,simulation=="Without cheating",select=c("site","cheater_cost","efficience","cost","interf","pers_moy","feas","res","pers_p","pers_a")),by=c("site","cheater_cost","efficience","cost","interf"))
b2$pers=b2$pers_moy.x-b2$pers_moy.y
b2$pers_p=b2$pers_p.x-b2$pers_p.y
b2$pers_a=b2$pers_a.x-b2$pers_a.y
b2$feas=b2$feas.x-b2$feas.y
b2$res=b2$res.x-b2$res.y
b2$Country=factor(b2$Country,levels=c("Ecuador","Costa-Rica"))
ggplot()+
#geom_jitter(data=res,aes(x=as.factor(site),y=pers_tot,col=simulation),alpha=0.5,size=0.5,position=position_jitterdodge(jitter.width=0.1,dodge.width=0.2))+
geom_boxplot(data=b2,aes(x=as.factor(efficience),y=pers_p,col=as.factor(Country)),position=position_dodge2())+geom_hline(yintercept=0,color="red")+
theme_classic()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")))+
facet_grid(cols=vars(cost),rows=vars(cheater_cost,interf),labeller = label_bquote(cols=Lambda == .(cost)),scales="free")+
scale_color_manual(values=c("#0B4F6C","#CBB9A8"))+ylab("Effect of cheating on network persistence")+scale_fill_manual(values=c("#16E0BD","#98838F"))+xlab("Mutualistic strength")+
theme(strip.placement = "outside")


ggplot()+
#geom_jitter(data=res,aes(x=as.factor(site),y=pers_tot,col=simulation),alpha=0.5,size=0.5,position=position_jitterdodge(jitter.width=0.1,dodge.width=0.2))+
geom_point(data=b2,aes(x=resi,y=feas,col=Country),position=position_dodge2())+geom_hline(yintercept=0,color="red")+
theme_classic()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")))+
facet_grid(cols=vars(cost),rows=vars(cheater_cost,interf),labeller = label_bquote(cols=Lambda == .(cost)),scales="free")+
scale_color_manual(values=c("#0B4F6C","#CBB9A8"))+ylab("Effect of cheating on network persistence")+scale_fill_manual(values=c("#16E0BD","#98838F"))+xlab("Mutualistic strength")+
theme(strip.placement = "outside")



ggplot()+
#geom_jitter(data=res,aes(x=as.factor(site),y=pers_tot,col=simulation),alpha=0.5,size=0.5,position=position_jitterdodge(jitter.width=0.1,dodge.width=0.2))+
geom_point(data=subset(b2,self_cros=="no"),aes(x=min_transect_elev,y=pers,color=as.factor(Country)),position=position_dodge2())+geom_hline(yintercept=0,color="red")+
theme_classic()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")))+
facet_grid(cols=vars(cost),rows=vars(efficience),labeller = label_bquote(cols=Lambda == .(cost)),scales="free")+
ylab("Network persistence")+ylab("Network persistence")+scale_color_manual(values=c("#0B4F6C","#CBB9A8"))+
xlab("Sites")+
theme(strip.placement = "outside")

ggplot()+
#geom_jitter(data=res,aes(x=as.factor(site),y=pers_tot,col=simulation),alpha=0.5,size=0.5,position=position_jitterdodge(jitter.width=0.1,dodge.width=0.2))+
geom_point(data=subset(b2,self_cros=="no"),aes(x=min_transect_elev,y=pers,color=as.factor(Country)),position=position_dodge2())+geom_hline(yintercept=0,color="red")+
theme_classic()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")))+
facet_grid(cols=vars(cost),rows=vars(efficience),labeller = label_bquote(cols=Lambda == .(cost)),scales="free")+
ylab("Network persistence")+ylab("Network persistence")+scale_color_manual(values=c("#0B4F6C","#CBB9A8"))+
xlab("Sites")+
theme(strip.placement = "outside")



b=subset(res,self_cros==1) %>% group_by(Country,site,simulation,cost) %>% summarise(pers_moy=mean(pers_tot),sde=sd(pers_tot)/sqrt(length(pers_tot)),feas=length(pers_tot[pers_tot==1])/length(pers_tot))

model=glmmTMB(pers_tot~simulation+Country+cost+self_cros+efficience+(1|site/essai),family=gaussian,data=res)

ggplot()+
#geom_jitter(data=res,aes(x=as.factor(site),y=pers_tot,col=simulation),alpha=0.5,size=0.5,position=position_jitterdodge(jitter.width=0.1,dodge.width=0.2))+
geom_errorbar(data=b,aes(x=as.factor(site),y=pers_moy,ymin=pers_moy-sde,ymax=pers_moy+sde,color=simulation),width=0,position=position_dodge(width=0.2))+
geom_point(data=b,aes(x=as.factor(site),y=pers_moy,fill=simulation),col="white",shape=21,size=2,position=position_dodge(width=0.2))+
theme_classic()+coord_flip()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")))+
facet_grid(cols=vars(cost),rows=vars(Country),labeller = label_bquote(cols=Lambda == .(cost)), space ="free",scales="free")+
scale_color_manual(values=c("#16E0BD","#98838F"))+ylab("Network persistence")+scale_fill_manual(values=c("#16E0BD","#98838F"))+ylab("Network persistence")+
xlab("Sites")+
theme(strip.placement = "outside")


ggplot()+
#geom_jitter(data=res,aes(x=as.factor(site),y=pers_tot,col=simulation),alpha=0.5,size=0.5,position=position_jitterdodge(jitter.width=0.1,dodge.width=0.2))+
geom_point(data=b,aes(x=as.factor(site),y=feas,fill=simulation),col="white",shape=21,size=2,position=position_dodge(width=0.2))+
theme_classic()+coord_flip()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")))+
facet_grid(cols=vars(cost),rows=vars(Country),labeller = label_bquote(cols=Lambda == .(cost)), space ="free",scales="free")+
scale_color_manual(values=c("#16E0BD","#98838F"))+ylab("Network persistence")+scale_fill_manual(values=c("#16E0BD","#98838F"))+ylab("Network persistence")+
xlab("Sites")+
theme(strip.placement = "outside")

b=res %>% group_by(Country,site,simulation,cost) %>% summarise(pers_moy=mean(pers_tot),resilience=mean(resilience,na.rm=T))

ggplot(data=b,aes(x=as.factor(Country),y=pers_moy,col=simulation))+geom_boxplot(position=position_dodge(width=0.9))+
theme_classic()+geom_jitter(alpha=0.5,width=0.1)+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")))+
facet_grid(cols=vars(cost),labeller = label_bquote(cols=Lambda == .(cost)), space ="free",scales="free")+
scale_color_manual(values=c("#16E0BD","#98838F"))+ylab("Network persistence")+xlab("")+
theme(strip.placement = "outside")



b=subset(res,pers_tot>=0.5) %>% group_by(Country,site,simulation,cost) %>% summarise(res_moy=mean(resilience),sde=sd(resilience)/sqrt(length(resilience)),feas=length(pers_tot[pers_tot==1]))

ggplot()+
#geom_jitter(data=res,aes(x=as.factor(site),y=resilience,col=simulation),alpha=0.5,size=0.5,position=position_jitterdodge(jitter.width=0.1,dodge.width=0.2))+
geom_errorbar(data=b,aes(x=as.factor(site),y=res_moy,ymin=res_moy-sde,ymax=res_moy+sde,color=simulation),width=0,position=position_dodge(width=0.2))+
geom_point(data=b,aes(x=as.factor(site),y=res_moy,fill=simulation),col="white",shape=21,size=2,position=position_dodge(width=0.2))+
theme_classic()+coord_flip()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")))+
facet_grid(cols=vars(cost),rows=vars(Country),labeller = label_bquote(cols=Lambda == .(cost)), space ="free",scales="free")+
scale_color_manual(values=c("#16E0BD","#98838F"))+ylab("Network persistence")+scale_fill_manual(values=c("#16E0BD","#98838F"))+ylab("Network persistence")+
xlab("Sites")+
theme(strip.placement = "outside")






