library(ggplot2)
library(gridExtra)
library(data.table)
library(dplyr)
library(ggforce)
library(ggeffects)
library(ggExtra)
library(viridis)
library(cowplot)
library(scales)
library(car)
source("C:/Users/Duchenne/Documents/cheating/script/theme_border.r")
setwd(dir="C:/Users/Duchenne/Documents/cheating/")

#LOAD SITE METADATARAIT
sites=fread("table_s2.csv")
resi=fread("resi.txt")

#LOAD RESULTS SIMULATIONS
res=fread("equilibriums_analyse_empir.txt")
res=merge(res,sites,by="site")

res$simulation[res$simulation=="tout"]="Complete network"
res$simulation[res$simulation=="without_cheat"]="Without cheating"
res$pers_p=res$nbsp_p/res$nbsp_p_dep
res$pers_a=res$nbsp_a/res$nbsp_a_dep

#### EFFECT OF CHEATING ON PERSISTENCE:
b2=merge(subset(res,simulation=="Complete network"),subset(res,simulation=="Without cheating",select=c("site","interf","cost","pers_tot","efficience","essai")),
by=c("site","cost","efficience","essai","interf"))
b2$pers=b2$pers_tot.x-b2$pers_tot.y

b=b2 %>% group_by(interf,cost,efficience,Country,site) %>% summarise(pers_moy=mean(pers),sde=sd(pers)/sqrt(length(pers)))
b$Country=gsub("-"," ",b$Country,fixed=T)
b$Country=factor(b$Country,levels=c("Ecuador","Costa Rica"))
b$site=factor(b$site,levels=c("AMIG","BOQU","CUSI","GIGA","LONG","MILL","NIMB","NUBE","QUEB","RIOM","SANM","TOLO","Alaspungo","Alaspungo_disturbed","Verdecocha","Yanacocha","Yanacocha_disturbed"))

b3=b2 %>% group_by(interf,cost,efficience) %>% summarise(pers_moy=mean(pers),sde=sd(pers)/sqrt(length(pers)))

dodge=0.5

#PLOT
pl1=ggplot(data=b,aes(x=as.factor(cost),y=pers_moy))+
geom_hline(yintercept=0,color="red")+
geom_errorbar(aes(ymin=pers_moy-1.96*sde,ymax=pers_moy+1.96*sde,col=as.factor(Country),group=site),width=0,position=position_dodge(width=dodge),alpha=0.5)+
geom_point(aes(fill=as.factor(Country),group=site),position=position_dodge(width=dodge),shape=21,col="white",alpha=0.5)+
geom_pointrange(data=b3,aes(x=as.factor(cost),y=pers_moy,ymin=pers_moy-1.96*sde,ymax=pers_moy+1.96*sde),col="black",shape=4)+
theme_classic()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")),axis.text.x=element_text(angle=45,hjust=1))+
facet_grid(cols=vars(efficience),rows=vars(interf),labeller = label_bquote(cols= italic(alpha)== .(efficience),rows= italic(c)== .(interf)),scales="free")+
scale_color_manual(values=c("#0B4F6C","#CBB9A8"))+scale_fill_manual(values=c("#0B4F6C","#CBB9A8"))+
ylab("Effect of cheating on network persistence")+xlab(expression(paste("Cost associated with mutualism (",Lambda,")")))+
theme(strip.placement = "outside")+labs(col="Country",fill="Country")

setwd(dir="C:/Users/Duchenne/Documents/cheating")
pdf("figure_5.pdf",width=10,height=6)
pl1
dev.off();

##### POSTER:

ggplot(data=subset(b,efficience<1.5 & interf==1 & cost<0.3),aes(x=as.factor(cost),y=pers_moy,fill=cost))+
geom_hline(yintercept=0,color="red")+
geom_violin(width=1.4,alpha=0.4) +
geom_boxplot(width=0.1, color="grey", alpha=0.8)+
theme_classic()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")),axis.text.x=element_text(angle=0),legend.position="none")+
scale_fill_viridis()+
ylab("Effect of cheating on community persistence")+xlab("Cost associated with mutualism")+
theme(strip.placement = "outside")

#################################################### ANALYSING NULL MODEL RESULTS

res=b2[,c("site","cost","efficience","essai","interf","min_transect_elev","Country","pers","pers_tot.x","pers_tot.y")]

#LOAD RESULTS SIMULATIONS NULL MODEL
resnull=fread("equilibriums_analyse_empir_null.txt")

#MERGE BOTH
resf=merge(res,resnull,by=c("site","cost","efficience","essai","interf"))

#CONVERT TO LOGIT SCALE 
resf$persnull_l=car::logit(resf$pers_tot,adjust=0.01)
resf$pers_l=car::logit(resf$pers_tot.x,adjust=0.01)
resf$pers_wc_l=car::logit(resf$pers_tot.y,adjust=0.01)


#CALCULATE ZSCORE FOR EACH TRIAL
b4=resf %>% dplyr::group_by(site,cost,efficience,essai,interf,Country,min_transect_elev) %>% dplyr::summarise(null_mean=mean(persnull_l),null_sd=sd(persnull_l),true=mean(pers_l),absolute_eff=mean(pers),
without_cheating=mean(pers_wc_l),nbsp=mean(nbsp_a_dep+nbsp_p_dep))
b4=b4 %>% dplyr::group_by(site,cost,efficience,interf,Country,min_transect_elev) %>% dplyr::mutate(min_null_sd=min(null_sd[null_sd!=0],na.rm=T))
b4$zscore=(b4$true-b4$null_mean)/b4$null_sd
b4$zscore_wc=(b4$without_cheating-b4$null_mean)/b4$null_sd
b4$zscore[b4$zscore==-Inf]=-5 ### IF Z-score is -Inf or +Inf, set it to -5 or +5, respectively
b4$zscore[b4$zscore==Inf]=5 ### IF Z-score is -Inf or +Inf, set it to -5 or +5, respectively

b4=merge(b4,resi,by="site")

#AVERAGE ZSCORE PER SITE
b5=subset(b4,!is.infinite(zscore) & null_sd>0) %>% dplyr::group_by(site,cost,efficience,interf,Country,min_transect_elev,cheat) %>% dplyr::summarise(zscore_persite=mean(zscore,na.rm=T),nb=length(zscore),
absolute_eff=mean(absolute_eff),zscore_wc_persite=mean(zscore_wc))


#AVERAGE ZSCORE ACROSS SITES
b6=b5 %>% group_by(interf,cost,efficience) %>% summarise(zscore_moy=mean(zscore_persite,na.rm=T),sde=sd(zscore_persite,na.rm=T)/sqrt(length(zscore_persite[!is.na(zscore_persite)])),
zscore_wc_moy=mean(zscore_wc_persite,na.rm=T),sde_wc=sd(zscore_wc_persite,na.rm=T)/sqrt(length(zscore_wc_persite[!is.na(zscore_wc_persite)])),)

ggplot(data=b5,aes(x=as.factor(cost),y=zscore_persite))+
geom_point(position=position_dodge(width=dodge),aes(fill=Country,color=Country,group=site,size=cheat),alpha=0.6)+geom_hline(yintercept=c(-1.96,1.96),linetype="dashed")+
geom_hline(yintercept=0)+
theme_classic()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")),axis.text.x=element_text(angle=45,hjust=1))+
facet_grid(cols=vars(efficience),rows=vars(interf),labeller = label_bquote(cols= italic(alpha)== .(efficience),rows= italic(c)== .(interf)))+
scale_color_manual(values=c("#CBB9A8","#0B4F6C"))+scale_fill_manual(values=c("#CBB9A8","#0B4F6C"))+
theme(strip.placement = "outside")+labs(col="Country",fill="Country",size="Overall level\nof cheating")+
geom_pointrange(data=b6,aes(x=as.factor(cost),y=zscore_moy,ymin=zscore_moy-1.96*sde,ymax=zscore_moy+1.96*sde),col="black",shape=4)+
ylab("Z-score of persistence due to observed cheating patterns\n(relative to randomized ones)")+xlab(expression(paste("Cost associated with mutualism (",Lambda,")")))


dodge=0.3
setwd(dir="C:/Users/Duchenne/Documents/cheating/")
png("Figure_S8.png",width=1000,height=1200,res=120)
ggplot(data=b5,aes(x=as.factor(cost),y=zscore_persite))+
geom_point(position=position_dodge(width=dodge),aes(fill=Country,color=Country,group=site,size=cheat),alpha=0.6)+geom_hline(yintercept=c(-1.96,1.96),linetype="dashed")+
geom_hline(yintercept=0)+
theme_classic()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),panel.border = theme_border(type = c("bottom","left")),axis.text.x=element_text(angle=45,hjust=1))+
facet_grid(cols=vars(efficience),rows=vars(interf),labeller = label_bquote(cols= italic(alpha)== .(efficience),rows= italic(c)== .(interf)))+
scale_color_manual(values=c("#CBB9A8","#0B4F6C"))+scale_fill_manual(values=c("#CBB9A8","#0B4F6C"))+
theme(strip.placement = "outside")+labs(col="Country",fill="Country",size="Overall level\nof cheating")+
geom_pointrange(data=b6,aes(x=as.factor(cost),y=zscore_moy,ymin=zscore_moy-1.96*sde,ymax=zscore_moy+1.96*sde),col="black",shape=4)+
ylab("Z-score of persistence due to observed cheating patterns\n(relative to randomized ones)")+xlab(expression(paste("Cost associated with mutualism (",Lambda,")")))
dev.off();