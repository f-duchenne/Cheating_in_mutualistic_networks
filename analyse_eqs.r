library(data.table)
library(dplyr)
setwd(dir="C:/Users/Duchenne/Documents/cheating/eq")

lili=list.files()
lili=lili[grep("netcar",lili)]
res=NULL
for(i in lili){
#bidon=fread(paste0("netcar_",i,".txt"))
bidon=fread(i)
bidon2=merge(as.data.frame(bidon[1:nrow(subset(bidon,prop_cheaters==0)),which(!(names(bidon) %in% c("prop_cheaters","prop_innovative","scenario"))),with=F]),
expand.grid(prop_cheaters=unique(bidon$prop_cheaters),prop_cheating=0,prop_innovative=unique(bidon$prop_innovative),cost=unique(bidon$cost),scenario=unique(bidon$scenario),connectance=unique(bidon$connectance)),
by=c("prop_cheating","connectance","cost"),all=T)
if (length(which(is.na(bidon2$pers_tot)))>0) {
break
}
bidon2=bidon2[,names(bidon)]
res=rbind(res,rbind(bidon2,bidon[-(1:nrow(subset(bidon,prop_cheaters==0))),]))
}

fwrite(res,"data_for_analyse.txt")

###########################################
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
setwd(dir="C:/Users/Duchenne/Documents/cheating/eq")

res=fread("data_for_analyse.txt")

res$asy=res$nbsp_p_dep/res$nbsp_a_dep
res$feas=0
res$feas[res$pers_tot==1]=1
res=res %>% dplyr::group_by(essai,cost,connectance) %>% dplyr::mutate(pers_ini=mean(pers_tot[prop_cheaters==0]),resilience_ini=mean(-1*valprop[prop_cheaters==0],na.rm=T))
res$resilience=-1*res$valprop

############ FIGURE 2 ###########
b=subset(res,prop_cheaters>0 & connectance==0.4) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(pers=mean(pers_tot),resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T)) 

pl1=ggplot(data=subset(b,prop_cheaters %in% c(0.1,0.3,0.5) & cost<0.6),aes(x=prop_cheating,y=pers,col=prop_innovative,group=paste0(prop_innovative,scenario),linetype=scenario))+geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)), scales="free")+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Network persistence")+scale_y_continuous(breaks=seq(0,1,0.2),labels = scales::percent_format(accuracy=1))+ggtitle("a")+
scale_color_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])+scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))

b=subset(res,prop_cheaters>0 & connectance==0.4) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(pers=mean(pers_tot),resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T)) 

pl2=ggplot(data=subset(b,prop_cheaters %in% c(0.1,0.3,0.5) & prop_innovative %in% c(0,1) & cost<0.6 & prop_cheating<0.5),aes(x=prop_cheating,y=pers,col=prop_innovative,group=paste0(prop_innovative,scenario),linetype=scenario))+
geom_hline(data=subset(b,prop_cheaters %in% c(0.1,0.3,0.5) & prop_innovative==0 & cost<0.6 & prop_cheating==0),aes(yintercept=pers),linetype="dashed",color="lightgrey")+
geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none")+
labs(color=expression(Psi))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)), scales="free")+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Network persistence")+scale_color_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])+ scale_y_continuous(labels = scales::percent_format(accuracy=1),n.breaks=4)+ggtitle("b")+
scale_x_continuous(breaks=seq(0,0.4,0.2),labels=c(0,0.2,0.4))


leg <- ggpubr::as_ggplot(cowplot::get_legend(pl1))
pl1=pl1+theme(legend.position="none")

grid.arrange(pl1,pl2,leg,layout_matrix =rbind(c(1,3),c(2,3)),widths=c(4,1))

setwd(dir="C:/Users/Duchenne/Documents/cheating")
pdf("fig2.pdf",width=6,height=7)
grid.arrange(pl1,pl2,leg,layout_matrix =rbind(c(1,3),c(2,3)),widths=c(4,1))
dev.off();

############ FIGURE S1 ###########
b=subset(res,prop_cheaters>0 & connectance==0.4) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(pers_cheaters=mean(nb_cheaters/nb_cheaters_dep,na.rm=T)) 

pls=ggplot(data=subset(b,prop_cheaters %in% c(0.1,0.3,0.5) & cost<0.6),aes(x=prop_cheating,y=pers_cheaters,col=prop_innovative,group=paste0(prop_innovative,scenario),linetype=scenario))+
geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)))+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Persistence of cheaters")+scale_y_continuous(labels = scales::percent)+ggtitle("")+
scale_color_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])+scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))

png("fig_s1.png",width=900,height=700,res=150)
pls
dev.off();

############ FIGURE 3 ###########
b=subset(res,prop_cheaters>0 & connectance==0.4) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(pers=mean(pers_tot),resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T)) 

pl1=ggplot(data=subset(b,prop_cheaters %in% c(0.1,0.3,0.5) & cost<0.6),aes(x=prop_cheating,y=resilience,col=prop_innovative,group=paste0(prop_innovative,scenario),linetype=scenario))+
geom_hline(data=subset(b,prop_cheaters %in% c(0.1,0.3,0.5) & prop_innovative==0 & cost<0.6 & prop_cheating==0),aes(yintercept=resilience),linetype="dashed",color="lightgrey")+
geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)), scales="free")+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Network resilience")+
scale_color_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])+scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))


setwd(dir="C:/Users/Duchenne/Documents/cheating")
pdf("fig_3.pdf",width=6,height=3.5)
pl1
dev.off();

############ FIGURE S2 ###########
b=subset(res,prop_cheaters>0 & connectance==0.2) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(pers=mean(pers_tot),resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T)) 

pl1=ggplot(data=subset(b,prop_cheaters %in% c(0.1,0.3,0.5) & cost<0.6),aes(x=prop_cheating,y=pers,col=prop_innovative,group=paste0(prop_innovative,scenario),linetype=scenario))+geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)), scales="free")+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Network persistence")+scale_y_continuous(breaks=seq(0,1,0.2),labels = scales::percent_format(accuracy=1))+ggtitle("a")+
scale_color_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])+scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))


b=subset(res,prop_cheaters>0 & connectance==0.2) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(pers=mean(pers_tot),resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T)) 

pl2=ggplot(data=subset(b,prop_cheaters %in% c(0.1,0.3,0.5) & prop_innovative %in% c(0,1) & cost<0.6 & prop_cheating<0.5),aes(x=prop_cheating,y=pers,col=prop_innovative,group=paste0(prop_innovative,scenario),linetype=scenario))+
geom_hline(data=subset(b,prop_cheaters %in% c(0.1,0.3,0.5) & prop_innovative==0 & cost<0.6 & prop_cheating==0),aes(yintercept=pers),linetype="dashed",color="lightgrey")+
geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none")+
labs(color=expression(Psi))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)), scales="free")+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Network persistence")+scale_color_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])+ scale_y_continuous(labels = scales::percent_format(accuracy=1),n.breaks=4)+ggtitle("b")+
scale_x_continuous(breaks=seq(0,0.4,0.2),labels=c(0,0.2,0.4))

leg <- ggpubr::as_ggplot(cowplot::get_legend(pl1))
pl1=pl1+theme(legend.position="none")

grid.arrange(pl1,pl2,leg,layout_matrix =rbind(c(1,3),c(2,3)),widths=c(4,1))

setwd(dir="C:/Users/Duchenne/Documents/cheating")
png("fig_S2.png",width=900,height=1400,res=150)
grid.arrange(pl1,pl2,leg,layout_matrix =rbind(c(1,3),c(2,3)),widths=c(4,1))
dev.off();

############ FIGURE S3 ###########
b=subset(res,prop_cheaters>0 & connectance==0.2) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(pers=mean(pers_tot),resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T)) 

pl1=ggplot(data=subset(b,prop_cheaters %in% c(0.1,0.3,0.5) & cost<0.6),aes(x=prop_cheating,y=pers,col=prop_innovative,group=paste0(prop_innovative,scenario),linetype=scenario))+geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)))+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Network resilience")+ggtitle("")+
scale_color_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])+scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))

setwd(dir="C:/Users/Duchenne/Documents/cheating")
png("fig_S3.png",width=900,height=700,res=150)
pl1
dev.off();

############ FIGURE S4 ###########

res$NODF_c="intermediate"
res$NODF_c[res$NODF_Iini>quantile(res$NODF_Iini[res$connectance==0.4],prob=0.75)]="high"
res$NODF_c[res$NODF_Iini<quantile(res$NODF_Iini[res$connectance==0.4],prob=0.25)]="low"
b=subset(res,prop_cheaters %in% c(0.1,0.3,0.5) & connectance==0.4 & prop_innovative==1) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,NODF_c,scenario) %>%
summarise(pers=mean(pers_tot),resilience=mean(-1*valprop,na.rm=T),pers_eff=mean(pers_tot)-mean(pers_ini))

pl1=ggplot(data=b,aes(x=prop_cheating,y=pers_eff,col=NODF_c,group=paste0(NODF_c,scenario),linetype=scenario))+
geom_hline(yintercept=0,linetype="dashed",color="lightgrey")+
geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)),scales="free")+
labs(color="NODF")+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Network presistence, relative to a case without cheating")+ggtitle("")+
scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))+scale_color_manual(values=c("firebrick3","deeppink2","pink"))

setwd(dir="C:/Users/Duchenne/Documents/cheating")
png("fig_S4.png",width=1000,height=800,res=150)
pl1
dev.off();

############ FIGURE S5 ###########

b=subset(res,prop_cheaters %in% c(0.1,0.3,0.5) & connectance==0.4 & nbsp_a>4 & nbsp_p>4) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(pers=mean(pers_tot),resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),
contrib=mean(cheaters_to_poll,na.rm=T),n=length(cheaters_to_poll[!is.na(cheaters_to_poll)])) 

png("fig_S4.png",width=1000,height=800,res=150)
pl1=ggplot(data=subset(b,cost==0),aes(x=prop_cheating,y=contrib,col=prop_innovative,group=paste0(prop_innovative,scenario),linetype=scenario))+geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)),scales="free")+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Average total effect of cheaters on other pollinators")+ggtitle("")+
scale_color_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])+scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))

sum(b$n)

setwd(dir="C:/Users/Duchenne/Documents/cheating")
png("fig_S5.png",width=1000,height=600,res=150)
pl1
dev.off();








res$pers_eff=res$pers_tot-res$pers_ini
b2=subset(res,prop_cheaters %in% c(0.1,0.3,0.5) & connectance==0.4 & prop_innovative==1 & scenario=="specialists") %>% group_by(prop_innovative,cost,prop_cheaters,NODF_Iini,scenario,essai) %>%
summarise(seuil=max(prop_cheating[pers_eff>0]))

ggplot(data=b2,aes(y=NODF_Iini,x=as.factor(seuil),shape=scenario))+
geom_boxplot()+theme_bw()+coord_flip()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)),scales="free")+
labs(color=expression(phi1))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Network resilience")+ggtitle("a")+
scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))+scale_color_manual(values=c("dodgerblue3","gold3"))




res$mod_c="low"
res$mod_c[res$mod_Iini>mean(res$mod_Iini[res$connectance==0.4])]="high"
b=subset(res,prop_cheaters %in% c(0.1,0.3,0.5) & connectance==0.4 & scenario=="specialists" & prop_innovative==1) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,mod_c) %>% summarise(pers=mean(pers_tot),resilience=mean(-1*valprop,na.rm=T)) 

ggplot(data=b,aes(x=prop_cheating,y=pers,col=mod_c,group=mod_c))+geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)))+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Network resilience")+ggtitle("a")+
scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))


b=subset(res,prop_cheaters>0) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(pers_c0.2=mean(pers_tot[connectance==0.2]),pers_c0.4=mean(pers_tot[connectance==0.4]))
ggplot(data=subset(b,prop_cheaters %in% c(0.1,0.3,0.5) & cost<0.6),aes(x=pers_c0.4,y=pers_c0.2,col=prop_innovative,group=paste0(prop_innovative,scenario),shape=scenario))+geom_point()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)),scales="free")+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Network persistence")+scale_y_continuous(labels = scales::percent)+ggtitle("a")+
scale_color_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])+scale_x_continuous(labels = scales::percent)






b=subset(res,prop_cheaters>0) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(resilience_c0.2=mean(resilience[connectance==0.2],na.rm=T),resilience_c0.4=mean(resilience[connectance==0.4],na.rm=T))
ggplot(data=subset(b,prop_cheaters %in% c(0.1,0.3,0.5) & cost<0.6 & prop_cheating==0.2),aes(x=resilience_c0.4,y=resilience_c0.2,col=prop_innovative,group=paste0(prop_innovative,scenario),shape=scenario))+geom_point()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)),scales="free")+
labs(color=expression(Psi))+
xlab("Network resilience")+ylab("Network resilience")+scale_y_continuous()+ggtitle("a")+
scale_color_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])+scale_x_continuous()
















res$eff_pers=0
res$eff_pers[res$pers_ini<res$pers_tot]=1
res$eff_res=0
res$eff_res[res$resilience_ini<(-1*res$valprop)]=1
res$eff="nothing"
res$eff[res$eff_pers==1 & res$eff_res==0]="persistence only"
res$eff[res$eff_pers==0 & res$eff_res==1]="resilience only"
res$eff[res$eff_pers==1 & res$eff_res==1]="persistence & resilience"

b=subset(res,prop_cheaters %in% c(0.1,0.3,0.5) & connectance==0.4) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(moy=mean(eff_pers+eff_res,na.rm=T)) 

pl1=ggplot(data=subset(b,cost==0 & scenario=="specialists"),aes(x=prop_cheating,y=prop_innovative,fill=moy))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_wrap(~prop_cheaters,ncol=1)+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab(expression(paste("Proportion of innovative cheating (",Psi,")")))+ggtitle("")+
scale_fill_viridis()+scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))+labs(fill="Positive effect on:")+
scale_y_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))



b=subset(res,prop_cheaters==0.5 & connectance==0.4 & cost==0) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(pollinators=mean(cheaters_to_poll-no_cheaters_to_poll,na.rm=T),
plants=mean(cheaters_to_plant/no_cheaters_to_plant,na.rm=T)) 

ggplot(data=b,aes(x=prop_cheating,y=pollinators,linetype=scenario,col=prop_innovative,group=paste0(prop_innovative,scenario)))+
geom_line()+
theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_wrap(~scenario)+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Total effect of cheaters to")+ggtitle("a")+
scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))+
scale_color_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])






















Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
b=subset(res,prop_cheaters>0 & connectance==0.4) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(moy=Mode(eff)) 
b$moy=factor(b$moy,levels=c("nothing","resilience only","persistence only","persistence & resilience"))

pl1=ggplot(data=subset(b,prop_cheaters %in% c(0.1,0.3,0.5) & cost==0 & scenario=="specialists"),aes(x=prop_cheating,y=prop_innovative,fill=moy))+geom_tile()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(cols=vars(prop_cheaters))+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Proportion of innovative cheating")+ggtitle("")+
scale_fill_viridis(discrete=T)+scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))+labs(fill="Positive effect on:")


b=subset(res,prop_cheaters %in% c(0.1,0.3,0.5) & cost==0 & scenario=="specialists") %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(persistence_only=length(eff[eff=="persistence only"]),
nothing=length(eff[eff=="nothing"]),resilience_only=length(eff[eff=="resilience only"]),persistence_resilience=length(eff[eff=="persistence & resilience"])) 

b[,6:9]=apply(b[,6:9],2,as.numeric)

ggplot(data=b,aes(x=prop_cheating,y=prop_innovative))+geom_scatterpie(data=b,aes(x=prop_cheating,y=prop_innovative),cols=c("persistence_only","nothing","resilience_only","persistence_resilience"),color=NA,radius=50)+theme_bw()


theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(cols=vars(prop_cheaters))+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Proportion of innovative cheating")+ggtitle("")+
scale_fill_viridis(discrete=T)+scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))+labs(fill="Positive effect on:")





ggplot(data=b,aes(x=prop_cheating,y=resilience,col=prop_innovative,group=prop_innovative))+geom_line()+theme_bw()+theme(panel.grid=element_blank())+facet_grid(rows=vars(cost),cols=vars(prop_cheaters))

ggplot(data=b,aes(x=pers,y=resilience,col=prop_innovative))+geom_point()+theme_bw()+theme(panel.grid=element_blank())+facet_grid(rows=vars(cost),cols=vars(prop_cheaters))



ggplot(data=b,aes(x=prop_cheating,y=contrib,col=prop_innovative,group=prop_innovative))+geom_line()+theme_bw()+theme(panel.grid=element_blank())+facet_grid(rows=vars(cost),cols=vars(prop_cheaters))+
xlab("Cheating frequency of cheaters (Ω)")+ylab("Persistence")+scale_color_viridis()

































setwd(dir="C:/Users/Duchenne/Documents/cheating/initial")
NO=NULL
for(jj in 1:500){
connectance=0.2
load(paste0("replicate_",jj,".RData"))
Iini=IPOLL
Iini[Iini<quantile(c(Iini),probs=(1-connectance))]=0
Iini[Iini>0]=1
check_a=apply(Iini,2,sum)
check_p=apply(Iini,1,sum)
if(min(check_a)==0){
for(i in which(check_a==0)){
Iini[which.max(IPOLL[,i]),i]=1
}
}
if(min(check_p)==0){
for(i in which(check_p==0)){
Iini[i,which.max(IPOLL[i,])]=1
}
}

bidon1=data.frame(essai=jj,NODF=networklevel(Iini,index="NODF")[[1]],connectance=connectance)
connectance=0.4
load(paste0("replicate_",jj,".RData"))
Iini=IPOLL
Iini[Iini<quantile(c(Iini),probs=(1-connectance))]=0
Iini[Iini>0]=1
check_a=apply(Iini,2,sum)
check_p=apply(Iini,1,sum)
if(min(check_a)==0){
for(i in which(check_a==0)){
Iini[which.max(IPOLL[,i]),i]=1
}
}
if(min(check_p)==0){
for(i in which(check_p==0)){
Iini[i,which.max(IPOLL[i,])]=1
}
}

bidon2=data.frame(essai=jj,NODF=networklevel(Iini,index="NODF")[[1]],connectance=connectance)

NO=rbind(NO,bidon1,bidon2)
}
res=merge(res,NO,by=c("essai","connectance"))
