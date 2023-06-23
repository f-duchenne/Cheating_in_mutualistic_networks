library(data.table)
library(dplyr)
setwd(dir="C:/Users/Duchenne/Documents/cheating/eq_tot")

lili=list.files()
lili=lili[grep("netcar",lili)]
res=NULL
for(i in lili){
#bidon=fread(paste0("netcar_",i,".txt"))
bidon=fread(i)
bidon=bidon[,c("prop_cheating","prop_innovative","cost","prop_cheaters","scenario","connectance","essai","pers_tot","nbsp_a","nbsp_p","nb_cheaters","nb_cheaters_dep","valprop","C_Iini","NODF_Iini","mod_Iini",
"ai_contrib_aa","cheaters_to_poll","cheaters_to_plant"),with=F]
res=rbind(res,bidon)
}


bidon=data.table(essai=1:500,avg_ra=NA,avg_rp=NA)
for(j in 1:500){
load(paste0("C:/Users/Duchenne/Documents/cheating/initial/replicate_",j,".RData"))
bidon$avg_ra[bidon$essai==j]=mean(r[1:20])
bidon$avg_rp[bidon$essai==j]=mean(r[21:40])
}

dim(res)
res=merge(res,bidon,by="essai")
dim(res)

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
library(gMOIP)
library(rgl)
library(cxhull)
library(magick)
library(partR2)
library(ggpubr)
setwd(dir="C:/Users/Duchenne/Documents/cheating/eq_tot")

res=fread("data_for_analyse.txt")

res$feas=0
res$feas[res$pers_tot==1]=1
res$resilience=-1*res$valprop
nrow(subset(res,resilience<0))
res=res %>% group_by(prop_innovative,cost,prop_cheaters,scenario,connectance,essai) %>% mutate(pers0=pers_tot[prop_cheating==0],resilience0=resilience[prop_cheating==0])

############ FIGURE 2 ###########
b=subset(res,prop_cheaters>0 & connectance==0.4) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(persr=mean(pers_tot-pers0),persr_sde=sd(pers_tot-pers0)/sqrt(length(persr)),
resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T)) 

bidon=subset(b,cost<=0.15 & prop_cheaters %in% c(0.1,0.5))
zero_pos=unique(scales::rescale(bidon$persr, to = c(0, 1))[bidon$persr==0])


pl1=
ggplot(data=subset(b,cost<=0.15 & prop_cheaters==0.1),aes(x=prop_cheating,y=prop_innovative,fill=persr))+
theme_bw()+
geom_raster(interpolate=FALSE)+
geom_tile(data=subset(b,cost<=0.15 & prop_cheaters==0.1 & persr>0), alpha = 0.0, color = "black", size = 1, linejoin = "round")+
geom_tile(data=subset(b,cost<=0.15 & prop_cheaters==0.1 & persr>0),aes(fill=persr),col=NA)+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),strip.background=element_rect(fill=NA,color=NA),
panel.border = element_rect(color = "black", fill = NA, size = 1),legend.position="bottom",legend.text = element_text(angle=45,hjust=1))+
facet_grid(rows=vars(cost),cols=vars(scenario),labeller = label_bquote(rows=Lambda == .(cost)))+
labs(fill="Average effect\non persistence\n\n")+ 
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab(expression(paste("Proportion of innovative cheating (",Psi,")")))+
ggtitle("a",subtitle=expression(paste("10% of cheaters (",bar(Delta) == 0.1,")")))+
scale_fill_gradientn(colors=c("firebrick2","#ffecfb","white","lightsteelblue1","midnightblue"), values=c(0,zero_pos-0.02,zero_pos,zero_pos+0.02,1),
n.breaks=4,limits=c(min(bidon$persr),max(bidon$persr)),labels = scales::percent_format(accuracy=1))+
scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))+scale_y_continuous(breaks=seq(0,1,0.2))+
coord_fixed(ratio=1,expand=F)+guides(fill = guide_colorbar(ticks.colour="black"))


pl2=
ggplot(data=subset(b,cost<=0.15 & prop_cheaters==0.5),aes(x=prop_cheating,y=prop_innovative,fill=persr))+
theme_bw()+
geom_raster(interpolate=FALSE)+
geom_tile(data=subset(b,cost<=0.15 & prop_cheaters==0.5 & persr>0), alpha = 0.0, color = "black", size = 1, linejoin = "round")+
geom_tile(data=subset(b,cost<=0.15 & prop_cheaters==0.5 & persr>0),aes(fill=persr),col=NA)+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),strip.background=element_rect(fill=NA,color=NA),
panel.border = element_rect(color = "black", fill = NA, size = 1),legend.position="none",legend.text = element_text(angle=45,hjust=1))+
facet_grid(rows=vars(cost),cols=vars(scenario),labeller = label_bquote(rows=Lambda == .(cost)))+
labs(fill="Average effect\non persistence")+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab(expression(paste("Proportion of innovative cheating (",Psi,")")))+
ggtitle("b",subtitle=expression(paste("50% of cheaters (",bar(Delta) == 0.5,")")))+
scale_fill_gradientn(colors=c("firebrick2","#ffecfb","white","lightsteelblue1","midnightblue"), values=c(0,zero_pos-0.02,zero_pos,zero_pos+0.02,1),
n.breaks=4,labels = function(x) sprintf("%.2f", x),limits=c(min(bidon$persr),max(bidon$persr)))+
scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))+scale_y_continuous(breaks=seq(0,1,0.2))+
coord_fixed(ratio=1,expand=F)+guides(fill = guide_colorbar(ticks.colour="black"))


leg <- ggpubr::as_ggplot(cowplot::get_legend(pl1))
pl1=pl1+theme(legend.position="none")


pl1=grid.arrange(pl1,right=text_grob("Cost associated with mutualism", size = 10,rot=270,hjust = 0.40,vjust=2),ncol=1)
pl2=grid.arrange(pl2,right=text_grob("Cost associated with mutualism", size = 10,rot=270,hjust = 0.40,vjust=2),ncol=1)

######## HYPERVOLUMES
b=subset(res,connectance==0.4) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(persr=mean(pers_tot-pers0),persr_sde=sd(pers_tot-pers0)/sqrt(length(persr)),
resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T)) 

ini3D(argsAxes3d=list(edges =c('y+', 'x', 'x', 'z')),argsPlot3d=list(xlim=c(0.1,1),zlim=c(0.1,1),ylim =c(0,1)))
bidon=subset(b,cost==0 & persr>0 & scenario=="specialists")
X=bidon$prop_cheating
Y=bidon$prop_innovative
Z=bidon$prop_cheaters
mat=as.matrix(unique(cbind(X,Y,Z)))

ls1=plotHull3D(mat, drawPoints = FALSE,drawLines =FALSE,argsPolygon3d = list(color = "grey",alpha=0.2)) # a line

# bidon=subset(b,cost==0.15 & persr>0)
# X=bidon$prop_cheating
# Y=bidon$prop_innovative
# Z=bidon$prop_cheaters
# mat2=as.matrix(unique(cbind(X,Y,Z)))

# ls2=plotHull3D(mat2, drawPoints = FALSE,drawLines =FALSE,argsPlot3d=list(add=TRUE),argsPolygon3d = list(color = "slategray4",alpha=0.8)) # a line

bidon=subset(b,cost==0.3 & persr>0 & scenario=="specialists")
X=bidon$prop_cheating
Y=bidon$prop_innovative
Z=bidon$prop_cheaters
mat3=as.matrix(unique(cbind(X,Y,Z)))
cxhull(mat3)
ls3=plotHull3D(mat3, drawPoints = FALSE,drawLines =FALSE,argsPlot3d=list(add=TRUE),argsPolygon3d = list(color = "midnightblue",alpha = 1)) # a line

finalize3D(argsTitle3d=list(xlab=expression(Omega),ylab="",zlab=expression(bar(Delta)),cex=2),argsAxes3d=list(edges =c('y+', 'x', 'x', 'z')))
mtext3d(expression(Psi),edge="y+",line=2,las=2,cex=2)

liste=unique(b[,c("cost","scenario")])
liste$volume=NA
for(i in 1:nrow(liste)){
bidon=subset(b,cost==liste$cost[i] & persr>0 & scenario==liste$scenario[i])
X=bidon$prop_cheating
Y=bidon$prop_innovative
Z=bidon$prop_cheaters
mat=as.matrix(unique(cbind(X,Y,Z)))
plotHull3D(mat3, drawPoints = TRUE,drawLines =TRUE) # a line
liste$volume[i]=if(nrow(bidon)>2){cxhull(mat)$volume}else{0}
}

pl3=ggplot(data=liste,aes(x=scenario,fill=as.factor(cost),y=volume))+geom_bar(stat="identity",position=position_dodge(),col="white")+scale_fill_manual(values=c("grey","slategray4","midnightblue"))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
axis.title.x=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
coord_cartesian(expand=F)+labs(fill=expression(paste("Cost (",Lambda,")")))+ylab("Volume of the parameter space in which\ncheating increases persistence")+ggtitle("d")

setwd(dir="C:/Users/Duchenne/Documents/cheating")
td <- image_read("essai.png",density=900) 
pl4 <- image_ggplot(td,interpolate=T)+theme(plot.title=element_text(size=14,face="bold",hjust = 0))+ggtitle("c")


grid.arrange(pl1,pl2,leg,pl3,pl4,layout_matrix =rbind(c(1,2),c(3,3),c(4,5)),widths=c(1,1),heights=c(5,1,4))

setwd(dir="C:/Users/Duchenne/Documents/cheating")
pdf("fig2.pdf",width=8,height=8)
grid.arrange(pl1,pl2,leg,pl4,pl3,layout_matrix =rbind(c(1,2),c(3,3),c(4,5)),widths=c(1,1),heights=c(5,1,3.5))
dev.off();

############ FIGURE S1 ###########
b=subset(res,prop_cheaters>0 & connectance==0.4) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(pers=mean(pers_tot),resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T)) 

pl1=ggplot(data=subset(b,cost<0.6),aes(x=prop_cheating,y=pers,col=prop_innovative,group=paste0(prop_innovative,scenario),linetype=scenario))+geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)), scales="free")+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Network persistence")+scale_y_continuous(breaks=seq(0,1,0.2),labels = scales::percent_format(accuracy=1))+ggtitle("a")+
scale_color_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])+scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))

b=subset(res,connectance==0.4 & prop_cheaters %in% c(0.1,0.3,0.5,0.7,0.9)) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(pers=mean(pers_tot),resilience=mean(-1*valprop,na.rm=T),
feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T))

pl2=ggplot(data=subset(b,prop_innovative %in% c(0,1) & cost<0.6 & prop_cheating<0.5),aes(x=prop_cheating,y=pers,col=prop_innovative,group=paste0(prop_innovative,scenario),linetype=scenario))+
geom_hline(data=subset(b,prop_innovative==0 & cost<0.6 & prop_cheating==0),aes(yintercept=pers),linetype="dashed",color="lightgrey")+
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
pdf("fig_S1.pdf",width=8,height=7)
grid.arrange(pl1,pl2,leg,layout_matrix =rbind(c(1,3),c(2,3)),widths=c(4,1))
dev.off();


############ FIGURE S2 ###########
b=subset(res,prop_cheaters>0 & connectance==0.4 & prop_cheaters %in% c(0.1,0.3,0.5,0.7,0.9)) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(pers_cheaters=mean(nb_cheaters/nb_cheaters_dep,na.rm=T)) 

pls=ggplot(data=b,aes(x=prop_cheating,y=pers_cheaters,col=prop_innovative,group=paste0(prop_innovative,scenario),linetype=scenario))+
geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)))+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Persistence of cheaters")+scale_y_continuous(labels = scales::percent)+ggtitle("")+
scale_color_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])+scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))

png("fig_s2.png",width=1300,height=700,res=150)
pls
dev.off();

############ FIGURE S3 ###########

b=subset(res,connectance==0.4 & prop_cheaters %in% c(0.1,0.3,0.5,0.7,0.9)) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(pers=mean(pers_tot),resilience=mean(-1*valprop,na.rm=T),
feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T))
b=b %>% group_by(cost,prop_cheaters,scenario) %>% mutate(ref=mean(pers[prop_cheating==0]))
b=b %>% group_by(prop_cheating,cost,prop_cheaters,scenario) %>% mutate(crois=max(prop_cheating[(pers-ref)>=0]),gain=pers[which.max(pers)]-ref,opti=prop_innovative[which.max(pers)])
b=b %>% group_by(prop_innovative,cost,prop_cheaters,scenario) %>% mutate(crois=max(prop_cheating[(pers-ref)>=0]))
b=b %>% group_by(cost) %>% mutate(pers_norm=pers/max(pers),maxi=pers[which.max(pers)])
b$cheat_lev=b$crois*b$prop_cheaters

fwrite(b,"theory_results_for_combine.txt")

b2=b %>% group_by(cost,prop_innovative,scenario) %>% summarise(cheat_lev_moy=mean(cheat_lev),cheat_lev_sde=sd(cheat_lev)/sqrt(length(cheat_lev)))

pls2=ggplot(data=b2,aes(x=prop_innovative,y=cheat_lev_moy,linetype=scenario))+
geom_line()+theme_bw()+facet_grid(rows=vars(cost),labeller = label_bquote(rows=Lambda == .(cost)))+
geom_pointrange(aes(ymin=cheat_lev_moy-cheat_lev_sde,ymax=cheat_lev_moy+cheat_lev_sde),linetype=1)+
scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
xlab(expression(paste("Proportion of innovative cheating (",Psi,")")))+ylab("Maximum tolerance to overall level of cheating tolerance")+
labs(color=expression(bar(Delta)))

setwd(dir="C:/Users/Duchenne/Documents/cheating")
png("fig_s3.png",width=1000,height=900,res=150)
pls2
dev.off();


############ FIGURE S4 ###########

b=subset(res,prop_cheaters %in% c(0.1) & connectance==0.4 & nb_cheaters<(nbsp_a) & nbsp_p>1 & resilience>0) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(pers=mean(pers_tot),
resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(cheaters_to_poll,na.rm=T),n=length(cheaters_to_poll[!is.na(cheaters_to_poll)]),contrib2=mean(cheaters_to_plant,na.rm=T))

pl1=ggplot(data=subset(b,n>0),aes(x=prop_cheating,y=contrib,col=prop_innovative,group=paste0(prop_innovative,scenario),linetype=scenario))+geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)),scales="free")+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Average total effect of cheaters on other pollinators")+ggtitle("")+
scale_color_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])+scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1),lim=c(0,1))

pl2=ggplot(data=subset(b,n>0),aes(x=prop_cheating,y=contrib2,col=prop_innovative,group=paste0(prop_innovative,scenario),linetype=scenario))+geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)),scales="free")+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Average total effect of cheaters on other pollinators")+ggtitle("")+
scale_color_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])+scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1),lim=c(0,1))

sum(b$n)

grid.arrange(pl1,pl2,ncol=2)

setwd(dir="C:/Users/Duchenne/Documents/cheating")
png("fig_S4.png",width=800,height=800,res=150)
pl1
dev.off();

############ FIGURE S5 ###########

res$NODF_c="intermediate"
res$NODF_c[res$NODF_Iini>quantile(res$NODF_Iini[res$connectance==0.4],prob=0.75)]="high"
res$NODF_c[res$NODF_Iini<quantile(res$NODF_Iini[res$connectance==0.4],prob=0.25)]="low"
b=subset(res,prop_cheaters %in% c(0.1,0.3,0.5,0.7,0.9) & connectance==0.4 & prop_innovative==1) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,NODF_c,scenario) %>%
summarise(pers=mean(pers_tot),pers_eff=mean(pers_tot)-mean(pers0))

pl1=ggplot(data=b,aes(x=prop_cheating,y=pers_eff,col=NODF_c,group=paste0(NODF_c,scenario),linetype=scenario))+
geom_hline(yintercept=0,linetype="dashed",color="lightgrey")+
geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)),scales="free")+
labs(color="NODF")+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Network presistence, relative to a case without cheating")+ggtitle("")+
scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))+scale_color_manual(values=c("firebrick3","deeppink2","pink"))

setwd(dir="C:/Users/Duchenne/Documents/cheating")
png("fig_S5.png",width=1300,height=1300,res=170)
pl1
dev.off();


############ FIGURE S6 ###########
b=subset(res,prop_cheaters>0 & connectance==0.2 & prop_cheaters %in% c(0.1,0.3,0.5,0.7,0.9)) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% 
summarise(pers=mean(pers_tot),resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T)) 

pl1=ggplot(data=b,aes(x=prop_cheating,y=pers,col=prop_innovative,group=paste0(prop_innovative,scenario),linetype=scenario))+geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
facet_grid(rows=vars(cost),cols=vars(prop_cheaters),labeller = label_bquote(cols=bar(Delta) == .(prop_cheaters),rows=Lambda == .(cost)), scales="free")+
labs(color=expression(Psi))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab("Network persistence")+scale_y_continuous(breaks=seq(0,1,0.2),labels = scales::percent_format(accuracy=1))+ggtitle("a")+
scale_color_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])+scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))


b=subset(res,prop_cheaters>0 & connectance==0.2 & prop_cheaters %in% c(0.1,0.3,0.5,0.7,0.9)) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>%
summarise(pers=mean(pers_tot),resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T)) 

pl2=ggplot(data=subset(b,prop_innovative %in% c(0,1) & prop_cheating<0.5),aes(x=prop_cheating,y=pers,col=prop_innovative,group=paste0(prop_innovative,scenario),linetype=scenario))+
geom_hline(data=subset(b, prop_innovative==0 & prop_cheating==0),aes(yintercept=pers),linetype="dashed",color="lightgrey")+
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
png("fig_S6.png",width=1300,height=1400,res=150)
grid.arrange(pl1,pl2,leg,layout_matrix =rbind(c(1,3),c(2,3)),widths=c(4,1))
dev.off();

####################
######## FIG 3 ##########
res$resiliencer=res$pers_tot-res$pers0
res=res %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% mutate(pers_moy=mean(pers_tot),pers_sd=sd(pers_tot))
1-(var(res$pers_tot-res$pers_moy))/var(res$pers_tot)
res$resi=(res$pers_tot-res$pers_moy)
res$resi2=round(res$resi,digits=3)

res=res %>% group_by(connectance) %>% mutate(NODF_moy=mean(NODF_Iini,na.rm=T),mod_moy=mean(mod_Iini,na.rm=T),resi_moy=mean(resi),resi_sd=sd(resi))
res$NODF_Iini2=res$NODF_Iini-res$NODF_moy
res$mod_Iini2=res$mod_Iini-res$mod_moy
res$C_Iini2=res$C_Iini


cor(res[,c("mod_Iini","NODF_Iini","connectance","NODF_Iini2","mod_Iini2","C_Iini2")])
model=lm(resi~avg_ra+avg_rp+mod_Iini2+NODF_Iini2+C_Iini2,data=res)
summary(model)
ano=Anova(model)
vec=(ano[,1]/sum(ano[,1]))
sum(vec[1:2])
sum(vec[3:5])
sum(vec[6])


newdat=data.frame(NODF_Iini2=mean(res$NODF_Iini2),C_Iini2=0.3,
avg_rp=c(seq(min(res$avg_rp),max(res$avg_rp),length.out=100),rep(mean(res$avg_rp),each=100)),
avg_ra=c(rep(mean(res$avg_ra),each=100),seq(min(res$avg_ra),max(res$avg_ra),length.out=100)),
guild=rep(c("plants","animals"),each=100),x=seq(min(res$avg_ra),max(res$avg_ra),length.out=100),mod_Iini2=mean(res$mod_Iini2))
pred=cbind(newdat,as.data.frame(predict(model,interval ="confidence",newdata=newdat)))

pl1=ggplot(data=pred,aes(x=x,color=as.factor(guild),fill=as.factor(guild),y=fit,ymax=upr,ymin=lwr))+
geom_ribbon(alpha=0.3)+
geom_line()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0))+ggtitle("a")+ylab("Residual effect of cheating on persistence")+xlab("Average growth rate")+scale_color_manual(values=c("#D0CD94","#3C787E"))+
scale_fill_manual(values=c("#D0CD94","#3C787E"))+labs(color="",fill="")+ scale_alpha(guide = 'none')

newdat=data.frame(NODF_Iini2=seq(min(res$NODF_Iini2),max(res$NODF_Iini2),length.out=100),C_Iini2=rep(c(0.2,0.3,0.4),each=100),avg_rp=mean(res$avg_rp),avg_ra=mean(res$avg_ra),mod_Iini2=mean(res$mod_Iini2))
pred=cbind(newdat,as.data.frame(predict(model,interval ="confidence",newdata=newdat)))

pl2=ggplot()+
geom_ribbon(data=pred,aes(x=NODF_Iini2,ymax=upr,ymin=lwr,fill=as.factor(C_Iini2)),alpha=0.4)+
geom_line(data=pred,aes(x=NODF_Iini2,y=fit,color=as.factor(C_Iini2)),size=1)+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none")+ggtitle("b")+ylab("Residual effect of cheating on persistence")+xlab("Corrected nestedness (NODF)")+scale_alpha(range = c(0.001, 0.1))+
scale_color_manual(values=c("lightpink","hotpink","deeppink4"))+
scale_fill_manual(values=c("lightpink","hotpink","deeppink4"))+labs(color=expression(phi1),fill=expression(phi1))+ scale_alpha(guide = 'none')

newdat=data.frame(NODF_Iini2=mean(res$NODF_Iini2),C_Iini2=rep(c(0.2,0.3,0.4),each=100),avg_rp=mean(res$avg_rp),avg_ra=mean(res$avg_ra),mod_Iini2=seq(min(res$mod_Iini2),max(res$mod_Iini2),length.out=100))
pred=cbind(newdat,as.data.frame(predict(model,interval ="confidence",newdata=newdat)))

pl3=ggplot()+
geom_ribbon(data=pred,aes(x=mod_Iini2,ymax=upr,ymin=lwr,fill=as.factor(C_Iini2)),alpha=0.4)+
geom_line(data=pred,aes(x=mod_Iini2,y=fit,color=as.factor(C_Iini2)),size=1)+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="right",axis.title.y=element_blank())+ggtitle("c")+ylab("Residual effect of cheating on persistence")+xlab("Corrected modularity")+scale_alpha(range = c(0.001, 0.1))+
scale_color_manual(values=c("lightpink","hotpink","deeppink4"))+
scale_fill_manual(values=c("lightpink","hotpink","deeppink4"))+labs(color=expression(phi1),fill=expression(phi1))+ scale_alpha(guide = 'none')

grid.arrange(pl1,pl2,pl3,ncol=3,widths=c(3,2.6,3))

setwd(dir="C:/Users/Duchenne/Documents/cheating")
pdf("fig3.pdf",width=9,height=3)
grid.arrange(pl1,pl2,pl3,ncol=3,widths=c(3,2.6,3))
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
