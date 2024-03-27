library(data.table)
library(dplyr)
setwd(dir="C:/Users/Duchenne/Documents/cheating/eq_review")

lili=list.files()
lili=lili[grep("netcar",lili)]
res=NULL
for(i in lili){
#bidon=fread(paste0("netcar_",i,".txt"))
bidon=fread(i)
bidon=bidon[,c("prop_cheating","prop_innovative","cost","prop_cheaters","scenario","connectance","essai","pers_tot","nbsp_a","nbsp_p","nb_cheaters","nb_cheaters_dep","valprop","C_Iini","NODF_Iini","mod_Iini",
"ai_contrib_aa","cheaters_to_poll","cheaters_to_plant","nbsp_p_dep","nbsp_a_dep","interfp"),with=F]
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
setwd(dir="C:/Users/Duchenne/Documents/cheating/eq_review")

res=fread("data_for_analyse.txt")

res$feas=0
res$feas[res$pers_tot==1]=1
res$resilience=-1*res$valprop
nrow(subset(res,resilience<0))
res=res %>% group_by(prop_innovative,cost,prop_cheaters,scenario,connectance,essai,nbsp_a_dep,interfp) %>%
mutate(pers0=pers_tot[prop_cheating==0],resilience0=resilience[prop_cheating==0])


####poster & presentation stuffs
res$overall_level_cheating=res$prop_cheaters*res$prop_cheating
b=subset(res,prop_cheaters>0) %>% group_by(overall_level_cheating,prop_innovative,cost,scenario,nbsp_a_dep,interfp) %>%
summarise(persr=mean(pers_tot-pers0),persr_sde=sd(pers_tot-pers0)/sqrt(length(persr)),
resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T)) 

summary(b$persr)

bidon=subset(b,cost<=0.15 & nbsp_a_dep==20 & interfp==1)
bidon$overall_level_cheating2=plyr::round_any(bidon$overall_level_cheating,0.1)
bidon$overall_level_cheating2[bidon$overall_level_cheating>0 & bidon$overall_level_cheating2==0]=0.1
zero_pos=unique(scales::rescale(bidon$persr, to = c(0, 1))[bidon$persr==0])

cost_plot=0

setwd(dir="C:/Users/Duchenne/Documents/cheating")
png("figpres.png",width=900,height=900,res=200)
ggplot(data=subset(bidon,cost==cost_plot & nbsp_a_dep==20),aes(x=overall_level_cheating2,y=prop_innovative,fill=persr))+
theme_bw()+
geom_raster(interpolate=FALSE)+
geom_tile(data=subset(bidon,persr>0 & cost==cost_plot & nbsp_a_dep==20 ), alpha = 0.0, color = "black", size = 1, linejoin = "round")+
geom_tile(data=subset(bidon,persr>0 & cost==cost_plot & nbsp_a_dep==20),aes(fill=persr),col=NA)+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),strip.background=element_rect(fill=NA,color=NA),
panel.border = element_rect(color = "black", fill = NA, size = 1),legend.position="top",legend.text = element_text(angle=45,hjust=1))+
facet_grid(cols=vars(scenario))+
labs(fill="Average effect\non persistence\n\n")+ 
xlab(expression(paste("Overall level of cheating")))+ylab(expression(paste("Proportion of innovative cheating ( ",Psi,")")))+
scale_fill_gradientn(colors=c("firebrick2","#ffecfb","white","lightsteelblue1","midnightblue"), values=c(0,zero_pos-0.02,zero_pos,zero_pos+0.02,1),
n.breaks=4,limits=c(min(bidon$persr),max(bidon$persr)),labels = scales::percent_format(accuracy=1))+
scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))+scale_y_continuous(breaks=seq(0,1,0.2))+
coord_fixed(ratio=1,expand=F)+guides(fill = guide_colorbar(ticks.colour="black"))
dev.off();

############ FIGURE 2 ###########
b=subset(res,prop_cheaters>0) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario,nbsp_a_dep,interfp) %>%
summarise(persr=mean(pers_tot-pers0),persr_sde=sd(pers_tot-pers0)/sqrt(length(persr)),
resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T)) 

summary(b$persr)

bidon=subset(b,cost<=0.15 & prop_cheaters %in% c(0.1,0.5) & nbsp_a_dep==20 & interfp==1)
zero_pos=unique(scales::rescale(bidon$persr, to = c(0, 1))[bidon$persr==0])

pl1=
ggplot(data=subset(bidon, prop_cheaters==0.1),aes(x=prop_cheating,y=prop_innovative,fill=persr))+
theme_bw()+
geom_raster(interpolate=FALSE)+
geom_tile(data=subset(bidon,prop_cheaters==0.1 & persr>0), alpha = 0.0, color = "black", size = 1, linejoin = "round")+
geom_tile(data=subset(bidon,prop_cheaters==0.1 & persr>0),aes(fill=persr),col=NA)+
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
ggplot(data=subset(bidon, prop_cheaters==0.5),aes(x=prop_cheating,y=prop_innovative,fill=persr))+
theme_bw()+
geom_raster(interpolate=FALSE)+
geom_tile(data=subset(bidon,prop_cheaters==0.5 & persr>0), alpha = 0.0, color = "black", size = 1, linejoin = "round")+
geom_tile(data=subset(bidon,prop_cheaters==0.5 & persr>0),aes(fill=persr),col=NA)+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),strip.background=element_rect(fill=NA,color=NA),
panel.border = element_rect(color = "black", fill = NA, size = 1),legend.position="bottom",legend.text = element_text(angle=45,hjust=1))+
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
b=subset(res,connectance==0.4) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario,interfp,nbsp_a_dep) %>%
summarise(persr=mean(pers_tot-pers0),persr_sde=sd(pers_tot-pers0)/sqrt(length(persr)),
resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T)) 

ini3D(argsAxes3d=list(edges =c('y+', 'x', 'x', 'z')),argsPlot3d=list(xlim=c(0.1,1),zlim=c(0.1,1),ylim =c(0,1)))
bidon=subset(b,cost==0 & persr>0 & nbsp_a_dep==20 & interfp==1 & scenario=="specialists")
X=bidon$prop_cheating
Y=bidon$prop_innovative
Z=bidon$prop_cheaters
mat=as.matrix(unique(cbind(X,Y,Z)))

ls1=plotHull3D(mat, drawPoints = FALSE,drawLines =FALSE,argsPolygon3d = list(color = "grey",alpha=0.2)) # a line

bidon=subset(b,cost==0.3 & persr>0 & scenario=="specialists" & nbsp_a_dep==20 & interfp==1)
X=bidon$prop_cheating
Y=bidon$prop_innovative
Z=bidon$prop_cheaters
mat3=as.matrix(unique(cbind(X,Y,Z)))
cxhull(mat3)
ls3=plotHull3D(mat3, drawPoints = FALSE,drawLines =FALSE,argsPlot3d=list(add=TRUE),
argsPolygon3d = list(color = "midnightblue",alpha = 1)) # a line

finalize3D(argsTitle3d=list(xlab=expression(Omega),ylab="",zlab=expression(bar(Delta)),cex=2),
argsAxes3d=list(edges =c('y+', 'x', 'x', 'z')))
mtext3d(expression(Psi),edge="y+",line=2,las=2,cex=2)

liste=unique(b[,c("cost","scenario","nbsp_a_dep")])
liste$volume=NA
for(i in 1:nrow(liste)){
bidon=subset(b,cost==liste$cost[i] & persr>0 & nbsp_a_dep==liste$nbsp_a_dep[i] & interfp==1 & scenario==liste$scenario[i])
X=bidon$prop_cheating
Y=bidon$prop_innovative
Z=bidon$prop_cheaters
mat=as.matrix(unique(cbind(X,Y,Z)))
#plotHull3D(mat, drawPoints = TRUE,drawLines =TRUE) # a line
liste$volume[i]=if(length(unique(X))>2 & length(unique(Y))>2 & length(unique(Z))>2){cxhull(mat)$volume}else{0}
}

pl3=ggplot(data=subset(liste,nbsp_a_dep==20),aes(x=scenario,fill=as.factor(cost),y=volume))+geom_bar(stat="identity",position=position_dodge(),col="white")+scale_fill_manual(values=c("grey","slategray4","midnightblue"))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
axis.title.x=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
coord_cartesian(expand=F)+labs(fill=expression(paste("Cost (",Lambda,")")))+ylab("Volume of the parameter space in which\ncheating increases persistence")+ggtitle("d")

setwd(dir="C:/Users/Duchenne/Documents/cheating")
td <- image_read("essai2.png",density=900) 
pl4 <- image_ggplot(td,interpolate=T)+theme(plot.title=element_text(size=14,face="bold",hjust = 0))+ggtitle("c")

grid.arrange(pl1,pl2,leg,pl3,pl4,layout_matrix =rbind(c(1,2),c(3,3),c(4,5)),widths=c(1,1),heights=c(5,1,4))

setwd(dir="C:/Users/Duchenne/Documents/cheating")
pdf("fig2.pdf",width=8,height=8)
grid.arrange(pl1,pl2,leg,pl4,pl3,layout_matrix =rbind(c(1,2),c(3,3),c(4,5)),widths=c(1,1),heights=c(5,1,3.5))
dev.off();


compar=dcast(data=b,prop_cheating+prop_innovative+cost+prop_cheaters+scenario+interfp~paste0("d",nbsp_a_dep),value.var="persr")
compar=subset(compar,!is.na(d10))

pl5=ggplot(data=compar,aes(x=d20,y=d10))+
geom_rect(aes(xmin=-0.06, xmax = max(compar$d20,na.rm=T)+0.005,ymin = min(subset(compar,d20>-0.06 & d20<max(compar$d20,na.rm=T)+0.005)$d10),
ymax = max(subset(compar,d20>-0.06 & d20<max(compar$d20,na.rm=T)+0.005)$d10)+0.025),fill="grey93",color="red")+
geom_abline(intercept=0,slope=1)+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
geom_point(alpha=0.5)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0))+
coord_fixed(ratio=1,expand=T)+facet_zoom(xlim = c(-0.05, max(compar$d20[!is.na(compar$d10)],na.rm=T)),
ylim=c(min(compar$d10[compar$d20>=(-0.05)],na.rm=T),max(compar$d10[compar$d20>=(-0.05)],na.rm=T)),
horizontal=F,shrink=TRUE,zoom.size=1,
show.area=F)+ggtitle("b")+xlab(expression(paste("Effect of cheating on persistence when ",n[sp]==40)))+
ylab(expression(paste("Effect of cheating on persistence when ",n[sp]==20)))

liste$nbsp_dep=liste$nbsp_a_dep*2
liste$nbsp_dep=factor(liste$nbsp_dep,levels=c("40","20"))
pl6=ggplot(data=subset(liste,cost==0.15),aes(x=scenario,fill=nbsp_dep,y=volume))+
geom_bar(stat="identity",position=position_dodge(),col="white")+
scale_fill_manual(values=c("black","grey"))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
axis.title.x=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
coord_cartesian(expand=F)+labs(fill=expression(n[sp]))+
ylab("Volume of the parameter space in which\ncheating increases persistence")+ggtitle("c")


bidon=subset(b,cost==0.15 & prop_cheaters==0.1 & nbsp_a_dep==20 & scenario=="specialists")
zero_pos=unique(scales::rescale(bidon$persr, to = c(0, 1))[bidon$persr==0])
bidon$interf="no competition among\npollinators for partners"
bidon$interf[bidon$interfp==1]="competition among\npollinators for partners"

pl7=
ggplot(data=bidon,aes(x=prop_cheating,y=prop_innovative,fill=persr))+
theme_bw()+
geom_raster(interpolate=FALSE)+
geom_tile(data=subset(bidon,persr>0), alpha = 0.0, color = "black", size = 1, linejoin = "round")+
geom_tile(data=subset(bidon,persr>0),aes(fill=persr),col=NA)+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
panel.border = element_rect(color = "black", fill = NA, size = 1),
legend.position="bottom",legend.text = element_text(angle=45,hjust=1))+
facet_grid(rows=vars(cost),cols=vars(interf),labeller = label_bquote(rows=Lambda == .(cost)))+
labs(fill="Average effect\non persistence\n\n")+ 
xlab(expression(paste("Cheating frequency (",Omega,")")))+
ylab(expression(paste("Proportion of innovative cheating (",Psi,")")))+
ggtitle("a",subtitle=expression(paste("10% of cheaters (",bar(Delta) == 0.1,")")))+
scale_fill_gradientn(colors=c("firebrick2","#ffecfb","white","lightsteelblue1","midnightblue"),
values=c(0,zero_pos-0.02,zero_pos,zero_pos+0.02,1),
n.breaks=4,limits=c(min(bidon$persr),max(bidon$persr)),labels = scales::percent_format(accuracy=1))+
scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))+scale_y_continuous(breaks=seq(0,1,0.2))+
coord_fixed(ratio=1,expand=F)+guides(fill = guide_colorbar(ticks.colour="black"))

setwd(dir="C:/Users/Duchenne/Documents/cheating")
pdf("fig3.pdf",width=10,height=8)
grid.arrange(pl7,pl6,pl5,layout_matrix=rbind(c(1,3),c(2,3)),widths=c(1.2,1),heights=c(1.2,1))
dev.off();


############ FIGURE S1 ###########
b=subset(res,prop_cheaters>0 & connectance==0.4  & prop_cheaters %in% c(0.1,0.3,0.5,0.7,0.9)) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>%
summarise(pers=mean(pers_tot),resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T)) 

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
b=subset(res,prop_cheaters>0 & connectance==0.4 & prop_cheaters %in% c(0.1,0.3,0.5,0.7,0.9)) %>%
group_by(prop_cheating,prop_innovative,cost,prop_cheaters,scenario) %>% summarise(pers_cheaters=mean(nb_cheaters/nb_cheaters_dep,na.rm=T)) 

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


############ FIGURE S5 ###########
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
png("fig_S5.png",width=1200,height=400,res=120)
grid.arrange(pl1,pl2,pl3,ncol=3,widths=c(3,2.6,3))
dev.off();

############ FIGURE S4 ###########
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
png("fig_S4.png",width=1300,height=1400,res=150)
grid.arrange(pl1,pl2,leg,layout_matrix =rbind(c(1,3),c(2,3)),widths=c(4,1))
dev.off();






























