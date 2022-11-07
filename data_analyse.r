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
library(ggplotify)
library(ggtern)
library(ks)
setwd(dir="C:/Users/Duchenne/Documents/cheating")

dat=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_2022-11-07/Interactions_data_Costa-Rica.txt")
dat$year=year(dat$date)

cameras=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_2022-11-07/Cameras_data_Costa-Rica.txt",na.strings = c("",NA))
cameras$month=month(cameras$start_date)
cameras$year=year(cameras$start_date)
dim(dat)
dat=merge(dat,cameras,by=c("waypoint","site"),all.x=T,all.y=F) #if you want to exclude camera which did not detect any hummingbird
dim(dat)

getmode <- function(v) {
 uniqv <- unique(na.omit(v))
 uniqv[which.max(tabulate(match(v, uniqv)))]
}

site=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_2022-11-07/Site_metadata_Costa-Rica.txt")


tr=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/plant_traits_2022-11-07/Plant_traits.txt",na.strings = c("",NA))
tr2=tr %>% group_by(plant_species,Symmetry,Type,AttractivePart) %>% summarise(Tube_length=mean(Tube_length,na.rm=T),Anther_length=mean(Anther_length,na.rm=T), Stigma_length=mean(Stigma_length,na.rm=T),
Opening_corrolla=mean(Opening_corrolla,na.rm=T),Curvature_middle=mean(Curvature_middle,na.rm=T))
tr2[tr2$plant_species %in% tr2$plant_species[duplicated(tr2$plant_species)],]
tr2=subset(tr,!is.na(plant_species))  %>% group_by(plant_species,plant_family,plant_genus) %>% summarise(Tube_length=mean(Tube_length,na.rm=T),Anther_length=mean(Anther_length,na.rm=T), Stigma_length=mean(Stigma_length,na.rm=T),
Opening_corrolla=mean(Opening_corrolla,na.rm=T),Curvature_middle=mean(Curvature_middle,na.rm=T),Type=getmode(Type))


###### PLOT NETWORKS:

bi=dat %>% group_by(plant_species,site) %>% mutate(duration_plant=sum(duration_sampling_hours[!duplicated(waypoint)],na.rm=T))
b=bi %>% group_by(plant_species,hummingbird_species,site) %>% summarise(nb_inter=length(date[piercing=="no"]),duration=unique(duration_plant),nb_cheat=length(date[piercing=="yes"]))
b$value=(b$nb_cheat+b$nb_inter)/(b$duration/24)
b$cheating=b$nb_cheat/(b$nb_cheat+b$nb_inter)
b=subset(b,!is.na(value) & !is.na(plant_species) & !is.na(hummingbird_species))
rbPal <- colorRampPalette(c("black",'red'))
b$Col <- rbPal(10)[as.numeric(cut(b$cheating,breaks = 10))]
png("figS6.png",width=1600,height=1400,res=120)
par(mfrow=c(3,4))
for(i in unique(dat$site)){
bidon=subset(b,site==i)
mat=dcast(bidon,plant_species~hummingbird_species,value.var="value",fill=0)

plotweb(sqrt(mat[,-1]),col.interaction=b$Col,col.high="dodgerblue3",col.low="forestgreen",bor.col.interaction=NA,high.lablength=0,low.lablength=0,bor.col.high="dodgerblue3",bor.col.low="forestgreen")
title(paste0(i))
}
dev.off();

bidon=subset(b,site=="TOLO")
mat=dcast(bidon,plant_species~hummingbird_species,value.var="value",fill=0)

gr1=plot_grid(base2grob(~plotweb(sqrt(mat[,-1]),col.interaction=b$Col,col.high="dodgerblue3",col.low="forestgreen",bor.col.interaction=NA,high.lablength=0,low.lablength=0,bor.col.high="dodgerblue3",bor.col.low="forestgreen")))+
ggtitle("a")+theme(plot.title=element_text(size=14,face="bold",hjust = 0))

#### HYPOTHESE ONE AND TWO: CHEATERS ARE SPECIALISTS
dat=dat %>% group_by(plant_species,site) %>% mutate(duration_plant=sum(duration_sampling_hours[!duplicated(waypoint)],na.rm=T))
b=b=subset(dat,!is.na(duration_plant) & !is.na(plant_species) & !is.na(hummingbird_species)) %>% group_by(plant_species,hummingbird_species,site) %>% summarise(nb_inter=length(date[piercing=="no"]),duration=unique(duration_plant),nb_cheat=length(date[piercing=="yes"]))
b$value=(b$nb_cheat+b$nb_inter)/(b$duration/24)
b$mut=(b$nb_inter)/(b$duration/24)
b$cheating=b$nb_cheat/(b$nb_cheat+b$nb_inter)
b$cheat=(b$nb_cheat)/(b$duration/24)
b=b %>% group_by(hummingbird_species) %>% mutate(cheating_moy=mean((cheat/(mut+cheat)),na.rm=T),div=length(unique(plant_species[mut>0.1])),total_inter=sum(cheat+mut))

b2=b %>% group_by(hummingbird_species,site) %>% summarise(div=length(unique(plant_species[mut>0.5])),cheat=sum(cheat)/(sum(cheat)+sum(mut)),div2=vegan::diversity(mut),total_inter=sum(cheat)+sum(mut))

bsup=b2 %>% group_by(site) %>% summarise(prop_cheaters=length(hummingbird_species[cheat>0.05])/length(hummingbird_species))
bsup=merge(bsup,site,by="site")

model=glm(prop_cheaters~min_transect_elev,family="quasibinomial",data=bsup)
summary(model)
pre=ggpredict(model,"min_transect_elev[600:3110]")

pl1=ggplot()+geom_point(data=bsup,aes(x=min_transect_elev,y=prop_cheaters))+
theme_bw()+
geom_ribbon(data=pre,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.2)+geom_line(data=pre,aes(x=x,y=predicted))+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="none",panel.border = element_blank(),axis.line= element_line(),
plot.title=element_text(size=14,face="bold",hjust = 0))+ylab("Proportion of cheaters")+
xlab("Elevation (m)")+scale_y_continuous(breaks=c(0,0.5,1),labels=c(0,0.5,1),limits = c(0,1))+
ggtitle("b")

model=glmer(cheat~div+(1|site),family="binomial",data=b2,weight=total_inter)
summary(model)
pre=ggpredict(model,"div")

pl2=ggplot()+geom_point(data=b2,aes(x=div,y=cheat,size=total_inter))+geom_ribbon(data=pre,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.2)+geom_line(data=pre,aes(x=x,y=predicted))+
theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",panel.border = element_blank(),axis.line= element_line())+ggtitle("c")+
xlab("Number of partners")+ylab("Frequency of cheating")+scale_y_continuous(breaks=c(0,0.5,1),labels=c(0,0.5,1),limits = c(0,1))
#scale_color_manual(values=c("red","#D95F02","#7570B3","#E6AB02","#A6761D","#666666",rgb(219/620,194/620,207/620,1),"#9fa2b2","#3c7a89","#2e4756","#16262e","#0a210f","#14591d","#99aa38","#068D9D","#acd2ed"))


#### HYPOTHESE THREE: CHEATING IS INNOVATIVE
b3=b %>% group_by(hummingbird_species,plant_species,cheating_moy,div,total_inter) %>% summarise(dens_mut=sum(mut),dens_cheat=sum(cheat))
b3=b3 %>% group_by(hummingbird_species,cheating_moy,div,total_inter) %>% mutate(some_mut=sum(dens_mut),some_cheat=sum(dens_cheat))

#INFER missing tube_length
bidon=subset(tr2,!is.na(Anther_length) & !is.na(Stigma_length) & !is.na(Type) & !is.na(Tube_length))
model_tube=randomForest(Tube_length~plant_family+plant_genus+Type+Anther_length+Stigma_length,data=bidon)
tr2$Tube_length_pre=predict(model_tube,newdata=tr2)
plot(Tube_length_pre~Tube_length,data=tr2)
tr2$Tube_length_pre[!is.na(tr2$Tube_length)]=tr2$Tube_length[!is.na(tr2$Tube_length)]
#INFER missing curvature
bidon=subset(tr2,!is.na(Anther_length) & !is.na(Stigma_length) & !is.na(Type)  & !is.na(Curvature_middle))
model_curv=randomForest(Curvature_middle~plant_family+Type+plant_genus+Anther_length+Stigma_length,data=bidon)
tr2$Curvature_middle_pre=predict(model_curv,newdata=tr2)
plot(Curvature_middle_pre~Curvature_middle,data=tr2)
tr2$Curvature_middle_pre[!is.na(tr2$Curvature_middle)]=tr2$Curvature_middle[!is.na(tr2$Curvature_middle)]

b3=merge(b3,tr2,by=c("plant_species"),all.x=T,all.y=F) #if you want to exclude camera which did not detect any hummingbird

list_cheaters=unique(subset(b3,some_cheat>0 & !is.na(Tube_length))$hummingbird_species)


#ONE DIMENSION
resf=NULL
for(i in 1:length(list_cheaters)){
bidon=subset(b3,hummingbird_species==list_cheaters[i] & !is.na(Tube_length))
if(nrow(bidon)>1){
weight1=bidon$dens_mut/bidon$some_mut
weight1[is.na(weight1)]=0
weight2=bidon$dens_cheat/bidon$some_cheat
weight2[is.na(weight2)]=0
resi=data.frame(hummingbird_species=list_cheaters[i],dens_mut=density(bidon$Tube_length,weight=weight1,from=0, to=max(b3$Tube_length,na.rm=T))$y,
dens_cheat=density(bidon$Tube_length,weight=weight2,from=0, to=max(b3$Tube_length,na.rm=T))$y,
Tube_length=density(bidon$Tube_length,weight=weight2,from=0, to=max(b3$Tube_length,na.rm=T))$x,cheating_moy=unique(bidon$cheating_moy),div=unique(bidon$div),total_inter=unique(bidon$total_inter))
resf=rbind(resf,resi)
}
}

resf$mini=apply(resf[,c("dens_mut","dens_cheat")],1,min)
b4=resf %>% group_by(hummingbird_species,cheating_moy,div,total_inter) %>% summarise(prop_innovative=1-sum(mini)*(resf$Tube_length[2]-resf$Tube_length[1]))

ggplot(data=b4,aes(x=1,y=prop_innovative))+geom_boxplot()+geom_point(alpha=0.5,aes(color=cheating_moy,size=total_inter))+scale_color_gradientn(colors=c("black","firebrick3"))
ggplot(data=resf)+geom_line(aes(x=Tube_length,y=dens_mut_pred))+geom_line(aes(x=Tube_length,y=dens_cheat),col="firebrick3")+facet_wrap(~hummingbird_species)


#TWO DIMENSIONS
resf=NULL
list_plot=list()
for(i in 1:length(list_cheaters)){
bidon=subset(b3,hummingbird_species==list_cheaters[i] & !is.na(Tube_length_pre)  & !is.na(Curvature_middle_pre))
if(nrow(bidon)>2){
weight1=bidon$dens_mut/bidon$some_mut
weight1[is.na(weight1)]=0
weight2=bidon$dens_cheat/bidon$some_cheat
weight2[is.na(weight2)]=0


#calculate bandwidth to use
mati=Hpi(bidon[,c("Tube_length_pre","Curvature_middle_pre")])
### Kernel functional estimate for 1- to 6-dimensional data:
obj=kde(bidon[,c("Tube_length_pre","Curvature_middle_pre")],
w=weight1,gridsize=50,H=mati,xmin=c(0,0),xmax=c(max(b3$Tube_length_pre,na.rm=T),max(b3$Curvature_middle_pre,na.rm=T)))
image2D(obj$estimate)
resi=melt(obj$estimate)
names(resi)[1:3]=c("Tube_length","Curvature_middle","dens_mut")
resi$Tube_length=obj$eval.points[[1]]
resi$Curvature_middle=rep(obj$eval.points[[2]],each=length(obj$eval.points[[1]]))

obj=kde(bidon[,c("Tube_length_pre","Curvature_middle_pre")],
w=weight2,gridsize=50,H=mati,xmin=c(0,0),xmax=c(max(b3$Tube_length_pre,na.rm=T),max(b3$Curvature_middle_pre,na.rm=T)))
resi$dens_cheat=melt(obj$estimate)$value
resi$hummingbird_species=list_cheaters[i]
resi$cheating_moy=unique(bidon$cheating_moy)
resi$div=unique(bidon$div)
resi$total_inter=unique(bidon$total_inter)


list_plot[[(length(list_plot)+1)]]=ggplot()+geom_tile(data=resi,aes(x=Tube_length,y=Curvature_middle,fill=dens_mut))+scale_fill_gradientn(colours = terrain.colors(10))+
geom_point(data=bidon,aes(x=Tube_length_pre,y=as.numeric(Curvature_middle_pre),size=dens_cheat,alpha=0.6))+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="none",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=12))+xlab("")+
ylab("")+coord_cartesian(expand = FALSE)

#check:
sum(resi$dens_mut)*(resi$Tube_length[2]-resi$Tube_length[1])*(resi$Curvature_middle[2]-resi$Curvature_middle[1])

resf=rbind(resf,resi)
}
}

list_plot = list_plot[lapply(list_plot,length)>0]
grid.arrange(grobs = list_plot , ncol = 4,bottom ="Corolla length (cm)",left="Curvature measure")

resf$mini=apply(resf[,c("dens_mut","dens_cheat")],1,min)
b4=resf %>% dplyr::group_by(hummingbird_species,cheating_moy,div,total_inter) %>% dplyr::summarise(prop_innovative=1-sum(mini)*(Tube_length[2]-Tube_length[1])*(obj$eval.points[[2]][2]-obj$eval.points[[2]][1]))

pl3=ggplot(data=b4,aes(x=1,y=prop_innovative))+geom_violin(trim=F)+geom_boxplot(width=0.2)+geom_jitter(alpha=0.5,aes(color=cheating_moy,size=total_inter),width=0.02)+scale_color_gradientn(colors=c("black","firebrick3"))+
theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="none",panel.border = element_blank(),,axis.line= element_line(),
plot.title=element_text(size=14,face="bold",hjust = 0),
axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_text())+ylab("Proportion of innovative cheating")+
xlab("\n")+scale_y_continuous(breaks=c(0,0.5,1),labels=c(0,0.5,1),limits = c(0,1))+
ggtitle("d")


#### HYPOTHESE THREE: CHEATING FREQUENCY IS LOW
dat=dat %>% group_by(plant_species,site) %>% mutate(duration_plant=sum(duration_sampling_hours[!duplicated(waypoint)],na.rm=T))
b=b=subset(dat,!is.na(duration_plant) & !is.na(plant_species) & !is.na(hummingbird_species)) %>% group_by(plant_species,hummingbird_species,site) %>% summarise(nb_inter=length(date[piercing=="no"]),duration=unique(duration_plant),nb_cheat=length(date[piercing=="yes"]))
b$value=(b$nb_cheat+b$nb_inter)/(b$duration/24)
b$mut=(b$nb_inter)/(b$duration/24)
b$cheating=b$nb_cheat/(b$nb_cheat+b$nb_inter)
b$cheat=(b$nb_cheat)/(b$duration/24)
b=b %>% group_by(hummingbird_species) %>% mutate(cheating_moy=mean((cheat/(mut+cheat)),na.rm=T))

b2=b %>% group_by(site) %>% summarise(div=length(unique(plant_species[mut>0.5])),cheat=sum(cheat)/(sum(cheat)+sum(mut)))
b2=merge(b2,site,by="site")

pl4=ggplot(data=b2,aes(x=1,y=cheat))+geom_violin(trim=F)+geom_boxplot(width=0.2)+geom_jitter(alpha=0.5,aes(color=cheat),width=0.02)+scale_color_gradientn(colors=c("black","firebrick3"))+
theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="none",panel.border = element_blank(),axis.line= element_line(),
plot.title=element_text(size=14,face="bold",hjust = 0),
axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_text())+ylab("Frequency of cheating")+
xlab("\n")+scale_y_continuous(breaks=c(0,0.5,1),labels=c(0,0.5,1),limits = c(0,1))+
ggtitle("e")

model=glm(cheat~min_transect_elev,family="quasibinomial",data=b2)
summary(model)
pre=ggpredict(model,"min_transect_elev[600:3110]")

pl4=ggplot()+geom_point(data=b2,aes(x=min_transect_elev,y=cheat))+
theme_bw()+
geom_ribbon(data=pre,aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.2)+geom_line(data=pre,aes(x=x,y=predicted))+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="none",panel.border = element_blank(),axis.line= element_line(),
plot.title=element_text(size=14,face="bold",hjust = 0))+ylab("Overall level of cheating")+
xlab("Elevation (m)")+scale_y_continuous(breaks=c(0,0.5,1),labels=c(0,0.5,1),limits = c(0,1))+
ggtitle("e")



grid.arrange(gr1,pl1,pl2,pl3,pl4,ncol=5,widths=c(2,2,2,1,2))

setwd(dir="C:/Users/Duchenne/Documents/cheating")
pdf("figure_4.pdf",width=10,height=3)
grid.arrange(gr1,pl1,pl2,pl3,pl4,ncol=5,widths=c(2,2,2,1,2))
dev.off();























#### HYPOTHESE FOUR: CHEATING FREQUENCY IS HIGHER IN NESTED NETWORK LOW
b2$NODF=NA
for(i in b2$site){
bidon=subset(b,site==i)
bidon$tot=bidon$mut+bidon$cheat
mat=dcast(bidon,plant_species~hummingbird_species,value.var="tot",fill=0)
b2$NODF[b2$site==i]=networklevel(mat[,-1],index="weighted NODF")
}

ggplot(data=b2,aes(x=NODF,y=cheat,color=site))+geom_point()




