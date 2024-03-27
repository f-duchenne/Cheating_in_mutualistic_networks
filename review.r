library(ggplot2)
tab1=data.frame(x=-1/2*rbeta(1000000,1,0.7),b=0.7)
tab2=data.frame(x=-1/2*rbeta(1000000,1,1),b=1)
tab3=data.frame(x=-1/2*rbeta(1000000,1,1.5),b=1.5)
tab4=data.frame(x=-1/2*rbeta(1000000,1,3),b=3)
tabf=rbind(tab1,tab2,tab3,tab4)

ggplot(data=tabf,aes(x=x,col=as.factor(b)))+geom_density(fill=NA,alpha=0,adjust =2)+
labs(color="shape parameter (b)")+theme_bw()+theme(panel.grid=element_blank())+xlab("Growth rate value")



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
library(sp)
library(spatialEco)
setwd(dir="C:/Users/Duchenne/Documents/cheating")

EPHI_version="2023-08-24"

#MERGE THEM WITH CAMERA INFORMATION:
cameras_ec=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/Ecuador_",EPHI_version,"/Cameras_data_Ecuador.txt"),na.strings = c("",NA))
cameras_ec$month=month(cameras_ec$start_date)
cameras_ec$year=year(cameras_ec$start_date)
cameras_ec$plant_species[cameras_ec$plant_pecies==""]=NA

#MERGE THEM WITH CAMERA INFORMATION:
cameras_cr=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_",EPHI_version,"/Cameras_data_Costa-Rica.txt"),na.strings = c("",NA))
cameras_cr$end_date=as.IDate(cameras_cr$end_date,"%Y-%m%-%d")
cameras_cr$start_date=as.IDate(cameras_cr$start_date,"%Y-%m%-%d")
cameras_cr$month=month(cameras_cr$start_date)
cameras_cr$year=year(cameras_cr$start_date)
cameras_cr$plant_species[cameras_cr$plant_pecies==""]=NA

#COMBINE BOTH DATASETS
cameras=rbind(cameras_ec,cameras_cr)
cameras_agg=cameras %>% group_by(site,plant_species) %>% summarise(n=length(unique(waypoint)))


transect_cr=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_",EPHI_version,"/Transect_data_Costa-Rica.txt"),na.strings = c("",NA))
transect_ec=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/Ecuador_",EPHI_version,"/Transect_data_Ecuador.txt"),na.strings = c("",NA))
transect_ec=transect_ec[transect_ec$site %in% c("Alaspungo","Alaspungo_disturbed","Yanacocha_disturbed","Verdecocha","Yanacocha"),]
#COMBINE BOTH DATASETS
transect=rbind(transect_ec,transect_cr)

transect_agg=transect %>% group_by(site,plant_species) %>% count()


tab=merge(transect_agg,cameras_agg,by=c("site","plant_species"),all.x=T)
tab=subset(tab,!is.na(plant_species))
tab$n.y[is.na(tab$n.y)]=0

tr2=fread("traits_data.csv",na.strings = c(""," ", NA))
tab=merge(tab,tr2,by="plant_species",all.x=T)
tab$n.y[is.na(tab$n.y)]=0
tab %>% group_by(site) %>% summarise(perc=length(unique(plant_species[n.y==0]))/length(unique(plant_species)))


tab$tubular=NA
tab$tubular[tab$Tube_length>0.1]="yes"
tab$tubular[tab$Tube_length<=0.1]="no"

#FIG. S8
ggplot(data=tab,aes(x=n.x,y=n.y,color=tubular))+geom_point()+facet_wrap(~site,ncol=3)+
scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),
n.breaks =2,minor_breaks=NULL)+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),
n.breaks =2,minor_breaks=NULL)+
theme_bw()+ annotation_logticks()+theme(panel.grid=element_blank())+xlab("Abundance on site")+
ylab("Sampling pressure (number of sampling events)")


fwrite(tab,"data_for_figure_S9.csv")
