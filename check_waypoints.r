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
library(MuMIn)
library(INLA)
library(inlabru)
library("inlaVP")


EPHI_version="2023-07-21"
anomalies=NULL

for(pays in c("Costa-Rica","Ecuador","Brazil")){

#LOAD HUMMINGBIRD DATA FROM COSTA-RICA:
dat=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Interactions_data_",pays,".txt"),na.strings = c("",NA))
dat[dat==""]=NA
dat$date=as.IDate(dat$date,"%Y/%m/%d") #be sure date columns is recognize as date

sites=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Site_metadata_",pays,".txt"),na.strings = c("",NA))

dat$waypoint_folder=ifelse(substr(dat$folder,1,1)=="/",sapply(strsplit(dat$folder,"/"),function(x){x[3]}),sapply(strsplit(dat$folder,"/"),function(x){x[2]}))

dat=merge(dat,sites,by=c("site"),all.x=T,all.y=F)
dim(dat)

anomalies=rbind(anomalies,subset(dat,waypoint_folder!=waypoint))

}

bidon=unique(anomalies[,c("waypoint","waypoint_folder","Country")])

fwrite(unique(anomalies[,c("waypoint","waypoint_folder","Country")]),"waypoint_to_check2.csv")



EPHI_version="2023-07-21"
anomalies=NULL

for(pays in c("Costa-Rica","Ecuador","Brazil")){

#LOAD HUMMINGBIRD DATA FROM COSTA-RICA:
dat=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Interactions_data_",pays,".txt"),na.strings = c("",NA))
dat[dat==""]=NA
dat$date=as.IDate(dat$date,"%Y/%m/%d") #be sure date columns is recognize as date

sites=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Site_metadata_",pays,".txt"),na.strings = c("",NA))

dat$waypoint_folder=sapply(strsplit(dat$folder,"/"),function(x){x[2]})

dat=merge(dat,sites,by=c("site"),all.x=T,all.y=F)
dim(dat)
cameras=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Cameras_data_",pays,".txt"),na.strings = c("",NA))
cameras$end_date=as.IDate(cameras$end_date,"%Y/%m/%d") #be sure date columns is recognize as date
cameras$start_date=as.IDate(cameras$start_date,"%Y/%m/%d") #be sure date columns is recognize as date
cameras$month=month(cameras$start_date) #extract month from date column
cameras$year=year(cameras$start_date) #extract year from date column

dim(dat)
dat=merge(dat,cameras,by=c("waypoint","site"),all.x=T,all.y=T) #we want to keep cameras that did not detect any hummingbird
dim(dat)

dat=as.data.frame(dat)
timo=dat$time
timo[is.na(timo)]="12:00:00"
dat$date_time=as.POSIXct(NA,"")
dat$date_time[!is.na(dat$date)]=as.POSIXct(paste(dat$date[!is.na(dat$date)],timo[!is.na(dat$date)],sep=" "),"%Y-%m-%d %H:%M:%S",tz="CET")

timo=dat$start_time
timo[is.na(timo)]="00:00:01"
dat$start_date_time=as.POSIXct(NA,"")
dat$start_date_time[!is.na(dat$start_date)]=as.POSIXct(paste(dat$start_date[!is.na(dat$start_date)],timo[!is.na(dat$start_date)],sep=" "),"%Y-%m-%d %H:%M:%S",tz="CET")
timo=dat$end_time
timo[is.na(timo)]="23:00:01"
dat$end_date_time=as.POSIXct(NA,"")
dat$end_date_time[!is.na(dat$end_date)]=as.POSIXct(paste(dat$end_date[!is.na(dat$end_date)],timo[!is.na(dat$end_date)],sep=" "),"%Y-%m-%d %H:%M:%S",tz="CET")


anomalies=rbind(anomalies,subset(dat,date_time<start_date_time | date_time>end_date_time))

}

bidon=anomalies %>% group_by(waypoint,start_date_time,end_date_time,Country,camera_problem) %>% summarise(first_pic=min(date_time,na.rm=T),last_pic=max(date_time,na.rm=T))
bidon$difference_start=difftime(bidon$start_date_time,bidon$first_pic,units ="hours")
bidon$difference_end=difftime(bidon$last_pic,bidon$end_date_time,units ="hours")
bidon$difference_start[bidon$difference_start<=0]=NA
bidon$difference_end[bidon$difference_end<=0]=NA

bidon=subset(bidon,as.numeric(difference_start)>0 | as.numeric(difference_end)>0)


setwd(dir="C:/Users/Duchenne/Documents/EPHI_data_clean")
fwrite(bidon,"waypoint_to_check_time_final.csv")




EPHI_version="2023-06-14"
anomalies=NULL
datf=NULL

for(pays in c("Costa-Rica","Ecuador","Brazil")){

#LOAD HUMMINGBIRD DATA FROM COSTA-RICA:
dat=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Interactions_data_",pays,".txt"),na.strings = c("",NA))
dat[dat==""]=NA
dat$date=as.IDate(dat$date,"%Y/%m/%d") #be sure date columns is recognize as date

sites=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Site_metadata_",pays,".txt"),na.strings = c("",NA))
datf=rbind(datf,dat)
}



dat$ID=1:nrow(dat)
vec=dat$ID[duplicated(dat[,c("filename","time","date")]) | duplicated(dat[,c("filename","time","date")],fromLast=TRUE)]

bidon=as.data.frame(dat[dat$ID %in% vec,])
bidon=bidon[order(bidon$filename,bidon$date,bidon$time),]

fwrite(bidon,"duplicate_pictures.csv")




EPHI_version="2023-06-14"
anomalies=NULL
datf=NULL

for(pays in c("Costa-Rica","Ecuador","Brazil")){

#LOAD HUMMINGBIRD DATA FROM COSTA-RICA:
dat=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Interactions_data_",pays,".txt"),na.strings = c("",NA))
dat[dat==""]=NA
dat$date=as.IDate(dat$date,"%Y/%m/%d") #be sure date columns is recognize as date

sites=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/",pays,"_",EPHI_version,"/Site_metadata_",pays,".txt"),na.strings = c("",NA))
datf=rbind(datf,dat)
}



dat$ID=1:nrow(dat)
bidon=susbet(datf,time<"")

bidon=as.data.frame(dat[dat$ID %in% vec,])
bidon=bidon[order(bidon$filename,bidon$date,bidon$time),]

fwrite(bidon,"duplicate_pictures.csv")