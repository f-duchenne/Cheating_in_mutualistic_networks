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

EPHI_version="2023-02-14"

#COMBINE FLOWER PIERCERS AND HUMMINGBIRD DATA FROM ECUADOR:
dat_fp=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/EPHI_FP_clean_flowerpiercers.csv")
dat_fp$date=as.IDate(format(as.Date(dat_fp$date,format="%d.%m.%Y"),"%Y-%m%-%d"))
dat_ec=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/Ecuador_",EPHI_version,"/Interactions_data_Ecuador.txt"))
dat_ec[dat_ec==""]=NA
dat_ec$piercing[is.na(dat_ec$piercing)]="no"
dat_ec$date=as.IDate(dat_ec$date,"%d/%m/%Y")
dat_ec=rbind(dat_fp[,names(dat_ec),with=F],dat_ec)
dat_ec$year=year(dat_ec$date)
#MERGE THEM WITH CAMERA INFORMATION:
cameras=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/Ecuador_",EPHI_version,"/Cameras_data_Ecuador.txt"),na.strings = c("",NA))
cameras$month=month(cameras$start_date)
cameras$year=year(cameras$start_date)
cameras$plant_species[cameras$plant_pecies==""]=NA
dim(dat_ec)
dat_ec=merge(subset(dat_ec,!is.na(waypoint) & !is.na(site)),cameras,by=c("waypoint","site"),all.x=T,all.y=T) #if you want to exclude camera which did not detect any hummingbird
dim(dat_ec)
dat_ec$Country="Ecuador"
dim(data)
dat_ec=dat_ec[dat_ec$site %in% unique(dat_fp$site),]
dim(data)

#LOAD HUMMINGBIRD DATA FROM COSTA-RICA:
dat_cr=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_",EPHI_version,"/Interactions_data_Costa-Rica.txt"))
dat_cr$date=as.IDate(dat_cr$date,"%Y-%m%-%d")
dat_cr$year=year(dat_cr$date)
dat_cr[dat_cr==""]=NA
#MERGE THEM WITH CAMERA INFORMATION:
cameras=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_",EPHI_version,"/Cameras_data_Costa-Rica.txt"),na.strings = c("",NA))
cameras$end_date=as.IDate(cameras$end_date,"%Y-%m%-%d")
cameras$start_date=as.IDate(cameras$start_date,"%Y-%m%-%d")
cameras$month=month(cameras$start_date)
cameras$year=year(cameras$start_date)
cameras$plant_species[cameras$plant_pecies==""]=NA
dim(dat_cr)
dat_cr=merge(subset(dat_cr,!is.na(waypoint) & !is.na(site)),cameras,by=c("waypoint","site"),all.x=T,all.y=T) #if you want to keep camera which did not detect any hummingbird
dim(dat_cr)
dat_cr$Country="Costa-Rica"

#COMBINE BOTH DATASETS
dat=rbind(dat_cr,dat_ec[,names(dat_cr),with=F])
dim(subset(dat,is.na(duration_sampling_hours)))
###### ESTIMATE MISSING DURATION FROM PICTURE TIMES:
dat$time[nchar(dat$time)==5 & !is.na(dat$time)]=paste0(dat$time[nchar(dat$time)==5 & !is.na(dat$time)],":00")
guessduration <- function(dt1,tm1,dt2,tm2) {
  if(!is.na(tm1)) {
    obj1=strptime(paste(dt1, tm1), "%Y-%m-%d %H:%M:%S")
  }else{
    obj1=strptime(paste(dt1, "10:00:00"), "%Y-%m-%d %H:%M:%S")
  }
  if(!is.na(tm2)) {
    obj2=strptime(paste(dt2, tm2), "%Y-%m-%d %H:%M:%S")
  }else{
    obj2=strptime(paste(dt2, "10:00:00"), "%Y-%m-%d %H:%M:%S")
  }
  return(difftime(obj2,obj1,units="hours"))
}


dat[,date_complete_pic:=strptime(paste(dat$date, dat$time,sep=" "), "%Y-%m-%d %H:%M:%S")]
dat=dat %>% group_by(waypoint,site) %>% mutate(start_date_pic=min(date),end_date_pic=max(date),
start_time_pic=format(min(date_complete_pic), format="%H:%M:%S"),end_time_pic=format(max(date_complete_pic), format="%H:%M:%S"))

dat$duration_sampling_hours[is.na(dat$duration_sampling_hours)]=mapply(guessduration,dat$start_date_pic[is.na(dat$duration_sampling_hours)],dat$start_time_pic[is.na(dat$duration_sampling_hours)],
dat$end_date_pic[is.na(dat$duration_sampling_hours)],dat$end_time_pic[is.na(dat$duration_sampling_hours)])
dim(subset(dat,is.na(duration_sampling_hours)))

#REMOVE CAMERA WITHOUT DURATION
dat=subset(dat,duration_sampling_hours>0)

##### EXPORTING DATA
setwd(dir="C:/Users/Duchenne/Documents/cheating")
fwrite(dat,"empirical_data.csv")

nrow(subset(dat,Country=="Ecuador" & !is.na(hummingbird_species))) #n interactions Ecuador
summary(subset(dat,Country=="Ecuador")$date) #dates
nrow(subset(dat,Country=="Costa-Rica" & !is.na(hummingbird_species))) #n interactions Costa-Rica
summary(subset(dat,Country=="Costa-Rica")$date) #dates


getmode <- function(v) {
 uniqv <- unique(na.omit(v))
 uniqv[which.max(tabulate(match(v, uniqv)))]
}

#LOAD TRAIT DATA AND COMBINE THEM TO HAVE ONE VALUE PER SPECIES
tr=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/plant_traits_",EPHI_version,"/Plant_traits.txt"),na.strings = c("",NA))
tr2=subset(tr,!is.na(plant_species))  %>% group_by(plant_species,plant_family,plant_genus) %>% summarise(Tube_length=mean(Tube_length,na.rm=T),Anther_length=mean(Anther_length,na.rm=T), Stigma_length=mean(Stigma_length,na.rm=T),
Opening_corrolla=mean(Opening_corrolla,na.rm=T),Curvature_middle=mean(Curvature_middle,na.rm=T),Type=getmode(Type))

setwd(dir="C:/Users/Duchenne/Documents/cheating")
fwrite(tr2,"traits_data.csv")


#LOAD SITE METADATA
site_cr=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_",EPHI_version,"/Site_metadata_Costa-Rica.txt"))
site_ec=fread(paste0("C:/Users/Duchenne/Documents/EPHI_data_clean/Ecuador_",EPHI_version,"/Site_metadata_Ecuador.txt"))
sites=rbind(site_cr,site_ec)

setwd(dir="C:/Users/Duchenne/Documents/cheating")
fwrite(as.data.frame(sites)[sites$site %in% unique(dat$site),],"table_s2.csv")