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

#COMBINE FLOWER PIERCERS AND HUMMINGBIRD DATA FROM ECUADOR:
dat_fp=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/EPHI_FP_clean_flowerpiercers.csv")
dat_fp$date=as.IDate(format(as.Date(dat_fp$date,format="%d.%m.%Y"),"%Y-%m%-%d"))
dat_ec=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/Ecuador_2022-11-07/Interactions_data_Ecuador.txt")
dat_ec[dat_ec==""]=NA
dat_ec$piercing[is.na(dat_ec$piercing)]="no"
dat_ec$date=as.IDate(dat_ec$date,"%d/%m/%Y")
dat_ec=rbind(dat_fp[,names(dat_ec),with=F],dat_ec)
dat_ec$year=year(dat_ec$date)

#MERGE THEM WITH CAMERA INFORMATION:
cameras=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/Ecuador_2022-11-07/Cameras_data_Ecuador.txt",na.strings = c("",NA))
cameras$month=month(cameras$start_date)
cameras$year=year(cameras$start_date)
dim(dat_ec)
dat_ec=merge(subset(dat_ec,!is.na(waypoint) & !is.na(site) & !is.na(year)),cameras,by=c("waypoint","site"),all.x=T,all.y=F) #if you want to exclude camera which did not detect any hummingbird
dim(dat_ec)
dat_ec$Country="Ecuador"
dim(data)
dat_ec=dat_ec[dat_ec$site %in% unique(dat_fp$site),]
dim(data)

#LOAD HUMMINGBIRD DATA FROM COSTA-RICA:
dat_cr=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_2022-11-07/Interactions_data_Costa-Rica.txt")
dat_cr$date=as.IDate(dat_cr$date,"%Y-%m%-%d")
dat_cr$year=year(dat_cr$date)
#MERGE THEM WITH CAMERA INFORMATION:
cameras=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_2022-11-07/Cameras_data_Costa-Rica.txt",na.strings = c("",NA))
cameras$end_date=as.IDate(cameras$end_date,"%Y-%m%-%d")
cameras$start_date=as.IDate(cameras$start_date,"%Y-%m%-%d")
cameras$month=month(cameras$start_date)
cameras$year=year(cameras$start_date)
dim(dat_cr)
dat_cr=merge(dat_cr,cameras,by=c("waypoint","site"),all.x=T,all.y=F) #if you want to exclude camera which did not detect any hummingbird
dim(dat_cr)
dat_cr$Country="Costa-Rica"

#COMBINE BOTH DATASETS
dat=rbind(dat_cr,dat_ec[,names(dat_cr),with=F])
dat=subset(dat,duration_sampling_hours>0)

getmode <- function(v) {
 uniqv <- unique(na.omit(v))
 uniqv[which.max(tabulate(match(v, uniqv)))]
}

#LOAD TRAIT DATA AND COMBINE THEM TO HAVE ONE VALUE PER SPECIES
tr=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/plant_traits_2022-11-07/Plant_traits.txt",na.strings = c("",NA))
tr2=tr %>% group_by(plant_species,Type,SexualSystem,Symmetry) %>% summarise(Tube_length=mean(Tube_length,na.rm=T),Anther_length=mean(Anther_length,na.rm=T), Stigma_length=mean(Stigma_length,na.rm=T),
Opening_corrolla=mean(Opening_corrolla,na.rm=T),Curvature_middle=mean(Curvature_middle,na.rm=T))
tr2[tr2$plant_species %in% tr2$plant_species[duplicated(tr2$plant_species)],]
tr2=subset(tr,!is.na(plant_species))  %>% group_by(plant_species,plant_family,plant_genus) %>% summarise(Tube_length=mean(Tube_length,na.rm=T),Anther_length=mean(Anther_length,na.rm=T), Stigma_length=mean(Stigma_length,na.rm=T),
Opening_corrolla=mean(Opening_corrolla,na.rm=T),Curvature_middle=mean(Curvature_middle,na.rm=T),Type=getmode(Type),sexualsys=getmode(SexualSystem))


#LOAD SITE METADATA
site_cr=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/Costa-Rica_2022-11-07/Site_metadata_Costa-Rica.txt")
site_ec=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/Ecuador_2022-11-07/Site_metadata_Ecuador.txt")
sites=rbind(site_cr,site_ec)

fwrite(as.data.frame(sites)[sites$site %in% unique(dat$site),],"table_s2.csv")

###### PLOT NETWORKS:

bi=dat %>% group_by(plant_species,site,Country) %>% mutate(duration_plant=sum(duration_sampling_hours[!duplicated(waypoint)],na.rm=T))
b=bi %>% group_by(plant_species,hummingbird_species,site,Country) %>% summarise(nb_inter=length(date[piercing=="no"]),duration=unique(duration_plant),nb_cheat=length(date[piercing=="yes"]))
b$value=(b$nb_cheat+b$nb_inter)/(b$duration/24)
b$cheating=b$nb_cheat/(b$nb_cheat+b$nb_inter)
b=subset(b,!is.na(value) & !is.na(plant_species) & !is.na(hummingbird_species))
rbPal <- colorRampPalette(c("black","pink",'red'))
b$Col <- rbPal(20)[as.numeric(cut(b$cheating,breaks = 20))]
png("figS7.png",width=1600,height=1400,res=120)
par(mfrow=c(3,4))
for(i in unique(dat$site[dat$Country=="Costa-Rica"])){
bidon=subset(b,site==i & Country=="Costa-Rica")
mat=dcast(bidon,plant_species~hummingbird_species,value.var="value",fill=0)

plotweb(sqrt(mat[,-1]),col.interaction=subset(b,site==i)$Col,col.high="dodgerblue3",col.low="forestgreen",bor.col.interaction=NA,high.lablength=0,low.lablength=0,bor.col.high="dodgerblue3",bor.col.low="forestgreen",y.width.low=0.05,y.width.high=0.05)
title(paste0(i))
}
dev.off();


png("figS8.png",width=1000,height=1400,res=120)
par(mfrow=c(3,2))
for(i in unique(dat$site[dat$Country=="Ecuador"])){
bidon=subset(b,site==i & Country=="Ecuador")
mat=dcast(bidon,plant_species~hummingbird_species,value.var="value",fill=0)
plotweb(sqrt(mat[,-1]),col.interaction=bidon$Col,col.high="dodgerblue3",col.low="forestgreen",bor.col.interaction=NA,high.lablength=0,low.lablength=0,bor.col.high="dodgerblue3",bor.col.low="forestgreen",y.width.low=0.05,y.width.high=0.05)
title(paste0(i))
}
dev.off();


# An example of network:
site_choisi="BOQU"
bidon=subset(b,site==site_choisi)
mat=dcast(bidon,plant_species~hummingbird_species,value.var="value",fill=0)
mat2=dcast(bidon,plant_species~hummingbird_species,value.var="cheating",fill=0)

rbPal <- colorRampPalette(c("gold3",'red'))
vec=c()
for(i in 1:nrow(mat)){
vec=c(vec,as.numeric(t(mat2)[-1,i]))
}

color=rbPal(20)[as.numeric(cut(vec,breaks = 20))]

plotweb(as.matrix(sqrt(mat[,-1])),method="normal",col.interaction=color)

plot_grid(base2grob(~plotweb(as.matrix(sqrt(mat[,-1])),method="normal",col.interaction=color,col.high="dodgerblue3",col.low="forestgreen",bor.col.interaction=NA,high.lablength=0,low.lablength=0,bor.col.high="dodgerblue3",bor.col.low="forestgreen",
y.width.low=0.05,y.width.high=0.05)))+
ggtitle("a")+theme(plot.title=element_text(size=14,face="bold",hjust = 0))

#### HYPOTHESE ONE: CHEATERS ARE SPECIALISTS
dat=dat %>% group_by(plant_species,site,Country) %>% mutate(duration_plant=sum(duration_sampling_hours[!duplicated(waypoint)],na.rm=T))
b=subset(dat,!is.na(duration_plant) & !is.na(plant_species) & !is.na(hummingbird_species)) %>% group_by(plant_species,hummingbird_species,site,Country) %>% summarise(nb_inter=length(date[piercing=="no"]),duration=unique(duration_plant),nb_cheat=length(date[piercing=="yes"]))
b$value=(b$nb_cheat+b$nb_inter)/(b$duration/24)
b$mut=(b$nb_inter)/(b$duration/24)
b$cheating=b$nb_cheat/(b$nb_cheat+b$nb_inter)
b$cheat=(b$nb_cheat)/(b$duration/24)
b=b %>% group_by(hummingbird_species,Country) %>% mutate(cheating_moy=mean((cheat/(mut+cheat)),na.rm=T),div=length(unique(plant_species[mut>0.1])),total_inter=sum(cheat+mut))

b1=b %>% group_by(hummingbird_species,site,Country) %>% summarise(div=length(unique(plant_species[mut>0])),cheat=sum(cheat)/(sum(cheat)+sum(mut)),div2=vegan::diversity(mut),total_inter=sum(cheat)+sum(mut))

b1$divlog=log(b1$div+1)
model=glmmTMB(cheat~divlog*Country+(1|site)+(1|hummingbird_species),family="binomial",data=b1)
summary(model)
ano1=Anova(model)
pre=ggpredict(model,c("divlog[all]","Country"))
pre[pre$x>max(b1$divlog[b1$Country=="Costa-Rica"]) & pre$group=="Costa-Rica",c("predicted","conf.low","conf.high")]=NA
pre$group=factor(pre$group,levels=c("Ecuador","Costa-Rica"))
b1$Country=factor(b1$Country,levels=c("Ecuador","Costa-Rica"))
pl2=ggplot()+geom_point(data=b1,aes(x=divlog,y=cheat,size=total_inter,color=Country))+geom_ribbon(data=pre,aes(x=x,ymin=conf.low,ymax=conf.high,fill=group),alpha=0.1)+geom_line(data=pre,aes(x=x,y=predicted,color=group))+
theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",panel.border = element_blank(),axis.line= element_line())+ggtitle("b")+
xlab("Number of partners")+ylab("Frequency of cheating")+scale_y_continuous(breaks=c(0,0.5,1),labels=c(0,0.5,1),limits = c(0,1))+scale_x_continuous(breaks=log(c(0,5,10,20,35)+1),labels=c(0,5,10,20,35))+
scale_color_manual(values=c("#0B4F6C","#CBB9A8"))+scale_fill_manual(values=c("#0B4F6C","#CBB9A8"))+scale_size(range = c(0.5,2))

##### HYPOTHESE TWO: CHEATING IS INNOVATIVE
b2=b %>% group_by(hummingbird_species,plant_species,cheating_moy,div,total_inter,Country,site) %>% summarise(dens_mut=sum(mut),dens_cheat=sum(cheat))
b2=b2 %>% group_by(hummingbird_species,cheating_moy,div,total_inter,Country,site) %>% mutate(some_mut=sum(dens_mut),some_cheat=sum(dens_cheat))

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

tr2$Curvature_middle_pre[!is.na(tr2$Curvature_middle_pre) & is.na(tr2$Curvature_middle)]

b2=merge(b2,tr2,by=c("plant_species"),all.x=T,all.y=F) #if you want to exclude camera which did not detect any hummingbird
length(b2$Tube_length_pre[!is.na(b2$Tube_length_pre) & is.na(b2$Tube_length)])
length(b2$Curvature_middle_pre[!is.na(b2$Curvature_middle_pre) & is.na(b2$Curvature_middle)])
length(b2$Curvature_middle_pre[is.na(b2$Curvature_middle_pre) | is.na(b2$Tube_length_pre)])

b2$innovative=NA

list_cheaters=unique(subset(b2,some_cheat>0 & !is.na(Tube_length) & !is.na(Curvature_middle_pre) & !is.na(Tube_length))[,c("hummingbird_species","Country")])
#TWO DIMENSIONS
resf=NULL
list_plot=list()
for(i in 1:nrow(list_cheaters)){
bidon=subset(b2,hummingbird_species==list_cheaters$hummingbird_species[i] & !is.na(Tube_length) & Country==list_cheaters$Country[i] & !is.na(Curvature_middle_pre))
if(nrow(bidon)>2){
bidon_sanscheat=subset(bidon,dens_mut>0)
poly1 <- sp::Polygon(bidon_sanscheat[chull(bidon_sanscheat[,c("Tube_length_pre","Curvature_middle_pre")]),c("Tube_length_pre","Curvature_middle_pre")])
Poly <- sp::Polygons(list(poly1), ID = "A")
SpatialPoly <- sp::SpatialPolygons(list(Poly))

point=SpatialPointsDataFrame(as.matrix(bidon[,c("Tube_length_pre","Curvature_middle_pre")]),data.frame(plant_species=bidon[,"plant_species"]))

points_poly=as.data.frame(point.in.poly(point,SpatialPoly))
points_poly$poly.ids[is.na(points_poly$poly.ids)]=0
bidon$innovative=1-points_poly$poly.ids
for(j in bidon$plant_species){b2$innovative[b2$hummingbird_species==list_cheaters$hummingbird_species[i] & b2$Country==list_cheaters$Country[i] &b2$plant_species==j]=bidon$innovative[bidon$plant_species==j]}


weight1=bidon$dens_mut/bidon$some_mut
weight1[is.na(weight1)]=0
weight2=bidon$dens_cheat/bidon$some_cheat
weight2[is.na(weight2)]=0

if(sum(weight2)>0){
#calculate bandwidth to use
mati=Hpi(bidon[,c("Tube_length_pre","Curvature_middle_pre")])
### Kernel functional estimate for 1- to 6-dimensional data:
obj=kde(bidon[,c("Tube_length_pre","Curvature_middle_pre")],
w=weight1,gridsize=50,H=mati,xmin=c(0,0),xmax=c(max(b2$Tube_length_pre,na.rm=T),max(b2$Curvature_middle_pre,na.rm=T)))
image2D(obj$estimate)
resi=melt(obj$estimate)
names(resi)[1:3]=c("Tube_length","Curvature_middle","dens_mut")
resi$Tube_length=obj$eval.points[[1]]
resi$Curvature_middle=rep(obj$eval.points[[2]],each=length(obj$eval.points[[1]]))

obj=kde(bidon[,c("Tube_length_pre","Curvature_middle_pre")],
w=weight2,gridsize=50,H=mati,xmin=c(0,0),xmax=c(max(b2$Tube_length_pre,na.rm=T),max(b2$Curvature_middle_pre,na.rm=T)))
resi$dens_cheat=melt(obj$estimate)$value
resi$hummingbird_species=list_cheaters$hummingbird_species[i]
resi$Country=list_cheaters$Country[i]
resi$cheating_moy=unique(bidon$cheating_moy)
resi$div=unique(bidon$div)
resi$total_inter=unique(bidon$total_inter)


list_plot[[(length(list_plot)+1)]]=ggplot()+geom_tile(data=resi,aes(x=Tube_length,y=Curvature_middle,fill=dens_mut))+scale_fill_gradientn(colours = terrain.colors(10))+
geom_point(data=bidon,aes(x=Tube_length_pre,y=as.numeric(Curvature_middle_pre),size=dens_cheat,alpha=0.6))+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="none",
plot.title=element_text(size=8,hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=12))+xlab("")+
ylab("")+coord_cartesian(expand = FALSE)+ggtitle(paste(list_cheaters$hummingbird_species[i],list_cheaters$Country[i],sep =" - "))

#check:
sum(resi$dens_mut)*(resi$Tube_length[2]-resi$Tube_length[1])*(resi$Curvature_middle[2]-resi$Curvature_middle[1])

resf=rbind(resf,resi)
}
}
}

list_plot = list_plot[lapply(list_plot,length)>0]
grid.arrange(grobs = list_plot , ncol = 5,bottom ="Corolla length (cm)",left="Curvature measure")

resf$mini=apply(resf[,c("dens_mut","dens_cheat")],1,min)
b2f=resf %>% dplyr::group_by(hummingbird_species,cheating_moy,div,total_inter,Country) %>% dplyr::summarise(prop_innovative=1-sum(mini)*(Tube_length[2]-Tube_length[1])*(obj$eval.points[[2]][2]-obj$eval.points[[2]][1]))

pl3=ggplot(data=b2f,aes(x=Country,y=prop_innovative,color=Country))+geom_boxplot(width=0.3)+geom_jitter(alpha=0.5,aes(size=total_inter),width=0.2)+
theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="none",panel.border = element_blank(),,axis.line= element_line(),
plot.title=element_text(size=14,face="bold",hjust = 0),
axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_text())+ylab("Proportion of innovative cheating")+
xlab("\n")+scale_y_continuous(breaks=c(0,0.5,1),labels=c(0,0.5,1),limits = c(0,1))+
ggtitle("c")+scale_color_manual(values=c("#0B4F6C","#CBB9A8"))+scale_fill_manual(values=c("#0B4F6C","#CBB9A8"))+scale_size(range = c(0.5,2))

###### HYPOTHESE THREE: LEVEL OF CHEATING IS LOW AND DECREASE WITH COST
b3=b1 %>% group_by(site) %>% summarise(prop_cheaters=length(hummingbird_species[cheat>0.05])/length(hummingbird_species))
b3=merge(b3,sites,by="site")

model=glm(prop_cheaters~min_transect_elev*Country,family="quasibinomial",data=b3)
summary(model)
ano2=Anova(model)
pre=ggpredict(model,c("min_transect_elev[605:3406]","Country"))
pre[pre$x>max(b3$min_transect_elev[b3$Country=="Costa-Rica"]) & pre$group=="Costa-Rica",c("predicted","conf.low","conf.high")]=NA
pre[pre$x<min(b3$min_transect_elev[b3$Country=="Ecuador"]) & pre$group=="Ecuador",c("predicted","conf.low","conf.high")]=NA
b3$Country=factor(b3$Country,levels=c("Ecuador","Costa-Rica"))
pre$group=factor(pre$group,levels=c("Ecuador","Costa-Rica"))
pl1=ggplot()+geom_point(data=b3,aes(x=min_transect_elev,y=prop_cheaters,color=Country))+
theme_bw()+
geom_ribbon(data=pre,aes(x=x,ymin=conf.low,ymax=conf.high,fill=group),alpha=0.2)+geom_line(data=pre,aes(x=x,y=predicted,color=group))+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="none",panel.border = element_blank(),axis.line= element_line(),
plot.title=element_text(size=14,face="bold",hjust = 0))+ylab("Proportion of cheaters")+
xlab("Elevation (m)")+scale_y_continuous(breaks=c(0,0.5,1),labels=c(0,0.5,1),limits = c(0,1))+
ggtitle("d")+scale_color_manual(values=c("#0B4F6C","#CBB9A8"))+scale_fill_manual(values=c("#0B4F6C","#CBB9A8"))

b4=b %>% group_by(site,Country) %>% summarise(cheat=sum(cheat)/(sum(cheat)+sum(mut)),nbsp_tot=length(unique(hummingbird_species))+length(unique(plant_species)))
b4=merge(b4,sites,by=c("site","Country"))

model=glm(cheat~min_transect_elev*Country,family="quasibinomial",data=b4)
summary(model)
ano3=Anova(model)
b4$resi=residuals(model)
fwrite(b4[,c("site","resi","cheat")],"resi.txt")
pre=ggpredict(model,c("min_transect_elev[605:3406]","Country"))
pre[pre$x>max(b4$min_transect_elev[b4$Country=="Costa-Rica"]) & pre$group=="Costa-Rica",c("predicted","conf.low","conf.high")]=NA
pre[pre$x<min(b4$min_transect_elev[b4$Country=="Ecuador"]) & pre$group=="Ecuador",c("predicted","conf.low","conf.high")]=NA
pre$group=factor(pre$group,levels=c("Ecuador","Costa-Rica"))
b4$Country=factor(b4$Country,levels=c("Ecuador","Costa-Rica"))
pl4=ggplot()+geom_point(data=b4,aes(x=min_transect_elev,y=cheat,color=Country))+
theme_bw()+
geom_ribbon(data=pre,aes(x=x,ymin=conf.low,ymax=conf.high,fill=group),alpha=0.2)+geom_line(data=pre,aes(x=x,y=predicted,color=group))+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="none",panel.border = element_blank(),axis.line= element_line(),
plot.title=element_text(size=14,face="bold",hjust = 0))+ylab("Overall level of cheating")+
xlab("Elevation (m)")+scale_y_continuous(breaks=c(0,0.5,1),labels=c(0,0.5,1),limits = c(0,1))+
ggtitle("e")+scale_color_manual(values=c("#0B4F6C","#CBB9A8"))+scale_fill_manual(values=c("#0B4F6C","#CBB9A8"))+scale_size(range = c(0.5,2))


##### OVERALL FIGURE AND TABLES
grid.arrange(gr1,pl2,pl3,pl1,pl4,ncol=5,widths=c(2,2,1,2,2))

setwd(dir="C:/Users/Duchenne/Documents/cheating")
pdf("figure_4.pdf",width=10,height=3)
grid.arrange(gr1,pl2,pl3,pl1,pl4,ncol=5,widths=c(2,2,1,2,2))
dev.off();


names(ano1)=names(ano2)
ano1$varia=rownames(ano1)
ano2$varia=rownames(ano2)
ano3$varia=rownames(ano3)
fwrite(rbind(ano1,ano2,ano3),"anova.csv")


#######PREPARE MATRIX FOR SIMULATIONS
for (i in unique(sites$site[sites$site %in% b$site])){
mut=subset(b2,site==i & dens_mut>0)
mut_mat=dcast(mut,plant_species~hummingbird_species,value.var="dens_mut",fill=0)
list_plant=mut_mat[,1]
list_hum=names(mut_mat[,-1])
mut_mat=as.matrix(mut_mat[,-1])

cheat=subset(b2,site==i & dens_mut>0)
cheat_mat=dcast(cheat,plant_species~hummingbird_species,value.var="dens_cheat",fill=0)
cheat_mat=as.matrix(cheat_mat[,-1])

for(j in 1:100){
#growth rates:
shape_a=runif(1,log(0.7),log(3))
shape_p=runif(1,log(0.7),log(3))
#r=c(runif(nbsp_h,-1,0),runif(nbsp_p,0,0.5))
r=c(-1*rbeta(length(list_hum),1,exp(shape_a))/2,-1*rbeta(length(list_plant),1,exp(shape_p))/2)
r[r>(-0.001)]=-0.001

bibi=tr2[tr2$plant_species %in% list_plant,]
bibi=bibi[match(list_plant, bibi$plant_species),]
r2=r
r[length(list_hum)+which(bibi$sexualsys %in% c("Homogamy"))]=-1*r[length(list_hum)+which(bibi$sexualsys %in% c("Homogamy"))]
r[grep("Diglossa",list_hum,fixed=T)]=-1*r[grep("Diglossa",list_hum,fixed=T)]

setwd(dir="C:/Users/Duchenne/Documents/cheating/initial_empir")
save(mut_mat,cheat_mat,list_plant,list_hum,r,r2,file=paste0(i,"_",j,".RData"),version = 2)
}
}












bi=dat %>% group_by(plant_species,site,Country) %>% mutate(duration_plant=sum(duration_sampling_hours[!duplicated(waypoint)],na.rm=T))
b3=merge(bi,tr2,by=c("plant_species"),all.x=T,all.y=F) #if you want to exclude camera which did not detect any hummingbird
b3=merge(b3,site,by="site")



theory=fread("theory_results_for_combine.txt")


#####INNOVATIVE PER SITE
b=subset(dat,!is.na(duration_plant) & !is.na(plant_species) & !is.na(hummingbird_species)) %>% group_by(plant_species,hummingbird_species,site,Country) %>% summarise(nb_inter=length(date[piercing=="no"]),duration=unique(duration_plant),nb_cheat=length(date[piercing=="yes"]))
b$value=(b$nb_cheat+b$nb_inter)/(b$duration/24)
b$mut=(b$nb_inter)/(b$duration/24)
b$cheating=b$nb_cheat/(b$nb_cheat+b$nb_inter)
b$cheat=(b$nb_cheat)/(b$duration/24)
b=b %>% group_by(hummingbird_species,Country) %>% mutate(cheating_moy=mean((cheat/(mut+cheat)),na.rm=T),div=length(unique(plant_species[mut>0.1])),total_inter=sum(cheat+mut))

b3=b %>% group_by(hummingbird_species,plant_species,cheating_moy,div,total_inter,Country,site) %>% summarise(dens_mut=sum(mut),dens_cheat=sum(cheat))
b3=b3 %>% group_by(hummingbird_species,cheating_moy,div,total_inter,Country,site) %>% mutate(some_mut=sum(dens_mut),some_cheat=sum(dens_cheat))

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

tr2$Curvature_middle_pre[!is.na(tr2$Curvature_middle_pre) & is.na(tr2$Curvature_middle)]

b3=merge(b3,tr2,by=c("plant_species"),all.x=T,all.y=F) #if you want to exclude camera which did not detect any hummingbird
length(b3$Tube_length_pre[!is.na(b3$Tube_length_pre) & is.na(b3$Tube_length)])
length(b3$Curvature_middle_pre[!is.na(b3$Curvature_middle_pre) & is.na(b3$Curvature_middle)])
length(b3$Curvature_middle_pre[is.na(b3$Curvature_middle_pre) | is.na(b3$Tube_length_pre)])


list_cheaters=unique(subset(b3,some_cheat>0 & !is.na(Tube_length) & !is.na(Curvature_middle_pre) & !is.na(Tube_length))[,c("hummingbird_species","Country","site")])

#TWO DIMENSIONS
resf=NULL
list_plot=list()
for(i in 1:nrow(list_cheaters)){
bidon=subset(b3,hummingbird_species==list_cheaters$hummingbird_species[i] & !is.na(Tube_length) & !is.na(Curvature_middle_pre))
bidon2=subset(b3,hummingbird_species==list_cheaters$hummingbird_species[i] & !is.na(Tube_length) & Country==list_cheaters$Country[i] & !is.na(Curvature_middle_pre) & site==list_cheaters$site[i])
if(nrow(bidon)>2){
weight1=bidon$dens_mut/bidon$some_mut
weight1[is.na(weight1)]=0
weight2=bidon2$dens_cheat/bidon2$some_cheat
weight2[is.na(weight2)]=0

if(sum(weight2)>0){
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

obj=kde(bidon2[,c("Tube_length_pre","Curvature_middle_pre")],
w=weight2,gridsize=50,H=mati,xmin=c(0,0),xmax=c(max(b3$Tube_length_pre,na.rm=T),max(b3$Curvature_middle_pre,na.rm=T)))
resi$dens_cheat=melt(obj$estimate)$value
resi$hummingbird_species=list_cheaters$hummingbird_species[i]
resi$Country=list_cheaters$Country[i]
resi$site=list_cheaters$site[i]
resi$cheating_moy=unique(bidon$cheating_moy)
resi$div=unique(bidon$div)
resi$total_inter=unique(bidon$total_inter)


list_plot[[(length(list_plot)+1)]]=ggplot()+geom_tile(data=resi,aes(x=Tube_length,y=Curvature_middle,fill=dens_mut))+scale_fill_gradientn(colours = terrain.colors(10))+
geom_point(data=bidon,aes(x=Tube_length_pre,y=as.numeric(Curvature_middle_pre),size=dens_cheat,alpha=0.6))+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="none",
plot.title=element_text(size=8,hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=12))+xlab("")+
ylab("")+coord_cartesian(expand = FALSE)+ggtitle(paste(list_cheaters$hummingbird_species[i],list_cheaters$Country[i],sep =" - "))

#check:
sum(resi$dens_mut)*(resi$Tube_length[2]-resi$Tube_length[1])*(resi$Curvature_middle[2]-resi$Curvature_middle[1])

resf=rbind(resf,resi)
}
}
}
resf$mini=apply(resf[,c("dens_mut","dens_cheat")],1,min)
b4=resf %>% dplyr::group_by(hummingbird_species,site) %>% dplyr::summarise(prop_innovative=1-sum(mini)*(Tube_length[2]-Tube_length[1])*(obj$eval.points[[2]][2]-obj$eval.points[[2]][1]))
b4f=b4 %>% dplyr::group_by(site) %>% dplyr::summarise(prop_innovative=mean(prop_innovative))


#CHEATING FREQUENCY & PROPORTION OF CHEATERS
b2=b %>% group_by(site) %>% summarise(cheat=sum(cheat)/(sum(cheat)+sum(mut)),prop_cheaters=length(hummingbird_species[cheating>0.05])/length(hummingbird_species))
b2=merge(b2,site,by="site")

point=merge(b2,b4f,by="site")


ggplot()+
geom_raster(data=theory,aes(x=prop_cheating,y=prop_cheaters,fill=opti,alpha=exp(pers_norm)^10))+theme_bw()+facet_grid(rows=vars(cost),cols=vars(scenario),labeller = label_bquote(rows=Lambda == .(cost)))+
scale_x_continuous(breaks=seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1))+
scale_y_continuous(breaks=seq(0.1,0.9,0.2),labels=c(0.1,0.3,0.5,0.7,0.9))+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
xlab(expression(paste("Cheating frequency (",Omega,")")))+ylab(expression(paste("Proportion of cheaters (",bar(Delta),")")))+scale_fill_gradientn(colours=scales::viridis_pal(option="turbo")(100)[1:80])+
labs(fill=expression(paste("Optimal ",Psi)))+coord_fixed(expand=F)+geom_point(data=subset(theory,pers==maxi),aes(x=prop_cheating,y=prop_cheaters),size=2,col="white",fill="black",shape=21)+scale_alpha(range=c(0,1),guide = 'none')+
geom_point(data=point,aes(x=cheat,y=prop_cheaters,fill=prop_innovative),color="white",shape=21)













#### HYPOTHESE FOUR: CHEATING FREQUENCY IS HIGHER IN NESTED NETWORK LOW
b2$NODF=NA
for(i in b2$site){
bidon=subset(b,site==i)
bidon$tot=bidon$mut+bidon$cheat
mat=dcast(bidon,plant_species~hummingbird_species,value.var="tot",fill=0)
b2$NODF[b2$site==i]=networklevel(mat[,-1],index="weighted NODF")
}

ggplot(data=b2,aes(x=NODF,y=cheat,color=site))+geom_point()




