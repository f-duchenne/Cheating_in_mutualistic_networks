###########################################
library(bipartite)
library(plot3D)
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
library(glmmTMB)
library(ggnewscale)
library(geiger)
library(ggthemes)
library(ggplotify)
library(randomForest)
library(ggtern)
library(ks)
library(sp)
# Set the working directory
setwd(dir="C:/Users/Duchenne/Documents/cheating")

# Load data from CSV files
dat=fread("empirical_data.csv",na.strings = c(""," ", NA))
tr2=fread("traits_data.csv",na.strings = c(""," ", NA))
sites=fread("table_s2.csv",na.strings = c(""," ", NA))

subset(dat,!is.na(hummingbird_species)) %>% group_by(Country) %>% count()

# Data manipulation and aggregation
# Calculate various metrics and create data subsets
bi=dat %>% group_by(plant_species,site,Country) %>% mutate(duration_plant=sum(duration_sampling_hours[!duplicated(waypoint)],na.rm=T)) #keep camera without interaction to calculate sampling pressure
b=subset(bi,!is.na(hummingbird_species)) %>%  # remove camera without interaction for next parts
group_by(plant_species,hummingbird_species,site,Country) %>% summarise(nb_inter=length(date[piercing=="no"]),duration=unique(duration_plant),nb_cheat=length(date[piercing=="yes"])) 
b$value=(b$nb_cheat+b$nb_inter)/(b$duration/24) # interaciton frequency
b$cheating=b$nb_cheat/(b$nb_cheat+b$nb_inter) #proportion of cheating for each pariwise interaction
b$cheat=(b$nb_cheat)/(b$duration/24) #frequency of cheating for each pariwise interaction
b$mut=(b$nb_inter)/(b$duration/24) #frequency of mutualism for each pariwise interaction
b=b %>% group_by(hummingbird_species,Country) %>% mutate(total_inter=sum(cheat+mut)) #total number of interaction per hummingbird species

b=subset(b,!is.na(value) & !is.na(plant_species) & !is.na(hummingbird_species)) # keep only interaction where both partners are identified

################################################################################################################### FIG. 4, PANEL B

################################################################################################################### FIG. 4, PANEL C
##### HYPOTHESE TWO: CHEATING IS INNOVATIVE
#summarise useful metric for that part:
b2=subset(b,!is.na(hummingbird_species)) %>% group_by(hummingbird_species,plant_species,total_inter,Country) %>% summarise(dens_mut=sum(mut),dens_cheat=sum(cheat),nrows=sum(nb_inter+nb_cheat))
b2=b2 %>% group_by(hummingbird_species,Country) %>% mutate(some_mut=sum(dens_mut),some_cheat=sum(dens_cheat))

#INFER missing tube_length
bidon=subset(tr2,!is.na(Anther_length) & !is.na(Stigmalength) & !is.na(Tube_length) &
!is.na(plant_genus) & !is.na(plant_family) & !is.na(Type))
model_tube=randomForest(Tube_length~plant_family+plant_genus+Type+Anther_length+Stigmalength,data=bidon)
tr2$Tube_length_pre=predict(model_tube,newdata=tr2)
plot(Tube_length_pre~Tube_length,data=tr2)
cor(tr2$Tube_length_pre,tr2$Tube_length,use="complete.obs")
tr2$Tube_length_pre[!is.na(tr2$Tube_length)]=tr2$Tube_length[!is.na(tr2$Tube_length)]
nrow(subset(tr2,!is.na(Tube_length_pre) & is.na(Tube_length)))
#INFER missing curvature
bidon=subset(tr2,!is.na(Anther_length) & !is.na(Stigmalength) & !is.na(Type)  & !is.na(Curvature_middle))
model_curv=randomForest(Curvature_middle~plant_family+Type+plant_genus+Anther_length+Stigmalength,data=bidon)
tr2$Curvature_middle_pre=predict(model_curv,newdata=tr2)
plot(Curvature_middle_pre~Curvature_middle,data=tr2)
cor(tr2$Curvature_middle_pre,tr2$Curvature_middle,use="complete.obs")
tr2$Curvature_middle_pre[!is.na(tr2$Curvature_middle)]=tr2$Curvature_middle[!is.na(tr2$Curvature_middle)]
nrow(subset(tr2,!is.na(Curvature_middle_pre) & is.na(Curvature_middle)))

#merging traits and data
b2=merge(b2,tr2,by=c("plant_species"),all.x=T,all.y=F)
length(unique(b2$plant_species[is.na(b2$Curvature_middle_pre) | is.na(b2$Tube_length_pre)]))
length(unique(b2$plant_species))

b2$innovative=NA

#INFERRING THE INNOVATIVE CHARACTERISTIC OF CHEATING USING A TWO DIMENSION NICHE
list_cheaters=unique(subset(b2,some_cheat>0 & !is.na(Tube_length) & !is.na(Curvature_middle_pre) & !is.na(Tube_length))[,c("hummingbird_species","Country")])
resf=NULL
list_plot=list()
for(i in 1:nrow(list_cheaters)){
bidon=subset(b2,hummingbird_species==list_cheaters$hummingbird_species[i] & !is.na(Tube_length) & Country==list_cheaters$Country[i] & !is.na(Curvature_middle_pre))
if(nrow(bidon)>2){

weight1=bidon$dens_mut/bidon$some_mut
weight1[is.na(weight1)]=0
weight2=bidon$dens_cheat/bidon$some_cheat
weight2[is.na(weight2)]=0

if(sum(weight2)>0){
#MUTUALISTIC NICHE OVER ALL SITES
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

#CALCULATE DENSITY FOR CHEATING ONLY PER SITE
obj=kde(bidon[,c("Tube_length_pre","Curvature_middle_pre")],
w=weight2,gridsize=50,H=mati,xmin=c(0,0),xmax=c(max(b2$Tube_length_pre,na.rm=T),max(b2$Curvature_middle_pre,na.rm=T)))
resi$dens_cheat=melt(obj$estimate)$value
resi$hummingbird_species=list_cheaters$hummingbird_species[i]
resi$Country=list_cheaters$Country[i]
resi$div=unique(bidon$div)
resi$total_inter=unique(bidon$total_inter)
resi$site=list_cheaters$site[i]


##### GRAPHICS
list_plot[[(length(list_plot)+1)]]=ggplot()+geom_tile(data=resi,aes(x=Tube_length,y=Curvature_middle,fill=dens_mut))+
scale_fill_gradientn(colours = terrain.colors(10))+
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

resf$mini=apply(resf[,c("dens_mut","dens_cheat")],1,min)# first step to calculate the overlap between mutualistic niche and cheating (calculate minimum)
# second step to calculate the overlap between mutualistic niche and cheating (integrate)
b2f=resf %>% dplyr::group_by(hummingbird_species,total_inter,Country) %>%
dplyr::summarise(prop_innovative=1-sum(mini)*(Tube_length[2]-Tube_length[1])*(obj$eval.points[[2]][2]-obj$eval.points[[2]][1]))


################################################################################################################### FIG. 4, PANEL D
b4=b %>% group_by(site,Country) %>% summarise(cheat=sum(cheat)/(sum(cheat)+sum(mut)),nbsp_tot=length(unique(hummingbird_species))+length(unique(plant_species)))

bff=merge(unique(b[,c("hummingbird_species","site","Country")]),b2f,by=c("hummingbird_species","Country"))
bff=merge(bff,b4,by=c("site","Country"))

bf=bff %>% group_by(site,Country, cheat,nbsp_tot,) %>% summarise(inno=mean(prop_innovative))

#### HYPOTHESE ONE: CHEATERS ARE SPECIALISTS
#summarise useful metric for that part:
b1=b %>% group_by(hummingbird_species,site,Country) %>% summarise(div=length(unique(plant_species[mut>0])),cheat=sum(cheat)/(sum(cheat)+sum(mut)),div2=vegan::diversity(mut),total_inter=sum(cheat)+sum(mut))
b1$divlog=log(b1$div+1) #log transform diversity

sites$scenario="specialists"
for(i in 1:nrow(sites)){
#test hypothesis one
model=glmmTMB(cheat~divlog,family="binomial",data=subset(b1,site==sites$site[i]))
if(fixef(model)$cond[2]>0){
sites$scenario[i]="generalists"
}
}

bf=merge(bf,sites,by=c("site","Country"))

setwd(dir="C:/Users/Duchenne/Documents/cheating")
fwrite(bf,"cheating_param_pres.csv")





#predict
pre=ggpredict(model,c("divlog[all]","Country"))
pre[pre$x>max(b1$divlog[b1$Country=="Costa-Rica"]) & pre$group=="Costa-Rica",c("predicted","conf.low","conf.high")]=NA
pre$group=factor(pre$group,levels=c("Ecuador","Costa-Rica"))
b1$Country=factor(b1$Country,levels=c("Ecuador","Costa-Rica"))

#plot the predicts and data
ggplot(data=b1,aes(x=divlog,y=cheat,size=total_inter))+
geom_point()+
stat_smooth(method="glm",method.args=list(family="quasibinomial"),col="black")+
theme_bw()+theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
legend.position="none",panel.border = element_blank(),axis.line= element_line())+ggtitle("b")+
xlab("Number of partners")+ylab("Frequency of cheating")+scale_y_continuous(breaks=c(0,0.5,1),labels=c(0,0.5,1),
limits = c(0,1))+scale_x_continuous(breaks=log(c(0,5,10,20,35)+1),labels=c(0,5,10,20,35))+
scale_size(range = c(1,2))