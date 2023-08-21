library(data.table)
library(dplyr)
setwd(dir="C:/Users/Duchenne/Documents/cheating/eq_empir/")

lili=list.files()
res=NULL
for(i in lili){
#bidon=fread(paste0("netcar_",i,".txt"))
bidon=fread(i)
res=rbind(res,bidon)
}

setwd(dir="C:/Users/Duchenne/Documents/cheating/")
fwrite(unique(res),"equilibriums_analyse_empir.txt")

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
library(grid)
library(tidyverse)
source("C:/Users/Duchenne/Documents/cheating/script/theme_border.r")
setwd(dir="C:/Users/Duchenne/Documents/cheating/")

#LOAD SITE METADATARAIT
sites=fread("table_s2.csv")
resi=fread("resi.txt")

res=fread("equilibriums_analyse_empir.txt")
res=merge(res,sites,by="site")
res$resilience=-1*res$valprop

res$simulation[res$simulation=="tout"]="Complete network"
res$simulation[res$simulation=="without_cheat"]="Without cheating"
res$pers_p=res$nbsp_p/res$nbsp_p_dep
res$pers_a=res$nbsp_a/res$nbsp_a_dep


#### SELECT SIMULATIONS:
res=subset(res, efficience<=2 & cheater_cost=="Escape")

b2=merge(subset(res,simulation=="Complete network"),subset(res,simulation=="Without cheating",select=c("site","cheater_cost","interf","cost","pers_tot","efficience","essai")),
by=c("site","cheater_cost","cost","efficience","essai","interf"))
b2$pers=b2$pers_tot.x-b2$pers_tot.y

b=b2 %>% group_by(interf,cost,efficience,Country,site) %>% summarise(pers_moy=mean(pers),sde=sd(pers)/sqrt(length(pers)))
b$Country=gsub("-"," ",b$Country,fixed=T)
b$Country=factor(b$Country,levels=c("Ecuador","Costa Rica"))
b$site=factor(b$site,levels=c("AMIG","BOQU","CUSI","GIGA","LONG","MILL","NIMB","NUBE","QUEB","RIOM","SANM","TOLO","Alaspungo","Alaspungo_disturbed","Verdecocha","Yanacocha","Yanacocha_disturbed"))

b3=b2 %>% group_by(interf,cost,efficience) %>% summarise(pers_moy=mean(pers),sde=sd(pers)/sqrt(length(pers)))

dodge=0.5

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











