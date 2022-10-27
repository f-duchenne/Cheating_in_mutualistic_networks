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

lili=list.files()
res=NULL
for(i in lili){
#bidon=fread(paste0("netcar_",i,".txt"))
bidon=fread(i)
bidon2=merge(as.data.frame(bidon[1,which(!(names(bidon) %in% c("prop_cheaters","prop_innovative","cost","random"))),with=F]),
expand.grid(prop_cheaters=unique(bidon$prop_cheaters),prop_cheating=0,prop_innovative=unique(bidon$prop_innovative),cost=unique(bidon$cost),random=unique(bidon$random)),
by=c("prop_cheating"),all=T)
bidon2=bidon2[,names(bidon)]
res=rbind(res,rbind(bidon2,bidon[-1,]))
}


bidon=subset(res,prop_cheaters==0 & prop_innovative==0 & cost==0 & prop_cheating==0 & random==1)

bidon$shape=sqrt(bidon$shape_p*bidon$shape_a)
plot(pers_ini~shape,data=bidon)

model=randomForest(pers_ini~shape_p+shape_a+asy+nbsp_p_dep+nbsp_a_dep,data=bidon,importance=T)
varImpPlot(model)

res$asy=res$nbsp_p_dep/res$nbsp_a_dep
res$feas=0
res$feas[res$pers_tot==1]=1
res=res %>% dplyr::group_by(essai) %>% dplyr::mutate(pers_ini=mean(pers_tot[prop_cheaters==0]))

b=subset(res,prop_cheaters>0) %>% group_by(prop_cheating,prop_innovative,cost,prop_cheaters) %>% summarise(pers=mean(pers_tot),resilience=mean(-1*valprop,na.rm=T),feasibility=mean(feas),contrib=mean(ai_contrib_aa,na.rm=T)) 
b$prop_cheaters2=paste0("Δ = ",b$prop_cheaters)
b$cost2=paste0("Λ = ",b$cost)

ggplot(data=b,aes(x=prop_cheating,y=pers,col=prop_innovative,group=prop_innovative,linetype=scenario))+geom_line()+theme_bw()+theme(panel.grid=element_blank())+facet_grid(rows=vars(cost2),cols=vars(prop_cheaters2))+labs(color="Ψ")+
xlab("Cheating frequency of cheaters (Ω)")+ylab("Persistence")+scale_color_viridis()

ggplot(data=b,aes(x=prop_cheating,y=resilience,col=prop_innovative,group=prop_innovative))+geom_line()+theme_bw()+theme(panel.grid=element_blank())+facet_grid(rows=vars(cost),cols=vars(prop_cheaters))

ggplot(data=b,aes(x=pers,y=resilience,col=prop_innovative))+geom_point()+theme_bw()+theme(panel.grid=element_blank())+facet_grid(rows=vars(cost),cols=vars(prop_cheaters))



ggplot(data=b,aes(x=prop_cheating,y=contrib,col=prop_innovative,group=prop_innovative))+geom_line()+theme_bw()+theme(panel.grid=element_blank())+facet_grid(rows=vars(cost),cols=vars(prop_cheaters))+
xlab("Cheating frequency of cheaters (Ω)")+ylab("Persistence")+scale_color_viridis()