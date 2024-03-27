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

tab=fread("data_for_figure_S9.csv")

#FIG. S9
ggplot(data=tab,aes(x=n.x,y=n.y,color=tubular))+geom_point()+facet_wrap(~site,ncol=3)+
scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),
n.breaks =2,minor_breaks=NULL)+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),
n.breaks =2,minor_breaks=NULL)+
theme_bw()+ annotation_logticks()+theme(panel.grid=element_blank())+xlab("Abundance on site")+
ylab("Sampling pressure (number of sampling events)")



