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

#exclude site without cheating:
res=subset(res,!(site %in% c("AMIG","MILL")))

setwd(dir="C:/Users/Duchenne/Documents/cheating/")
fwrite(unique(res),"equilibriums_analyse_empir.txt")



setwd(dir="C:/Users/Duchenne/Documents/cheating/eq_empir_null/")
lili=list.files()
res=NULL
for(i in lili){
#bidon=fread(paste0("netcar_",i,".txt"))
bidon=fread(i)
res=rbind(res,bidon)
}

#exclude site without cheating:
res=subset(res,!(site %in% c("AMIG","MILL")))
#exclude alpha >1.4
res=subset(res,efficience<2 & cost<0.1)

setwd(dir="C:/Users/Duchenne/Documents/cheating/")
fwrite(unique(res),"equilibriums_analyse_empir_null.txt")