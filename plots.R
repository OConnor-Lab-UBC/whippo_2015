library(vegan)
library(ggplot2)
library(grid)

data<-read.csv("comm_data.csv")
plot_data<-read.csv("plot_data.csv")
site_data<-read.csv("site_totals.csv")

simp<-diversity(data, index="simpson")

rare<-estimateR(data)

plot_data$Site = factor(plot_data$Site,levels(plot_data$Site)[c(5,9,1,8,6,2,3,7,4)])

aall<-ggplot(plot_data,aes(factor(Site),Abundance))
all_abund<-aall+geom_boxplot(aes(fill=factor(Time)))+xlab("Site")+theme_bw()+theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))+ylab("log(abundance)")

dd<-ggplot_build(all_abund)

dd$data[[1]]$xmin[7]<-2.875
dd$data[[1]]$xmax[7]<-3.125
dd$data[[1]]$xmin[11]<-4.875
dd$data[[1]]$xmax[11]<-5.125
dd$data[[1]]$xmin[12]<-5.875
dd$data[[1]]$xmax[12]<-6.125
dd$data[[1]]$xmin[19]<-8.875
dd$data[[1]]$xmax[19]<-9.125

grid.draw(ggplot_gtable(dd))

eall<-ggplot(plot_data,aes(factor(Site),ENS))
all_ens<-eall+geom_boxplot(aes(fill=factor(Time)))+xlab("Site")+theme_bw()+theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))+ylab("ENS")

ff<-ggplot_build(all_ens)

ff$data[[1]]$xmin[7]<-2.875
ff$data[[1]]$xmax[7]<-3.125
ff$data[[1]]$xmin[11]<-4.875
ff$data[[1]]$xmax[11]<-5.125
ff$data[[1]]$xmin[12]<-5.875
ff$data[[1]]$xmax[12]<-6.125
ff$data[[1]]$xmin[19]<-8.875
ff$data[[1]]$xmax[19]<-9.125


grid.draw(ggplot_gtable(ff))

rrall<-ggplot(plot_data,aes(factor(Site),Rarefied_Richness))
all_rich<-rrall+geom_boxplot(aes(fill=factor(Time)))+xlab("Site")+theme_bw()+theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))+ylab("Rarefied richness")

ee<-ggplot_build(all_rich)

ee$data[[1]]$xmin[7]<-2.875
ee$data[[1]]$xmax[7]<-3.125
ee$data[[1]]$xmin[11]<-4.875
ee$data[[1]]$xmax[11]<-5.125
ee$data[[1]]$xmin[12]<-5.875
ee$data[[1]]$xmax[12]<-6.125
ee$data[[1]]$xmin[19]<-8.875
ee$data[[1]]$xmax[19]<-9.125


grid.draw(ggplot_gtable(ee))

aplot<-subset(plot_data,Time=="A")
aplot$Abundance<-log(aplot$Abundance)
aplot$Site<- factor(aplot$Site,levels=unique(as.character(aplot$Site)))
cplot<-subset(plot_data,Time=="C")
cplot$Abundance<-log(cplot$Abundance)
cplot$Site<- factor(cplot$Site,levels=unique(as.character(cplot$Site)))
eplot<-subset(plot_data,Time=="E")
eplot$Abundance<-log(eplot$Abundance)
eplot$Site<- factor(eplot$Site,levels=unique(as.character(eplot$Site)))

aa<-ggplot(aplot,aes(factor(Site),Abundance))
a_abund<-aa+geom_boxplot()+theme_bw()+xlab("Site")

ca<-ggplot(cplot,aes(factor(Site),Abundance))
c_abund<-ca+geom_boxplot()+theme_bw()+xlab("Site")

ea<-ggplot(eplot,aes(factor(Site),Abundance))
e_abund<-ea+geom_boxplot()+theme_bw()+xlab("Site")

ae<-ggplot(aplot,aes(factor(Site),ENS))
a_ens<-ae+geom_boxplot()+theme_bw()+xlab("Site")

ce<-ggplot(cplot,aes(factor(Site),ENS))
c_ens<-ce+geom_boxplot()+theme_bw()+xlab("Site")

ee<-ggplot(eplot,aes(factor(Site),ENS))
e_ens<-ee+geom_boxplot()+theme_bw()+xlab("Site")



ad<-ggplot(aplot,aes(factor(Site),Simpson))
a_div<-ad+geom_boxplot()+theme_bw()+xlab("Site")

cd<-ggplot(cplot,aes(factor(Site),Simpson))
c_div<-cd+geom_boxplot()+theme_bw()+xlab("Site")

ed<-ggplot(eplot,aes(factor(Site),Simpson))
e_div<-ed+geom_boxplot()+theme_bw()+xlab("Site")

ar<-ggplot(aplot,aes(factor(Site),Rarefied_Richness))
a_rich<-ar+geom_boxplot()+theme_bw()+xlab("Site")+ylab("Rarefied Richness")

cr<-ggplot(cplot,aes(factor(Site),Rarefied_Richness))
c_rich<-cr+geom_boxplot()+theme_bw()+xlab("Site")+ylab("Rarefied Richness")

er<-ggplot(eplot,aes(factor(Site),Rarefied_Richness))
e_rich<-er+geom_boxplot()+theme_bw()+xlab("Site")+ylab("Rarefied Richness")

dca<-radfit(site_data[1,])
dcc<-radfit(site_data[7,])
dce<-radfit(site_data[12,])

wia<-radfit(site_data[2,])
wic<-radfit(site_data[8,])
wie<-radfit(site_data[13,])

rpa<-radfit(site_data[4,])
rpc<-radfit(site_data[9,])
rpe<-radfit(site_data[15,])

nba<-radfit(site_data[5,])
nbc<-radfit(site_data[10,])
nbe<-radfit(site_data[14,])

cba<-radfit(site_data[6,])
cbc<-radfit(site_data[11,])
cbe<-radfit(site_data[16,])

clump<-dispindmorisita(data[81:96,],unique.rm=T)
val<-clump$imst
na.omit(val)
write.csv(val,"val.csv")