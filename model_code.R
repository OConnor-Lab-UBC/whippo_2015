library(nlme)

#read in data

data<-read.csv("plot_data.csv")

#remove sites with only one dataset

pdata<-subset(data,Site!="EI")
pdata<-subset(pdata,Site!="CC")
pdata<-subset(pdata,Site!="BI")
pdata<-subset(pdata,Site!="BE")




nondir.model1<-lme(Abundance~Time,data=pdata,random=~1|Site)
dir.model1<-lme(Abundance~Time+Fw_dist_km,data=pdata,random=~1|Site)
dir.model2<-lme(Abundance~Time+Fw_dist_km+Time*Fw_dist_km,data=pdata,random=~1|Site)

new.model<-lm(Abundance~Time+Fw_dist_km+Time*Fw_dist_km,data=pdata)