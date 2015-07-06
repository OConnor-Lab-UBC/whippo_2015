library(nlme)
library(lme4)
library(ggplot2)
library(MuMIn)

#read in data

data<-read.csv("plot_data.csv")
site_data<-read.csv("site_avgs.csv")

data$Site = factor(data$Site,levels(data$Site)[c(5,9,1,8,6,2,3,7,4)])

site_data$Site = factor(site_data$Site,levels(site_data$Site)[c(5,9,1,8,6,2,3,7,4)])


#remove sites with only one dataset

pdata<-subset(data,Site!="EI")
pdata<-subset(pdata,Site!="CC")
pdata<-subset(pdata,Site!="BI")
pdata<-subset(pdata,Site!="BE")

sdata<-subset(site_data,Site!="EI")
sdata<-subset(sdata,Site!="CC")
sdata<-subset(sdata,Site!="BI")
sdata<-subset(sdata,Site!="BE")

#try a mixed effects model in which Site is a random effect and Time/Distance are
#fixed effects

ctrl <- lmeControl(opt='optim')

mem.a1<-lme(Log_Abundance~Time,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.a2<-lme(Log_Abundance~Time+Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.a3<-lme(Log_Abundance~Time+Fw_dist_km+Time*Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.a3.norand<-lme(Log_Abundance~Time+Fw_dist_km+Time*Fw_dist_km,data=pdata,random=~1|Site,control=ctrl,method="ML")
mem.a3.all<-lme(Log_Abundance~Time+Fw_dist_km+Time*Fw_dist_km,data=data,random=~Time|Site,control=ctrl,method="ML")

mem.a3.lmer<-lmer(Log_Abundance~Time+Fw_dist_km+Fw_dist_km*Time+(Time|Site),data=pdata)
mem.a3.lmer.ci<-confint(mem.a3.lmer)

mem.a3.all.lmer<-lmer(Log_Abundance~Time+Fw_dist_km+Fw_dist_km*Time+(Time|Site),data=data)
mem.a3.all.lmer.ci<-confint(mem.a3.lmer)

abund.sel<-model.sel(mem.a1,mem.a2,mem.a3)

mem.b1<-lme(log_ENS~Time,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.b2<-lme(log_ENS~Time+Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.b3<-lme(log_ENS~Time+Fw_dist_km+Time*Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.b4<-lme(log_ENS~Time+Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")

ens.sel<-model.sel(mem.b1,mem.b2,mem.b3)

mem.b3.lmer<-lmer(log_ENS~Time+Fw_dist_km+Time*Fw_dist_km+(Time|Site),data=pdata)
mem.b3.lmer.all<-lmer(log_ENS~Time+Fw_dist_km+Time*Fw_dist_km+(Time|Site),data=data)

mem.b3.lmer.ci<-confint(mem.b3.lmer)
mem.b3.lmer.all.ci<-confint(mem.b3.lmer.all)

mem.c1<-lme(Rarefied_Richness~Time,data=pdata,random=~1|Site,control=ctrl,method="ML")
mem.c2<-lme(Rarefied_Richness~Time+Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.c3<-lme(Rarefied_Richness~Time+Fw_dist_km+Time*Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.c4<-lme(Rarefied_Richness~Time+Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")

rich.sel<-model.sel(mem.c1,mem.c2,mem.c3)

mem.c2.lmer<-lmer(Rarefied_Richness~Time+Fw_dist_km+(Time|Site),data=pdata)
mem.c2.lmer.all<-lmer(Rarefied_Richness~Time+Fw_dist_km+(Time|Site),data=data)

mem.c2.lmer.ci<-confint(mem.c2.lmer)
mem.c2.lmer.all.ci<-confint(mem.c2.lmer.all)
#try rerunning these models with all sites included

mem.a11<-lme(Log_Abundance~Time+Fw_dist_km+Time*Fw_dist_km,data=data,random=~1|Site,control=ctrl)
mem.a21<-lme(Log_Abundance~Time+Fw_dist_km+Time*Fw_dist_km,data=data,random=~Time|Site,control=ctrl)

mem.b11<-lme(log_ENS~Time+Fw_dist_km+Time*Fw_dist_km,data=data,random=~1|Site,control=ctrl)
mem.b21<-lme(log_ENS~Time+Fw_dist_km+Time*Fw_dist_km,data=data,random=~Time|Site,control=ctrl)

mem.c11<-lme(Rarefied_Richness~Time+Fw_dist_km+Time*Fw_dist_km,data=data,random=~1|Site,control=ctrl)
mem.c21<-lme(Rarefied_Richness~Time+Fw_dist_km+Time*Fw_dist_km,data=data,random=~Time|Site,control=ctrl)

#plotting the data with fitted lines against distance from fw

attach(site_data)
ens_limits<-aes(ymax =ENS+ENS_std, ymin=ENS-ENS_std)
abund_limits<-aes(ymax=Log_Abundance+Log_Abundance_std,ymin=Log_Abundance-Log_Abundance_std)
rich_limits<-aes(ymax=Rarefied_Richness+Richness_std,ymin=Rarefied_Richness-Richness_std)
detach(site_data)

ens_a<-ggplot(site_data,aes(Fw_dist_km,ENS))
ens_b<-ens_a+geom_point(aes(colour=factor(Time)),size=5)+
  theme_bw()+geom_abline(intercept = 1.866, slope = 0.07608,colour="pink",size=2)+
  geom_abline(intercept=3.259,slope=-0.05315,colour="green",size=2)+
  geom_abline(intercept=1.6676,slope=0.08664,colour="blue",size=2)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))
  #geom_errorbar(ens_limits)

abund_a<-ggplot(site_data,aes(Fw_dist_km,Log_Abundance))
abund_b<-abund_a+geom_point(aes(colour=factor(Time)),size=5)+
  theme_bw()+geom_abline(intercept = 1.44, slope = 0.00734,colour="pink",size=2)+
  geom_abline(intercept=1.675,slope=0.0292,colour="green",size=2)+
  geom_abline(intercept=1.526,slope=0.0629,colour="blue",size=2)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))

  #geom_errorbar(abund_limits)

rich_a<-ggplot(site_data,aes(Fw_dist_km,Rarefied_Richness))
rich_b<-rich_a+geom_point(aes(colour=factor(Time)),size=5)+
  theme_bw()+geom_abline(intercept = 6.955, slope = 0.102,colour="pink",size=2)+
  geom_abline(intercept=11.13,slope=-0.0363,colour="green",size=2)+
  geom_abline(intercept=11.14,slope=0.136,colour="blue",size=2)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))
  #geom_errorbar(rich_limits)

#calculate best fit lines

sdata_A<-subset(site_data,Time=="A")
sdata_C<-subset(site_data,Time=="C")
sdata_E<-subset(site_data,Time=="E")
#use command coef(lm(ENS~Fw_dist_km,data=sdata_A))