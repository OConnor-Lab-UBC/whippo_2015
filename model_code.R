library(nlme)
library(lme4)
library(ggplot2)
library(MuMIn)
library(vegan)

#read in data

data<-read.csv("plot_data.csv")
site_data<-read.csv("site_avgs.csv")
commdata<-read.csv("comm_data.csv")
siteIDs<-read.csv("site_ID.csv")
s.pool<-read.csv("specpool_output.csv")
fdata<-read.csv("fdata.csv")

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

abund.sel<-model.sel(mem.a1,mem.a2,mem.a3)

mem.a3.lmer<-lmer(Log_Abundance~Time+Fw_dist_km+Fw_dist_km*Time+(Time|Site),data=pdata)
#mem.a3.lmer.ci<-confint(mem.a3.lmer)

mem.a3.all.lmer<-lmer(Log_Abundance~Time+Fw_dist_km+Fw_dist_km*Time+(Time|Site),data=data)
#mem.a3.all.lmer.ci<-confint(mem.a3.all.lmer)

mem.b1<-lme(log_ENS~Time,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.b2<-lme(log_ENS~Time+Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.b3<-lme(log_ENS~Time+Fw_dist_km+Time*Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.b4<-lme(log_ENS~Time+Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")

ens.sel<-model.sel(mem.b1,mem.b2,mem.b3)

mem.b1.lmer<-lmer(log_ENS~Time+(Time|Site),data=pdata)
mem.b1.lmer.all<-lmer(log_ENS~Time+(Time|Site),data=data)

#mem.b1.lmer.ci<-confint(mem.b1.lmer)
#mem.b1.lmer.all.ci<-confint(mem.b1.lmer.all)

mem.c1<-lme(Rarefied_Richness~Time,data=pdata,random=~1|Site,control=ctrl,method="ML")
mem.c2<-lme(Rarefied_Richness~Time+Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.c3<-lme(Rarefied_Richness~Time+Fw_dist_km+Time*Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.c4<-lme(Rarefied_Richness~Time+Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")

rich.sel<-model.sel(mem.c1,mem.c2,mem.c3)

mem.c2.lmer<-lmer(Rarefied_Richness~Time+Fw_dist_km+(Time|Site),data=pdata)
mem.c2.lmer.all<-lmer(Rarefied_Richness~Time+Fw_dist_km+(Time|Site),data=data)
mem.c2.lmer.all.ci<-confint(mem.c2.lmer.all)

#mem.c2.lmer.ci<-confint(mem.c2.lmer)
#mem.c2.lmer.all.ci<-confint(mem.c2.lmer.all)

mem.d1<-lme(Simpson~Time,data=pdata,random=~1|Site,control=ctrl,method="ML")
mem.d2<-lme(Simpson~Time+Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.d3<-lme(Simpson~Time+Fw_dist_km+Time*Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.d4<-lme(Simpson~Time,data=pdata,random=~Time|Site,control=ctrl,method="ML")

simp.sel<-model.sel(mem.d2,mem.d3,mem.d4)

mem.d3.lmer<-lmer(Simpson~Time+Fw_dist_km+Time*Fw_dist_km+(Time|Site),data=pdata)
#mem.d3.lmer.ci<-confint(mem.d3.lmer)

mem.d3.lmer.all<-lmer(Simpson~Time+Fw_dist_km+Time*Fw_dist_km+(Time|Site),data=data)
#mem.d3.lmer.all.ci<-confint(mem.d3.lmer.all)

#calculate shannon diversity to run more mems
sdiv<-diversity(commdata,index="shannon")

mem.e1<-lme(Shannon~Time,data=pdata,random=~1|Site,control=ctrl,method="ML")
mem.e2<-lme(Shannon~Time+Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.e3<-lme(Shannon~Time+Fw_dist_km+Time*Fw_dist_km,data=pdata,random=~Time|Site,control=ctrl,method="ML")
mem.e4<-lme(Shannon~Time,data=pdata,random=~Time|Site,control=ctrl,method="ML")

shan.sel<-model.sel(mem.e2,mem.e3,mem.e4)

mem.e2.lmer<-lmer(Shannon~Time+Fw_dist_km+(Time|Site),data=pdata)
#mem.e2.lmer.ci<-confint(mem.e2.lmer)

mem.e2.lmer.all<-lmer(Shannon~Time+Fw_dist_km+(Time|Site),data=data)
#mem.e2.lmer.all.ci<-confint(mem.e2.lmer.all)
#plotting the data with fitted lines against distance from fw

attach(site_data)
ens_limits<-aes(ymax =ENS+ENS_std, ymin=ENS-ENS_std)
abund_limits<-aes(ymax=Log_Abundance+Log_Abundance_std,ymin=Log_Abundance-Log_Abundance_std)
rich_limits<-aes(ymax=Rarefied_Richness+Richness_std,ymin=Rarefied_Richness-Richness_std)
detach(site_data)

ens_a<-ggplot(site_data,aes(Fw_dist_km,log_ENS))
ens_b<-ens_a+geom_point(aes(shape=factor(Time)),size=8)+
  theme_bw()+geom_abline(intercept = 0.309, slope = -0.009498,linetype=1,size=2)+
  geom_abline(intercept=0.54828828,slope=0.01125,linetype=2,size=2)+
  geom_abline(intercept=0.23882,slope=-0.013589,linetype=3,size=2)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"))+
  xlab("Distance from Alberni Inlet")+ylab("log(ENS)")+scale_x_reverse()

simp_a<-ggplot(site_data,aes(Fw_dist_km,Simpson))
simp_b<-simp_a+geom_point(aes(shape=factor(Time)),size=8)+
  theme_bw()+geom_abline(intercept = 0.490, slope = -0.00610,linetype=1,size=2)+
  geom_abline(intercept=.67692867,slope=0.01201652,linetype=2,size=2)+
  geom_abline(intercept=0.71637,slope=0.00682569,linetype=3,size=2)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"))+
  xlab("Distance from Alberni Inlet")+ylab("Simpson Index")+scale_x_reverse()

shan_a<-ggplot(site_data,aes(Fw_dist_km,Shannon))
shan_b<-shan_a+geom_point(aes(shape=factor(Time)),size=8)+
  theme_bw()+geom_abline(intercept = 0.971, slope = -0.017967,linetype=1,size=2)+
  geom_abline(intercept=1.46245265,slope=0.02889799,linetype=2,size=2)+
  geom_abline(intercept=0.54561232,slope=-0.04486018,linetype=3,size=2)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"))+
  xlab("Distance from Alberni Inlet")+ylab("Shannon Index")+scale_x_reverse()

abund_a<-ggplot(site_data,aes(Fw_dist_km,Log_Abundance))
abund_b<-abund_a+geom_point(aes(shape=factor(Time)),size=8)+
  theme_bw()+geom_abline(intercept = 1.44, slope = -0.00734,linetype=1,size=2)+
  geom_abline(intercept=1.84038333,slope=-0.01871312,linetype=2,size=2)+
  geom_abline(intercept=1.526,slope=-0.0629,linetype=3,size=2)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"))+
  xlab("Distance from Alberni Inlet")+ylab("log(Abundance)")+scale_x_reverse()


rich_a<-ggplot(site_data,aes(Fw_dist_km,Rarefied_Richness))
rich_b<-rich_a+geom_point(aes(shape=factor(Time)),size=8)+
  theme_bw()+geom_abline(intercept = 6.955, slope = -0.102,linetype=1,size=2)+
  geom_abline(intercept=11.13474420,slope=0.03633219,linetype=2,size=2)+
  geom_abline(intercept=11.14,slope=-0.136,linetype=3,size=2)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"))+
  xlab("Distance from Alberni Inlet")+ylab("Rarefied species richness")+scale_x_reverse()

#calculate best fit lines

sdata_A<-subset(site_data,Time==1)
sdata_C<-subset(site_data,Time==2.5)
sdata_C<-subset(sdata_C,Site!="EI")
sdata_C<-subset(sdata_C,Site!="BE")
sdata_C<-subset(sdata_C,Site!="CC")
sdata_C<-subset(sdata_C,Site!="BI")
sdata_E<-subset(site_data,Time==4)
#use command coef(lm(ENS~Fw_dist_km,data=sdata_A))

#test for correlations in the environmental data with position in the estuary

temp.1<-lm(avg_temp~Time, data=site_data)
temp.2<-lm(avg_temp~Time+Fw_dist_km,data=site_data)
temp.3<-lm(avg_temp~Fw_dist_km+Time+Fw_dist_km*Time,data=site_data)

temp.sel<-model.sel(temp.1,temp.2,temp.3)
temp.3.ci<-confint(temp.3)

sal.1<-lm(avg_sal~Time, data=site_data)
sal.2<-lm(avg_sal~Time+Fw_dist_km,data=site_data)
sal.3<-lm(avg_sal~Fw_dist_km+Time+Fw_dist_km*Time,data=site_data)

sal.sel<-model.sel(sal.1,sal.2,sal.3)
sal.2.ci<-confint(sal.2)

dens.1<-lm(avg_dens~Time, data=site_data)
dens.2<-lm(avg_dens~Time+Fw_dist_km,data=site_data)
dens.3<-lm(avg_dens~Fw_dist_km+Time+Fw_dist_km*Time,data=site_data)

dens.sel<-model.sel(dens.1,dens.2,dens.3)
dens.2.ci<-confint(dens.2)

LAI.1<-lm(avg_LAI~Time, data=site_data)
LAI.2<-lm(avg_LAI~Time+Fw_dist_km,data=site_data)
LAI.3<-lm(avg_LAI~Fw_dist_km+Time+Fw_dist_km*Time,data=site_data)

LAI.sel<-model.sel(LAI.1,LAI.2,LAI.3)
LAI.2.ci<-confint(LAI.2)

fish_ens<-lm(fish_ens~Fw_dist_km,data=site_data)

fish_ens.ci<-confint(fish_ens)
#calculate gamma for each meadow

gvalues<-specpool(commdata,siteIDs$Code)

#plot gamma values

gam.a<-subset(s.pool,Time=="A")
gam.c<-subset(s.pool,Time=="C")
gam.e<-subset(s.pool,Time=="E")

gamma_a1<-ggplot(s.pool[1:19,],aes(Fw_dist_km,chao))
gamma_b1<-gamma_a1+geom_point(aes(shape=factor(Time)),size=8)+
  theme_bw()+geom_abline(intercept = 33.07, slope = 0.26486,linetype=1,size=2)+
  geom_abline(intercept=22.60813254,slope=-0.04468044,linetype=2,size=2)+
  geom_abline(intercept=47.31,slope=0.5565,linetype=3,size=2)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"))+
  ylab("Species richness (Chao 2)")+xlab("Distance from Alberni Inlet")+scale_x_reverse()

gamma_a2<-ggplot(s.pool[1:19,],aes(Fw_dist_km,jack2))
gamma_b2<-gamma_a2+geom_point(aes(shape=factor(Time)),size=5)+
  theme_bw()+geom_abline(intercept = 33.07, slope = -0.26486,linetype=1,size=2)+
  geom_abline(intercept=19.011,slope=0.19367,linetype=2,size=2)+
  geom_abline(intercept=47.31,slope=0.5565,linetype=3,size=2)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"))+
  ylab("Species richness (2nd order jackknife)")+xlab("Distance from Alberni Inlet")

#bust out some linear models

gam1<-lm(chao~Time,data=s.pool)
gam2<-lm(chao~Time+Fw_dist_km,data=s.pool)
gam3<-lm(chao~Time+Fw_dist_km+Time*Fw_dist_km,data=s.pool)

chao.sel<-model.sel(gam1,gam2,gam3)
gam1_ci<-confint(gam1)

gam4<-lm(jack2~Time,data=s.pool)
gam5<-lm(jack2~Time+Fw_dist_km,data=s.pool)
gam6<-lm(jack2~Time+Fw_dist_km+Time*Fw_dist_km,data=s.pool)

jack.sel<-model.sel(gam4,gam5,gam6)

gam7<-lm(chao~Msize_m2,data=s.pool)
gam8<-lm(jack2~Msize_m2,data=s.pool)
gam9<-lm(jack2~Msize_m2+Time+Time*Msize_m2,data=s.pool)


#lets work in some fish data

#calculate ENS

nm.spec<-ncol(fdata)
nm.row<-nrow(fdata)

all_ens<-numeric(0)

for(i in 1:nm.row)
{
  total<-rowSums(fdata[i,])
  sum_p2<-0
  
  for(j in 1:nm.spec)
  {
    p2<-(fdata[i,j]/total)*(fdata[i,j]/total)
    sum_p2<-sum_p2+p2    
  }
  ens<-1/sum_p2
  all_ens<-c(all_ens,ens)
}


fish.1<-lm(Log_Abundance~fish_abund,data=site_data)
fish.2<-lm(fish_abund~Fw_dist_km,data=site_data)
fish.3<-lm(fish_ens~Fw_dist_km,data=site_data)

#plot temp data

temp_1<-ggplot(site_data,aes(Fw_dist_km,avg_temp))
temp_2<-temp_1+geom_point(aes(shape=factor(Time)),size=5)+
  theme_bw()+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"))+
  ylab("Temp")+xlab("Distance from Alberni Inlet")
