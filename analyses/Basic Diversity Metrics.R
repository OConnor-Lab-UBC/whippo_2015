###################################################################################
#                                                                                ##
# Basic Diversity Metrics                                                        ##
# Data are current as of 2018-03-28                                              ##
# Data source: O'Connor Lab - UBC                                                ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2018-03-28                                                        ##
#                                                                                ##
###################################################################################

# SUMMARY: 
# Basic diversity analyses used to address reviewer comments of "Epifaunal 
# diversity patterns among eelgrass meadows suggest metacommunity structure" at
# Ecosphere including rarified species richness, ENS, Shannon Diversity, Simpson
# Evenness, and Hellinger Distance Matrix


# Required Files (check that script is loading latest version):
# epicomm_201802.csv

# Associated Scripts:
# 2018-beta analysis.R 

# TO DO

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# RECENT CHANGES TO SCRIPT                                                        #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# MANIPULATE DATA                                                                 #   
#                                                                                 #
###################################################################################

###################################################################################
# RECENT CHANGES TO SCRIPT                                                        #
###################################################################################

# 2018-03-28 Forked with correct rawcomm processing code, clean up fig code
# 2018-03-26 Cleaned up for figures and added analyses from 2018-beta analysis.R
# 2018-03-24 Added MDS analysis
# 2018-03-05 Switched back to rawcomm dataset with treatment from Mary's script
# 2018-02-23 Added evenness code, and began Figure 3 panel.
# 2018-02-21 Created script. 

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

library(reshape2) # manipulate data
library(vegan) # diversity analyses
library(viridis) # color palette
library(ggpubr) # combining plots
library(lubridate) # manipulate data
library(car) # normality tests
library(lme4)
library(multcomp) #tukey post hoc
library(MuMIn) # lm analyses
library(tidyverse) # manipulate data

# function to scale hellinger matrix between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# read full community data ------------------------------------------------
data <- read.csv("./data/rawcomm.csv")
traits <- read.csv("./data/grazertraits3.csv")
sites <- read.csv("./data/site.info.csv")

# delete, add and redefine columns ----------------------------------------
traits <- traits[,-c(3,8:10)]

data$date1 <- mdy(data$date)

# data preparation --------------------------------------------------------
## Melt and recast so the datafile goes from long to wide (species as columns)
data.m <- melt(data, id = c(1,2,3,4,5,52, 53))

## Clean and correct species names. 
## in the data file, some names we initially used for taxa needed updating based on improvements in our ability to identify them. So these replacements reflect those updates. 
levels(data.m$variable)[levels(data.m$variable)== "Bittium.spp."] <- "Lirobittium.spp."
levels(data.m$variable)[levels(data.m$variable)== "Olivella.sp."] <- "Callianax.sp."
levels(data.m$variable)[levels(data.m$variable)== "Cypricercus."] <- "Cyprideis.beaconensis"
levels(data.m$variable)[levels(data.m$variable)== "Odontosyllis"] <- "Polychaete1"
# Idotea.resecata is also incorrect (should be Pentidotea.resecata, however it is kept as is to maintain functionality of code)

# Clean up time code labels so they are easier to model
levels(unique(data$Time.Code2))
levels(data.m$Time.Code)
data.m$Time.Code <- as.character(data.m$Time.Code)
data.m$Time.Code[data.m$Time.Code == "C "] <- "C"
data.m$Time.Code <- as.factor(data.m$Time.Code)
data.m$value <- as.numeric(data.m$value)
levels(data.m$Time.Code)

## Merge diversity file with site metadata
data.s <- merge(data.m, sites, by = "site")

library(plyr)
## Sum across size classes within plots (samples), takes several seconds to run
data.p <- ddply(data.s, .(site, date1, Sample, Time.Code2, variable, dfw,order.dfw,area,salinity, shoot.density, fetch.jc), summarise, sum(value))
detach(package:plyr)

data.p$time.ID <- paste(data.p$site, data.p$Time.Code2, sep = '.') #could look at finer time resolution by using Time.Code here
names(data.p) <- c("site", "Date", "Sample", "Time.Code2", "species", "dfw","order","area","salinity","shoot.density","fetch","abundance", "TimeID")

## Merge with traits and sort by taxa or functional groups
data.tr <- merge(data.p, traits[,-1], by.x = "species", by.y = "species.names", all.x = TRUE, all.y = FALSE)

## create datafile to be posted with paper: 
# write.csv(data.tr, "Whippodata.csv")

# Create datafiles for taxa and times -------------------------------------

## Remove all taxa that are not epifauna: 
levels(data.tr$eelgrss.epifauna)
epicomm_z <- data.tr %>% filter(eelgrss.epifauna == "yes" | eelgrss.epifauna == "sometimes")
epicomm_s <- epicomm_z %>%
  spread(species, abundance)

# remove unused columns
drop.cols <- c('Date', 'dfw', 'order', 'area', 'salinity', 'shoot.density', 'fetch', 'TimeID', 'function.', 'taxon', 'group', 'eelgrss.epifauna')
epicomm_s <- epicomm_s %>%
  select(-one_of(drop.cols))

# replace NAs with 0 
epicomm_s[is.na(epicomm_s)] <- 0

# summarize by sample
epicomm_s <- epicomm_s %>%
  group_by(site, Time.Code2, Sample) %>%
  summarize_all(funs(sum))



# separate into time periods for shannon and ens
EpiA <- epicomm_s %>%
  filter(Time.Code2 == "A")
EpiA$Time.Code2 <- "May"
EpiA$site <- factor(EpiA$site, levels = c("DC", "WI", "RP", "NB", "CB"))

EpiC <- epicomm_s %>%
  filter(Time.Code2 == "C")
EpiC$Time.Code2 <- "June/July"
EpiC$site <- factor(EpiC$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))

EpiE <- epicomm_s %>%
  filter(Time.Code2 == "E")
EpiE$Time.Code2 <- "August"
EpiE$site <- factor(EpiE$site, levels = c("DC", "WI", "RP", "NB", "CB"))



# separate each sitetime for analysis, have to for rarefy command and hell/bray distance

indiv_sites <- epicomm_s
# change time labels

levels(indiv_sites$Time.Code2)[levels(indiv_sites$Time.Code2)== "A"] <- "May"
levels(indiv_sites$Time.Code2)[levels(indiv_sites$Time.Code2)== "C"] <- "June/July"
levels(indiv_sites$Time.Code2)[levels(indiv_sites$Time.Code2)== "E"] <- "August"


DCA <- indiv_sites %>%
  filter(site == "DC" & Time.Code2 == "May")
DCC <- indiv_sites %>%
  filter(site == "DC" & Time.Code2 == "June/July")
DCE <- indiv_sites %>%
  filter(site == "DC" & Time.Code2 == "August")

WIA <- indiv_sites %>%
  filter(site == "WI" & Time.Code2 == "May")
WIC <- indiv_sites %>%
  filter(site == "WI" & Time.Code2 == "June/July")
WIE <- indiv_sites %>%
  filter(site == "WI" & Time.Code2 == "August")

RPA <- indiv_sites %>%
  filter(site == "RP" & Time.Code2 == "May")
RPC <- indiv_sites %>%
  filter(site == "RP" & Time.Code2 == "June/July")
RPE <- indiv_sites %>%
  filter(site == "RP" & Time.Code2 == "August")

NBA <- indiv_sites %>%
  filter(site == "NB" & Time.Code2 == "May")
NBC <- indiv_sites %>%
  filter(site == "NB" & Time.Code2 == "June/July")
NBE <- indiv_sites %>%
  filter(site == "NB" & Time.Code2 == "August")

CBA <- indiv_sites %>%
  filter(site == "CB" & Time.Code2 == "May")
CBC <- indiv_sites %>%
  filter(site == "CB" & Time.Code2 == "June/July")
CBE <- indiv_sites %>%
  filter(site == "CB" & Time.Code2 == "August")

BEC <- indiv_sites %>%
  filter(site == "BE" & Time.Code2 == "June/July")
EIC <- indiv_sites %>%
  filter(site == "EI" & Time.Code2 == "June/July")
BIC <- indiv_sites %>%
  filter(site == "BI" & Time.Code2 == "June/July")
CCC <- indiv_sites %>%
  filter(site == "CC" & Time.Code2 == "June/July")

#######################
#######################








###################################################################################
# DIVERSITY MEASURES                                                              #
###################################################################################

############### OBSERVED RICHNESS

Rich <- specnumber(epicomm_s[-c(1:3)])
obsrich <- epicomm_s[c(1:3)]
obsrich$obsrich <- Rich
obsrich$site <- factor(obsrich$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))


rich.C <- obsrich %>%
  filter(Time.Code2 == "C")
rich.C$site <- factor(rich.C$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))

prim_sites <- c("DC", "WI", "RP", "NB", "CB")
rich_prim <- obsrich %>%
  filter(site %in% prim_sites)
rich_prim$site <- factor(rich_prim$site, levels = c("DC", "WI", "RP", "NB", "CB"))
rich_prim <- rich_prim %>%
  group_by(site, Time.Code2) %>%
  summarise(mean(obsrich))
colnames(rich_prim)[which(names(rich_prim) == "mean(obsrich)")] <- "obsrich"

  
############### Richness lm across all sites in midsummer, model selection and post hoc test

richlm <- lm(rich.C$obsrich ~ rich.C$site)
rich0 <- lm(rich.C$obsrich ~ 1)
anova(richlm, rich0)
model.sel(richlm, rich0)
summary(richlm)
richaov <- aov(obsrich ~ site, data = rich.C)
richtuk <- glht(richaov, linfct = mcp(site = "Tukey"))
richtuk.cld <- cld(richtuk)
richtuk.cld <- print(richtuk.cld)




############### RAREFIED RICHNESS


DCA_rare <- DCA[,4:37] %>%
  rarefy(5)

DCA_rare <- DCA[,4:37] %>%
  rarefy(5)
DCC_rare <- DCC[,4:37] %>%
  rarefy(5)
DCE_rare <- DCE[,4:37] %>%
  rarefy(5)

WIA_rare <- WIA[,4:37] %>%
  rarefy(5)
WIC_rare <- WIC[,4:37] %>%
  rarefy(5)
WIE_rare <- WIE[,4:37] %>%
  rarefy(5)

RPA_rare <- RPA[,4:37] %>%
  rarefy(5)
RPC_rare <- RPC[,4:37] %>%
  rarefy(5)
RPE_rare <- RPE[,4:37] %>%
  rarefy(5)

NBA_rare <- NBA[,4:37] %>%
  rarefy(5)
NBC_rare <- NBC[,4:37] %>%
  rarefy(5)
NBE_rare <- NBE[,4:37] %>%
  rarefy(5)

CBA_rare <- CBA[,4:37] %>%
  rarefy(5)
CBC_rare <- CBC[,4:37] %>%
  rarefy(5)
CBE_rare <- CBE[,4:37] %>%
  rarefy(5)

BEC_rare <- BEC[,4:37] %>%
  rarefy(5)
EIC_rare <- EIC[,4:37] %>%
  rarefy(5)
BIC_rare <- BIC[,4:37] %>%
  rarefy(5)
CCC_rare <- CCC[,4:37] %>%
  rarefy(5)


# combine into single data frame FIGURE 4
# make all same length
length(DCA_rare) <- 17                      
length(DCC_rare) <- 17  
length(DCE_rare) <- 17 
length(WIA_rare) <- 17 
length(WIC_rare) <- 17  
length(WIE_rare) <- 17  
length(RPA_rare) <- 17  
length(RPC_rare) <- 17  
length(RPE_rare) <- 17  
length(NBA_rare) <- 17  
length(NBC_rare) <- 17  
length(NBE_rare) <- 17  
length(CBA_rare) <- 17  
length(CBC_rare) <- 17  
length(CBE_rare) <- 17 
# combine
rare_prim <- melt(data.frame(DCA_rare, DCC_rare, DCE_rare, WIA_rare, WIC_rare, WIE_rare, RPA_rare, RPC_rare, RPE_rare, NBA_rare, NBC_rare, NBE_rare, CBA_rare, CBC_rare, CBE_rare))
# rename columns, reduce sitetime values, and split into site and time
colnames(rare_prim) <- c('sitetime', 'value') 
rare_prim$sitetime <- as.character(rare_prim$sitetime)
rare_prim$sitetime <- substr(rare_prim$sitetime,1,3)
rare_prim <- transform(rare_prim, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
rare_prim <- rare_prim %>%
  select(-sitetime)
# reorder factors
rare_prim$site <- factor(rare_prim$site, levels = c("DC", "WI", "RP", "NB", "CB"))
# remove NA
rare_prim <- rare_prim %>%
  filter(!is.na(value))


# combine into single data frame FIGURE 2
# make all same length
length(DCC_rare) <- 17  
length(WIC_rare) <- 17  
length(BEC_rare) <- 17  
length(EIC_rare) <- 17  
length(RPC_rare) <- 17  
length(NBC_rare) <- 17  
length(CBC_rare) <- 17  
length(BIC_rare) <- 17 
length(CCC_rare) <- 17
# combine
rare_C <- melt(data.frame(DCC_rare, WIC_rare, BEC_rare, EIC_rare, RPC_rare, NBC_rare, NBE_rare, CBC_rare, BIC_rare, CCC_rare))
# rename columns, reduce sitetime values, and split into site and time
colnames(rare_C) <- c('sitetime', 'value') 
rare_C$sitetime <- as.character(rare_C$sitetime)
rare_C$sitetime <- substr(rare_C$sitetime,1,3)
rare_C <- transform(rare_C, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
rare_C <- rare_C %>%
  select(-sitetime)
# reorder factors
rare_C$site <- factor(rare_C$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))
# remove NA
rare_C <- rare_C %>%
  filter(!is.na(value))


############### ABUNDANCE

abund_all <- epicomm_s
abund_all$abundance <- rowSums(epicomm_s[4:37])  
abund_all <- abund_all[c(1:3, 38)]

abund_C <- abund_all %>%
  filter(Time.Code2 == "C")
abund_C$site <- factor(abund_C$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))

# abundance through time

abund_prim <- abund_all %>%
  filter(site %in% prim_sites) %>%
  group_by(site, Time.Code2) %>%
  summarise(mean(abundance))
colnames(abund_prim)[which(names(abund_prim) == "mean(abundance)")] <- "abundance"
abund_prim$site <- factor(abund_prim$site, levels = c("DC", "WI", "RP", "NB", "CB"))

############### SHANNON DIVERSITY & ENS


###### Figure 2

H_all <- diversity(epicomm_s[-(c(1:3))], index ="shannon")



#combine metrics
div.summary <- epicomm_s[c(1:3)]
div.summary$Shannon <- H_all
div.summary$ENS <- exp(div.summary$Shannon)


# reorder factors
div.summary$site <- factor(div.summary$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))

div.summaryC <- div.summary %>%
  filter(Time.Code2 == "C")
div.summaryC$site <- factor(div.summaryC$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))

shan_prim <- div.summary %>%
  filter(site %in% prim_sites)
shan_prim$site <- factor(shan_prim$site, levels = c("DC", "WI", "RP", "NB", "CB"))
ens_prim <- shan_prim %>%
  group_by(site, Time.Code2) %>%
  summarise(mean(ENS))
colnames(ens_prim)[which(names(ens_prim) == "mean(ENS)")] <- "ENS"
shan_time <- shan_prim %>%
  group_by(site, Time.Code2) %>%
  summarise(mean(Shannon))
colnames(shan_time)[which(names(shan_time) == "mean(Shannon)")] <- "Shannon"
shan_time$site <- factor(shan_time$site, levels = c("DC", "WI", "RP", "NB", "CB"))


############### abundance lm across all sites in midsummer, model selection and post hoc test

leveneTest(abund_C$abundance ~ abund_C$site)
abund_C$logabund <- log10(abund_C$abundance + 1)
leveneTest(abund_C$logabund ~ abund_C$site)

abunlm <- lm(abund_C$logabund ~ abund_C$site)
abun0 <- lm(abund_C$logabund ~ 1)
anova(abunlm, abun0)
model.sel(abunlm, abun0)
summary(abunlm)
abunaov <- aov(logabund ~ site, data = abund_C)
abuntuk <- glht(abunaov, linfct = mcp(site = "Tukey"))
abuntuk.cld <- cld(abuntuk)
abuntuk.cld <- print(abuntuk.cld)



############### ENS lm across all sites in midsummer, model selection and post hoc test

enslm <- lm(div.summaryC$ENS ~ div.summaryC$site)
ens0 <- lm(div.summaryC$ENS ~ 1)
anova(enslm, ens0)
model.sel(enslm, ens0)
summary(enslm)
ensaov <- aov(ENS ~ site, data = div.summaryC)
enstuk <- glht(ensaov, linfct = mcp(site = "Tukey"))
enstuk.cld <- cld(enstuk)
enstuk.cld <- print(enstuk.cld)

############### Shannon lm across all sites in midsummer, model selection and post hoc test

shanlm <- lm(div.summaryC$Shannon ~ div.summaryC$site)
shan0 <- lm(div.summaryC$Shannon ~ 1)
anova(shanlm, shan0)
model.sel(shanlm, shan0)
summary(shanlm)
shanaov <- aov(Shannon ~ site, data = div.summaryC)
shantuk <- glht(shanaov, linfct = mcp(site = "Tukey"))
shantuk.cld <- cld(shantuk)
shantuk.cld <- print(shantuk.cld)





############### STOPGAP RAREFIED ANALYSES FROM MARY

rarefied <- read.csv("./data/div.summary.csv")
rarefied <- rarefied %>%
  select(site, site.time, R2)
rarefied <- transform(rarefied, time = substr(site.time, 4, 4))

rareC <- rarefied %>%
  filter(time == "C")
rareC$site <- factor(rareC$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))

rare3 <- read.csv("./data/div.summary3times.csv")

new_rare_prim <- rare3 %>%
  filter(site %in% prim_sites)
new_rare_prim <- new_rare_prim %>%
  select(site, Time.Code2, R2)
new_rare_prim$site <- factor(new_rare_prim$site, levels = c("DC", "WI", "RP", "NB", "CB"))
new_rare_prim <- new_rare_prim %>%
  group_by(site, Time.Code2) %>%
  summarise(mean(R2))
colnames(new_rare_prim)[which(names(new_rare_prim) == "mean(R2)")] <- "R2"


############### RAREFIED lm across all sites in midsummer, model selection and post hoc test

rarelm <- lm(rareC$R2 ~ rareC$site)
rare0 <- lm(rareC$R2 ~ 1)
anova(rarelm, rare0)
model.sel(rarelm, rare0)
summary(rarelm)
rareaov <- aov(R2 ~ site, data = rareC)
raretuk <- glht(rareaov, linfct = mcp(site = "Tukey"))
raretuk.cld <- cld(raretuk)
raretuk.cld <- print(raretuk.cld)







############### HELLINGER DISTANCE

DCA_hell <- range01(vegdist(decostand(DCA[,4:37], "hellinger"), "euclidean"))
DCC_hell <- range01(vegdist(decostand(DCC[,4:37], "hellinger"), "euclidean"))
DCE_hell <- range01(vegdist(decostand(DCE[,4:37], "hellinger"), "euclidean"))

WIA_hell <- range01(vegdist(decostand(WIA[,4:37], "hellinger"), "euclidean"))
WIC_hell <- range01(vegdist(decostand(WIC[,4:37], "hellinger"), "euclidean"))
WIE_hell <- range01(vegdist(decostand(WIE[,4:37], "hellinger"), "euclidean"))

RPA_hell <- range01(vegdist(decostand(RPA[,4:37], "hellinger"), "euclidean"))
RPC_hell <- range01(vegdist(decostand(RPC[,4:37], "hellinger"), "euclidean"))
RPE_hell <- range01(vegdist(decostand(RPE[,4:37], "hellinger"), "euclidean"))

NBA_hell <- range01(vegdist(decostand(NBA[,4:37], "hellinger"), "euclidean"))
NBC_hell <- range01(vegdist(decostand(NBC[,4:37], "hellinger"), "euclidean"))
NBE_hell <- range01(vegdist(decostand(NBE[,4:37], "hellinger"), "euclidean"))

CBA_hell <- range01(vegdist(decostand(CBA[,4:37], "hellinger"), "euclidean"))
CBC_hell <- range01(vegdist(decostand(CBC[,4:37], "hellinger"), "euclidean"))
CBE_hell <- range01(vegdist(decostand(CBE[,4:37], "hellinger"), "euclidean"))

BEC_hell <- range01(vegdist(decostand(BEC[,4:37], "hellinger"), "euclidean"))
EIC_hell <- range01(vegdist(decostand(EIC[,4:37], "hellinger"), "euclidean"))
BIC_hell <- range01(vegdist(decostand(BIC[,4:37], "hellinger"), "euclidean"))
CCC_hell <- range01(vegdist(decostand(CCC[,4:37], "hellinger"), "euclidean"))


# combine into single data frame Figure 4
# make all same length
length(DCA_hell) <- 137                      
length(DCC_hell) <- 137
length(DCE_hell) <- 137 
length(WIA_hell) <- 137 
length(WIC_hell) <- 137  
length(WIE_hell) <- 137  
length(RPA_hell) <- 137  
length(RPC_hell) <- 137  
length(RPE_hell) <- 137  
length(NBA_hell) <- 137  
length(NBC_hell) <- 137  
length(NBE_hell) <- 137  
length(CBA_hell) <- 137  
length(CBC_hell) <- 137  
length(CBE_hell) <- 137
DCA_hell <- as.data.frame(DCA_hell)
DCC_hell <- as.data.frame(DCC_hell)
DCE_hell <- as.data.frame(DCE_hell)
WIA_hell <- as.data.frame(WIA_hell)
WIC_hell <- as.data.frame(WIC_hell)
WIE_hell <- as.data.frame(WIE_hell)
RPA_hell <- as.data.frame(RPA_hell)
RPC_hell <- as.data.frame(RPC_hell)
RPE_hell <- as.data.frame(RPE_hell)
NBA_hell <- as.data.frame(NBA_hell)
NBC_hell <- as.data.frame(NBC_hell)
NBE_hell <- as.data.frame(NBE_hell)
CBA_hell <- as.data.frame(CBA_hell)
CBC_hell <- as.data.frame(CBC_hell)
CBE_hell <- as.data.frame(CBE_hell)

# combine primaries
hell_prim <- melt(data.frame(DCA_hell, DCC_hell, DCE_hell, WIA_hell, WIC_hell, WIE_hell, RPA_hell, RPC_hell, RPE_hell, NBA_hell, NBC_hell, NBE_hell, CBA_hell, CBC_hell, CBE_hell))
# rename columns, reduce sitetime values, and split into site and time PRIMARY
colnames(hell_prim) <- c('sitetime', 'value') 
hell_prim$sitetime <- as.character(hell_prim$sitetime)
hell_prim$sitetime <- substr(hell_prim$sitetime,1,3)
hell_prim <- transform(hell_prim, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
hell_prim <- hell_prim %>%
  select(-sitetime) %>%
  filter(!is.na(value))
hell_prim$site <- factor(hell_prim$site, levels = c("DC", "WI", "RP", "NB", "CB"))





# combine into single data frame Figure 2
# make all same length 
length(BEC_hell) <- 137
length(EIC_hell) <- 137
length(BIC_hell) <- 137
length(CCC_hell) <- 137

BEC_hell <- as.data.frame(BEC_hell)
EIC_hell <- as.data.frame(EIC_hell)
BIC_hell <- as.data.frame(BIC_hell)
CCC_hell <- as.data.frame(CCC_hell)

# combine primaries
hell_C <- melt(data.frame(DCC_hell, WIC_hell, BEC_hell, EIC_hell, RPC_hell, NBC_hell, CBC_hell, BIC_hell, CCC_hell))
# rename columns, reduce sitetime values, and split into site and time PRIMARY
colnames(hell_C) <- c('sitetime', 'value') 
hell_C$sitetime <- as.character(hell_C$sitetime)
hell_C$sitetime <- substr(hell_C$sitetime,1,3)
hell_C <- transform(hell_C, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
hell_C <- hell_C %>%
  select(-sitetime) %>%
  filter(!is.na(value))
hell_C$site <- factor(hell_C$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))


############### Hellinger lm across all sites in midsummer, model selection and post hoc test

helllm <- lm(hell_C$value ~ hell_C$site)
hell0 <- lm(hell_C$value ~ 1)
anova(helllm, hell0)
model.sel(helllm, hell0)
summary(helllm)
hellaov <- aov(value ~ site, data = hell_C)
helltuk <- glht(hellaov, linfct = mcp(site = "Tukey"))
helltuk.cld <- cld(helltuk)
helltuk.cld <- print(helltuk.cld)




##################### MDS of community

epicomm_mds <- epicomm_s
epicomm_mds <- epicomm_mds %>%
  select(-Sample) %>%
  filter(site == c("DC", "WI", "RP", "NB", "CB")) %>%
  group_by(site, Time.Code2) %>%
  summarise_all(sum)
# reorder factor levels
epicomm_mds$site <- factor(epicomm_mds$site, levels = c("DC", "WI", "RP", "NB", "CB"))

# change time labels and prep for MDS
levels(epicomm_mds$Time.Code2)[levels(epicomm_mds$Time.Code2)== "A"] <- "May"
levels(epicomm_mds$Time.Code2)[levels(epicomm_mds$Time.Code2)== "C"] <- "June/July"
levels(epicomm_mds$Time.Code2)[levels(epicomm_mds$Time.Code2)== "E"] <- "August"

colnames(epicomm_mds)[which(names(epicomm_mds) == "Time.Code2")] <- "month"

# remove species only found at secondary sites
epicomm_mds <- epicomm_mds %>%
  select(-c(Nebalia.sp., Callianax.sp., Nemertea, Lirobittium.spp., Alia.carinata)) 

epi_mat <- epicomm_mds[,3:31]

epiMDS <- metaMDS(epi_mat)
epicomm_mds_points <- epiMDS$points
epicomm_mds_points <- data.frame(epicomm_mds_points)
plot_data_tax <- data.frame(epicomm_mds_points, epicomm_mds[,1:2])
library(plyr)
chulls_tax <- ddply(plot_data_tax, .(month), function(df) df[chull(df$MDS1, df$MDS2), ])
detach(package:plyr)



###################################################################################
# TIME SERIES BETA FIGURES                                                        #
###################################################################################

############### JACCARD DIVERSITY

DCA_jaccard <- DCA[,4:37] %>%
  vegdist(method = "jaccard")
DCC_jaccard <- DCC[,4:37] %>%
  vegdist(method = "jaccard")
DCE_jaccard <- DCE[,4:37] %>%
  vegdist(method = "jaccard")

WIA_jaccard <- WIA[,4:37] %>%
  vegdist(method = "jaccard")
WIC_jaccard <- WIC[,4:37] %>%
  vegdist(method = "jaccard")
WIE_jaccard <- WIE[,4:37] %>%
  vegdist(method = "jaccard")

RPA_jaccard <- RPA[,4:37] %>%
  vegdist(method = "jaccard")
RPC_jaccard <- RPC[,4:37] %>%
  vegdist(method = "jaccard")
RPE_jaccard <- RPE[,4:37] %>%
  vegdist(method = "jaccard")

NBA_jaccard <- NBA[,4:37] %>%
  vegdist(method = "jaccard")
NBC_jaccard <- NBC[,4:37] %>%
  vegdist(method = "jaccard")
NBE_jaccard <- NBE[,4:37] %>%
  vegdist(method = "jaccard")

CBA_jaccard <- CBA[,4:37] %>%
  vegdist(method = "jaccard")
CBC_jaccard <- CBC[,4:37] %>%
  vegdist(method = "jaccard")
CBE_jaccard <- CBE[,4:37] %>%
  vegdist(method = "jaccard")

BEC_jaccard <- BEC[,4:37] %>%
  vegdist(method = "jaccard")
BIC_jaccard <- BIC[,4:37] %>%
  vegdist(method = "jaccard")
CCC_jaccard <- CCC[,4:37] %>%
  vegdist(method = "jaccard")
EIC_jaccard <- EIC[,4:37] %>%
  vegdist(method = "jaccard")


# take mean of jaccard distance

DCA_jmean <- mean(DCA_jaccard)
DCC_jmean <- mean(DCC_jaccard)
DCE_jmean <- mean(DCE_jaccard)

WIA_jmean <- mean(WIA_jaccard)
WIC_jmean <- mean(WIC_jaccard)
WIE_jmean <- mean(WIE_jaccard)

RPA_jmean <- mean(RPA_jaccard)
RPC_jmean <- mean(RPC_jaccard)
RPE_jmean <- mean(RPE_jaccard)

NBA_jmean <- mean(NBA_jaccard)
NBC_jmean <- mean(NBC_jaccard)
NBE_jmean <- mean(NBE_jaccard)

CBA_jmean <- mean(CBA_jaccard)
CBC_jmean <- mean(CBC_jaccard)
CBE_jmean <- mean(CBE_jaccard)

# renamed time code for secondaries
BIC_jmean <- mean(BIC_jaccard)
BEC_jmean <- mean(BEC_jaccard)
CCC_jmean <- mean(CCC_jaccard)
EIC_jmean <- mean(EIC_jaccard)

all_jaccard <- melt(data.frame(DCA_jmean, DCC_jmean, DCE_jmean, WIA_jmean, WIC_jmean, WIE_jmean, RPA_jmean, RPC_jmean, RPE_jmean, NBA_jmean, NBC_jmean, NBE_jmean, CBA_jmean, CBC_jmean, CBE_jmean, BIC_jmean, BEC_jmean, CCC_jmean, EIC_jmean))
# rename columns, reduce sitetime values, and split into site and time
colnames(all_jaccard) <- c('sitetime', 'value') 
all_jaccard$sitetime <- as.character(all_jaccard$sitetime)
all_jaccard$sitetime <- substr(all_jaccard$sitetime,1,3)
all_jaccard <- transform(all_jaccard, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
all_jaccard <- all_jaccard %>%
  select(-sitetime)
# reorder factors
all_jaccard$site <- factor(all_jaccard$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))

# Primary sites only
jacc_prime <- all_jaccard %>%
  filter(site %in% prim_sites)
jacc_prime$site <- factor(jacc_prime$site, levels = c("DC", "WI", "RP", "NB", "CB"))

### MIDSUMMER ONLY JACCARD

length(DCC_jaccard) <- 137
length(WIC_jaccard) <- 137  
length(RPC_jaccard) <- 137  
length(NBC_jaccard) <- 137  
length(CBC_jaccard) <- 137  
length(CCC_jaccard) <- 137
length(BEC_jaccard) <- 137
length(EIC_jaccard) <- 137
length(BIC_jaccard) <- 137
DCC_jaccard <- as.data.frame(DCC_jaccard)
WIC_jaccard <- as.data.frame(WIC_jaccard)
RPC_jaccard <- as.data.frame(RPC_jaccard)
NBC_jaccard <- as.data.frame(NBC_jaccard)
CBC_jaccard <- as.data.frame(CBC_jaccard)
BEC_jaccard <- as.data.frame(BEC_jaccard)
EIC_jaccard <- as.data.frame(EIC_jaccard)
BIC_jaccard <- as.data.frame(BIC_jaccard)
CCC_jaccard <- as.data.frame(CCC_jaccard)


# combine time C jaccard
jaccard_C <- melt(data.frame(BEC_jaccard, DCC_jaccard, EIC_jaccard, BIC_jaccard, WIC_jaccard, CCC_jaccard, RPC_jaccard, NBC_jaccard, CBC_jaccard))
# rename columns, reduce sitetime values, and split into site and time C
colnames(jaccard_C) <- c('sitetime', 'value') 
jaccard_C$sitetime <- as.character(jaccard_C$sitetime)
jaccard_C$sitetime <- substr(jaccard_C$sitetime,1,3)
jaccard_C <- transform(jaccard_C, site = substr(sitetime, 1, 2))
jaccard_C <- jaccard_C %>%
  select(-sitetime) %>%
  filter(!is.na(value))
jaccard_C$site <- factor(jaccard_C$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))

############### JACCARD lm across all sites in midsummer, model selection and post hoc test

jaccardlm <- lm(jaccard_C$value ~ jaccard_C$site)
jaccard0 <- lm(jaccard_C$value ~ 1)
anova(jaccardlm, jaccard0)
model.sel(jaccardlm, jaccard0)
summary(jaccardlm)
jaccardaov <- aov(value ~ site, data = jaccard_C)
jaccardtuk <- glht(jaccardaov, linfct = mcp(site = "Tukey"))
jaccardtuk.cld <- cld(jaccardtuk)
jaccardtuk.cld <- print(jaccardtuk.cld)








############### BRAY DIVERSITY

DCA_bray <- DCA[,4:37] %>%
  vegdist(method = "bray")
DCC_bray <- DCC[,4:37] %>%
  vegdist(method = "bray")
DCE_bray <- DCE[,4:37] %>%
  vegdist(method = "bray")

WIA_bray <- WIA[,4:37] %>%
  vegdist(method = "bray")
WIC_bray <- WIC[,4:37] %>%
  vegdist(method = "bray")
WIE_bray <- WIE[,4:37] %>%
  vegdist(method = "bray")

RPA_bray <- RPA[,4:37] %>%
  vegdist(method = "bray")
RPC_bray <- RPC[,4:37] %>%
  vegdist(method = "bray")
RPE_bray <- RPE[,4:37] %>%
  vegdist(method = "bray")

NBA_bray <- NBA[,4:37] %>%
  vegdist(method = "bray")
NBC_bray <- NBC[,4:37] %>%
  vegdist(method = "bray")
NBE_bray <- NBE[,4:37] %>%
  vegdist(method = "bray")

CBA_bray <- CBA[,4:37] %>%
  vegdist(method = "bray")
CBC_bray <- CBC[,4:37] %>%
  vegdist(method = "bray")
CBE_bray <- CBE[,4:37] %>%
  vegdist(method = "bray")

BEC_bray <- BEC[,4:37] %>%
  vegdist(method = "bray")
BIC_bray <- BIC[,4:37] %>%
  vegdist(method = "bray")
CCC_bray <- CCC[,4:37] %>%
  vegdist(method = "bray")
EIC_bray <- EIC[,4:37] %>%
  vegdist(method = "bray")


# take mean of bray distance

DCA_bmean <- mean(DCA_bray, na.rm = TRUE)
DCC_bmean <- mean(DCC_bray, na.rm = TRUE)
DCE_bmean <- mean(DCE_bray, na.rm = TRUE)

WIA_bmean <- mean(WIA_bray, na.rm = TRUE)
WIC_bmean <- mean(WIC_bray, na.rm = TRUE)
WIE_bmean <- mean(WIE_bray, na.rm = TRUE)

RPA_bmean <- mean(RPA_bray, na.rm = TRUE)
RPC_bmean <- mean(RPC_bray, na.rm = TRUE)
RPE_bmean <- mean(RPE_bray, na.rm = TRUE)

NBA_bmean <- mean(NBA_bray, na.rm = TRUE)
NBC_bmean <- mean(NBC_bray, na.rm = TRUE)
NBE_bmean <- mean(NBE_bray, na.rm = TRUE)

CBA_bmean <- mean(CBA_bray, na.rm = TRUE)
CBC_bmean <- mean(CBC_bray, na.rm = TRUE)
CBE_bmean <- mean(CBE_bray, na.rm = TRUE)

# renamed time code for secondaries
BIC_bmean <- mean(BIC_bray)
BEC_bmean <- mean(BEC_bray)
CCC_bmean <- mean(CCC_bray)
EIC_bmean <- mean(EIC_bray)

all_bray <- melt(data.frame(DCA_bmean, DCC_bmean, DCE_bmean, WIA_bmean, WIC_bmean, WIE_bmean, RPA_bmean, RPC_bmean, RPE_bmean, NBA_bmean, NBC_bmean, NBE_bmean, CBA_bmean, CBC_bmean, CBE_bmean, BIC_bmean, BEC_bmean, CCC_bmean, EIC_bmean))
# rename columns, reduce sitetime values, and split into site and time
colnames(all_bray) <- c('sitetime', 'value') 
all_bray$sitetime <- as.character(all_bray$sitetime)
all_bray$sitetime <- substr(all_bray$sitetime,1,3)
all_bray <- transform(all_bray, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
all_bray <- all_bray %>%
  select(-sitetime)
# reorder factors
all_bray$site <- factor(all_bray$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))
target <- c("DC", "WI", "RP", "NB", "CB")
all_bray <- all_bray %>%
  filter(site %in% target)


### MIDSUMMER ONLY BRAY

length(DCC_bray) <- 137
length(WIC_bray) <- 137  
length(RPC_bray) <- 137  
length(NBC_bray) <- 137  
length(CBC_bray) <- 137  
length(CCC_bray) <- 137
length(BEC_bray) <- 137
length(EIC_bray) <- 137
length(BIC_bray) <- 137
DCC_bray <- as.data.frame(DCC_bray)
WIC_bray <- as.data.frame(WIC_bray)
RPC_bray <- as.data.frame(RPC_bray)
NBC_bray <- as.data.frame(NBC_bray)
CBC_bray <- as.data.frame(CBC_bray)
BEC_bray <- as.data.frame(BEC_bray)
EIC_bray <- as.data.frame(EIC_bray)
BIC_bray <- as.data.frame(BIC_bray)
CCC_bray <- as.data.frame(CCC_bray)


# combine time C bray
bray_C <- melt(data.frame(BEC_bray, DCC_bray, EIC_bray, BIC_bray, WIC_bray, CCC_bray, RPC_bray, NBC_bray, CBC_bray))
# rename columns, reduce sitetime values, and split into site and time C
colnames(bray_C) <- c('sitetime', 'value') 
bray_C$sitetime <- as.character(bray_C$sitetime)
bray_C$sitetime <- substr(bray_C$sitetime,1,3)
bray_C <- transform(bray_C, site = substr(sitetime, 1, 2))
bray_C <- bray_C %>%
  select(-sitetime) %>%
  filter(!is.na(value))
bray_C$site <- factor(bray_C$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))

############### Bray lm across all sites in midsummer, model selection and post hoc test

braylm <- lm(bray_C$value ~ bray_C$site)
bray0 <- lm(bray_C$value ~ 1)
anova(braylm, bray0)
model.sel(braylm, bray0)
summary(braylm)
brayaov <- aov(value ~ site, data = bray_C)
braytuk <- glht(brayaov, linfct = mcp(site = "Tukey"))
braytuk.cld <- cld(braytuk)
braytuk.cld <- print(braytuk.cld)











list(DCC_bray)
# Beta as raw alpha/gamma

DCA_rawbeta <- ncol(DCA[,4:37])/mean(specnumber(DCA[,4:37])) - 1
DCC_rawbeta <- ncol(DCC[,4:37])/mean(specnumber(DCC[,4:37])) - 1
DCE_rawbeta <- ncol(DCE[,4:37])/mean(specnumber(DCE[,4:37])) - 1

WIA_rawbeta <- ncol(WIA[,4:37])/mean(specnumber(WIA[,4:37])) - 1
WIC_rawbeta <- ncol(WIC[,4:37])/mean(specnumber(WIC[,4:37])) - 1
WIE_rawbeta <- ncol(WIE[,4:37])/mean(specnumber(WIE[,4:37])) - 1

RPA_rawbeta <- ncol(RPA[,4:37])/mean(specnumber(RPA[,4:37])) - 1
RPC_rawbeta <- ncol(RPC[,4:37])/mean(specnumber(RPC[,4:37])) - 1
RPE_rawbeta <- ncol(RPE[,4:37])/mean(specnumber(RPE[,4:37])) - 1

NBA_rawbeta <- ncol(NBA[,4:37])/mean(specnumber(NBA[,4:37])) - 1
NBC_rawbeta <- ncol(NBC[,4:37])/mean(specnumber(NBC[,4:37])) - 1
NBE_rawbeta <- ncol(NBE[,4:37])/mean(specnumber(NBE[,4:37])) - 1

CBA_rawbeta <- ncol(CBA[,4:37])/mean(specnumber(CBA[,4:37])) - 1
CBC_rawbeta <- ncol(CBC[,4:37])/mean(specnumber(CBC[,4:37])) - 1
CBE_rawbeta <- ncol(CBE[,4:37])/mean(specnumber(CBE[,4:37])) - 1

BIC_rawbeta <- ncol(BIC[,4:37])/mean(specnumber(BIC[,4:37])) - 1

BEC_rawbeta <- ncol(BEC[,4:37])/mean(specnumber(BEC[,4:37])) - 1

EIC_rawbeta <- ncol(EIC[,4:37])/mean(specnumber(EIC[,4:37])) - 1

CCC_rawbeta <- ncol(CCC[,4:37])/mean(specnumber(CCC[,4:37])) - 1

Sites <- rep(c("DC", "WI", "RP", "NB", "CB"),3)
Times <- c(rep("May", 5), rep("June/July",5), rep("August",5))
rawbeta <- c(DCA_rawbeta, WIA_rawbeta, RPA_rawbeta, NBA_rawbeta, CBA_rawbeta, DCC_rawbeta, WIC_rawbeta, RPC_rawbeta, NBC_rawbeta, CBC_rawbeta, DCE_rawbeta, WIE_rawbeta, RPE_rawbeta, NBE_rawbeta, CBE_rawbeta)



Raw_beta <- data.frame(Sites, Times, rawbeta)
Raw_beta$rawbeta <- as.numeric(as.character(Raw_beta$rawbeta))
Raw_beta$Sites <- factor(Raw_beta$Sites, levels = c("DC","WI","RP","NB","CB"))




################ CHECKERBOARD HEATMAP PREP

checkerboard_epicomm <- epicomm_s %>%
  select(-Sample) %>%
  filter(Time.Code2 == "C")
  
checkerboard_epicomm <- checkerboard_epicomm %>%
  ungroup() %>%
  select(-Time.Code2, -Margarites.helicinus, -Nebalia.sp., -Callianax.sp.) %>%
  group_by(site) %>%
  summarise_all(sum)

checkerboard_epicomm <- checkerboard_epicomm %>%
  gather(species, abun, -site)

checkerboard_epicomm$site <- factor(checkerboard_epicomm$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))

checkerboard_epicomm$abun <- sqrt(log10(checkerboard_epicomm$abun + 1))

checkerboard_epicomm$species <- as.factor(checkerboard_epicomm$species)

# rename species for plot
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Thick.nematode."] <- "Nematode sp. 2"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Pycnogonum.sp."] <- "Pycnogonum sp."
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Pugettia.richii"] <- "Pugettia richii"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Pugettia.richii"] <- "Pugettia richii"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Pontogeneia.rostrata"] <- "Pontogeneia rostrata"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Platynereis.bicanaliculata"] <- "Platynereis bicanaliculata"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Phyllaplysia.taylori"] <- "Phyllaplysia taylori"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Photis.brevipes"] <- "Photis brevipes"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Pagurus.quaylei"] <- "Pagurus quaylei"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Mytilus.trossulus"] <- "Mytilus trossulus"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Monocorophium.achersicum"] <- "Monocorophium achersicum"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Lottia.pelta"] <- "Lottia pelta"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Lirobittium.spp."] <- "Lirobittium spp."
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Leptochelia.dubia"] <- "Leptochelia dubia"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Lacuna.spp."] <- "Lacuna spp."
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Janua.pagastecheri"] <- "Janua pagastecheri"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Idotea.resecata"] <- "Pentidotea resecata"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Haminoea.sp."] <- "Haminoea spp."
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Halacarid.mite"] <- "Halacarid mite"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Eogammarus.confervicolus"] <- "Eogammarus confervicolus"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Cyprideis.beaconensis"] <- "Ostracoda sp."
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Copepod"] <- "Harpacticoid copepod"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Cirolana.harfordi"] <- "Cirolana harfordi"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Caprella.spp."] <- "Caprella spp."
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Balanus.spp."] <- "Balanus spp."
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Aoroides.columbiae"] <- "Aoroides columbiae"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Amphithoe.spp."] <- "Amphithoe spp."
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Amph.E..dorsal.teeth."] <- "Amphipoda sp."
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Alia.carinata"] <- "Alia carinata"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Nematode"] <- "Nematode sp. 1"

# change order of checkboard species to match table 2
checkerboard_epicomm$species <- factor(checkerboard_epicomm$species, levels = c("Caprella spp.", "Aoroides columbiae", "Pentidotea resecata", "Leptochelia dubia", "Photis brevipes", "Monocorophium achersicum", "Amphipoda sp.", "Pontogeneia rostrata", "Harpacticoid copepod", "Eogammarus confervicolus", "Amphithoe spp.", "Balanus spp.", "Cirolana harfordi", "Pugettia richii", "Pandalidae", "Pagurus quaylei", "Ostracoda sp.", "Phyllaplysia taylori", "Mytilus trossulus", "Lacuna spp.", "Lottia pelta", "Haminoea spp.", "Alia carinata", "Lirobittium spp.", "Platynereis bicanaliculata", "Janua pagastecheri", "Nematode sp. 1", "Nematode sp. 2", "Pycnogonum sp.", "Halacarid mite", "Nemertea"))

checkerboard_epicomm$species <- fct_rev(checkerboard_epicomm$species)












###################################################################################
# FIGURES                                                                         #
###################################################################################


########### FIGURE 2 - Diversity metrics for 9 sites in midsummer


# Observed Richness (ANOVA and TUKEY stats are above)

Rich_midsum_plot <- ggplot(data = rich.C, aes(x=site, y = obsrich, fill = site)) +
  geom_boxplot() +
  fill_palette(viridis(9, option = "viridis")) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text=element_text(size=18)) +
  labs(x="", y = "Observed Richness") +
  annotate("text", x = 1:9, y = 18.7, label = richtuk.cld)

R2_midsum_plot <- ggplot(data = rareC, aes(x=site, y = R2, fill = site)) +
  geom_boxplot() +
  fill_palette(viridis(9, option = "viridis")) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text=element_text(size=18)) +
  labs(x="", y = "Rarefied Richness") +
  annotate("text", x = 1:9, y = 4.1, label = raretuk.cld)


# ENS (ANOVA and TUKEY stats are above)

#ens_midsum_plot <- ggplot(div.summaryC, aes(x = site, y = R2, fill = site)) + 
#  geom_boxplot() + 
#  fill_palette(viridis(9, option = "viridis")) +
#  theme_minimal() +
#  theme(axis.text.x=element_blank()) +
#  labs(x="", y="ENS")+
#  annotate("text", x = 1:9, y = 8.7, label = enstuk.cld)

# SHANNON DIVERSITY (ANOVA and TUKEY stats are above)

shannon_midsum_plot <- ggplot(div.summaryC, aes(x = site, y = Shannon, fill = site)) + 
  geom_boxplot() + 
  fill_palette(viridis(9, option = "viridis")) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text=element_text(size=18)) +
  labs(x ="", y="Shannon Diversity")  +
  annotate("text", x = 1:9, y = 2.2, label = shantuk.cld)

# ABUNDANCE (ANOVA and TUKEY stats are above)

abund_midsum_plot <- ggplot(abund_C, aes(x = site, y = log10(abundance + 1), fill = site)) +
  geom_boxplot() +
  fill_palette(viridis(9, option = "viridis")) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text=element_text(size=18)) +
  labs(x ="", y="Log Abundance")  +
  annotate("text", x = 1:9, y = 3.3, label = abuntuk.cld)

# HELLINGER DISTANCE (ANOVA and TUKEY stats are above)

#hell_midsum_plot <- ggplot(hell_C, aes(x = site, y = value, fill = site)) + 
#  geom_boxplot() + 
#  fill_palette(viridis(9, option = "viridis")) +
#  theme_minimal() +
#  theme(axis.text.x=element_blank()) +
#  labs(x ="", y="Hellinger Distance") +
#  annotate("text", x = 1:9, y = 1.1, label = helltuk.cld)

# BRAY DISTANCE (ANOVA and TUKEY stats are above)

bray_midsum_plot <- ggplot(bray_C, aes(x = site, y = value, fill = site)) + 
  geom_boxplot() + 
  fill_palette(viridis(9, option = "viridis")) +
  theme_minimal() +
  #theme(axis.text.x=element_blank()) +
  theme(axis.text=element_text(size=18)) +
  labs(x ="", y="Bray-Curtis Distance") +
  annotate("text", x = 1:9, y = 1.05, label = braytuk.cld)

# JACCARD DISTANCE (ANOVA and TUKEY stats are above)

jaccard_midsum_plot <- ggplot(jaccard_C, aes(x = site, y = value, fill = site)) + 
  geom_boxplot() + 
  fill_palette(viridis(9, option = "viridis")) +
  theme_minimal() +
  #theme(axis.text.x=element_blank()) +
  theme(axis.text=element_text(size=18)) +
  labs(x ="", y="Jaccard Distance") +
  annotate("text", x = 1:9, y = 1.1, label = jaccardtuk.cld)

checkerboard_plot <- ggplot(checkerboard_epicomm, aes(site, species, fill = abun)) +
  geom_tile(colour = "gray30", stat = "identity") +
  scale_fill_viridis(option = "B") +
  theme_minimal() +
  guides(fill = guide_colorbar(label = TRUE, ticks = FALSE, title = "abundance", fill = c("high", "low")))


# FULL FIGURE 2

Figure2 <- ggarrange(ggarrange(Rich_midsum_plot, R2_midsum_plot, shannon_midsum_plot, abund_midsum_plot, bray_midsum_plot, jaccard_midsum_plot,
                     labels = c("A", "B", "C", "D", "E", "F"),
                     ncol = 2, nrow = 3,
                     common.legend = TRUE, legend = "right"),
                     ggarrange(checkerboard_plot, 
                               labels = "G",
                               ncol = 1, nrow = 1),                        
                     ncol = 1, nrow = 2,
                     common.legend = TRUE, legend = "right",
                     heights = c(2.1,1.5))
# annotate_figure(Figure2, bottom = text_grob("Figure 2: Measures of A) observed richness, B) shannon diversity, C) effective number of species (ENS), and  \n D) Hellinger distance across nine seagrass habitats types sampled in midsummer.", size = 10))

# ggsave to increase dpi to 300 for pub. 
ggsave("Figure2hires_v2", device = "png", width = 6, height = 7, units = 'in', dpi = 500)

# best size: ~800x1100 pix





########### FIGURE 4

########### JACCARD DISTANCE THROUGH TIME

jacc_plot <- ggplot(jacc_prime, aes(x = time, y = value, group = site)) + 
  geom_point(size=3, aes(colour = site, pch = time)) +
  geom_line(aes(color = site)) +
  scale_color_viridis(discrete=TRUE) +
  theme_minimal() +
  ylim(0.35, 0.9) +
  theme(axis.text.x=element_blank()) +
  labs(x=" May        June/July        August", y="Jacc. Dist.") 


# best size: ~600x400

########### BRAY DISTANCE THROUGH TIME

bray_plot <- ggplot(all_bray, aes(x = time, y = value, group = site)) + 
  geom_point(size=3, aes(colour = site, pch = time)) +
  geom_line(aes(color = site)) +
  scale_color_viridis(discrete=TRUE) +
  theme_minimal() +
  ylim(0.2, 0.8) +
  theme(axis.text.x=element_blank()) +
  labs(x=" May         June/July        August", y="B-C Dist.") 




# best size: ~600x400

# raw beta plot
#rawbeta_plot <- ggplot(na.omit(Raw_beta), aes(x = Times, y = rawbeta, group = Sites#)) +
#  geom_point(size=4, aes(color = Sites)) + 
#  geom_line(aes(color = Sites)) +
#  scale_color_viridis(discrete = TRUE, begin = 0.3) +
#  theme_minimal() +
#  theme(axis.text.x=element_blank()) +
#  labs(x="", y="Gamma/Mean(Alpha)") 

# OBSERVED RICHNESS

rich_plot <- ggplot(rich_prim, aes(x = Time.Code2, y = obsrich, group = site)) + 
  geom_point(size=3, aes(colour = site, pch = Time.Code2)) +
  geom_line(aes(color = site)) +
  scale_color_viridis(discrete=TRUE) +
  theme_minimal() +
  ylim(3.5, 15.5) +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Obs. Rich.") 


# RAREFIED RICHNESS

rare_plot <- ggplot(new_rare_prim, aes(x = Time.Code2, y = R2, group = site)) + 
  geom_point(size=3, aes(colour = site, pch = Time.Code2)) +
  geom_line(aes(color = site)) +
  scale_color_viridis(discrete=TRUE) +
  theme_minimal() +
  ylim(1.75, 2.35) +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Rare. Rich.") 


# ENS

ens_plot <- ggplot(ens_prim, aes(x = Time.Code2, y = ENS, group = site)) + 
  geom_point(size=3, aes(colour = site, pch = Time.Code2)) +
  geom_line(aes(color = site)) +
  scale_color_viridis(discrete=TRUE) +
  theme_minimal() +
  ylim(1.25, 6.5) +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="ENS") 

# SHANNON DIVERSITY

shannon_plot <- ggplot(shan_time, aes(x = Time.Code2, y = Shannon, group = site)) + 
  geom_point(size=3, aes(colour = site, pch = Time.Code2)) +
  geom_line(aes(color = site)) +
  scale_color_viridis(discrete=TRUE) +
  theme_minimal() +
  ylim(0.35, 1.8) +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Shan. Div.")


# ABUNDANCE

abun_plot <- ggplot(abund_prim, aes(x = Time.Code2, y = log10(abundance+1), group = site)) + 
  geom_point(size=3, aes(colour = site, pch = Time.Code2)) +
  geom_line(aes(color = site)) +
  scale_color_viridis(discrete=TRUE) +
  theme_minimal() +
  ylim(1.35, 3.3) +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Abund.") 


# Hellinger distance

#hell_plot <- ggplot(hell_prim, aes(x = time, y = value, fill = site)) + 
#  geom_boxplot() + 
#  fill_palette(viridis(5, begin = 0.3)) +
#  theme_minimal() +
#  theme(axis.text.x=element_blank()) +
#  labs(x=" May                 June/July                 August", y="Hellinger Distance")

# nMDS
mds_plot <- ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, pch = month, color = site)) +
  scale_color_viridis(discrete = TRUE) +
  theme_minimal() +
  geom_polygon(data=chulls_tax, aes(x=MDS1, y=MDS2, group=month), fill=NA, color = "grey") +
  geom_point(size = 3) +
  geom_vline(xintercept= -0.15, lty = 2, color = "light blue") 
  


# FULL FIGURE 4

# Tried to extract legend separately here to resize, but it didn't work,
# so I took it out.

#legend4 <- get_legend(mds_plot + theme(legend.position="right"))

Figure4 <- ggarrange(ggarrange(rich_plot, rare_plot, shannon_plot, abun_plot, bray_plot, jacc_plot, 
                               labels = c("A", "B", "C", "D", "E", "F"), 
                               font.label = list(size = 12),
                               ncol = 2, nrow = 3,
                               widths = .75,
                               legend = FALSE),
                     ggarrange(mds_plot,
                               labels = "G",
                               font.label = list(size = 12),
                               ncol = 1, nrow = 1),                        
                     ncol = 1, nrow = 2,
                     common.legend = TRUE, legend = "right",  
                     heights = c(2,1.5)) 

                  
#annotate_figure(Figure4, bottom = text_grob("Figure 5: Measures of A) observed richness, B) shannon diversity, and C) effective number of species (ENS) across five seagrass habitats types \n sampled in May, June/July, and August", size = 10))

# ggsave to increase dpi to 300 for pub. 
ggsave("Figure4hiresv2.png", device = "png", width = 6, height = 7, units = 'in', dpi = 500)

# best size: ~800x1100




#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#



###### SCRATCH PAD

# Proportion of grazers total

length(grep("grazer", epicomm_z$function.))
# 5454
length(epicomm_z$function.)
# 10302
5454/10302
# 0.5294118

# Proportion grazers in midsummer

julgraz <- epicomm_z %>%
  + filter(Time.Code2 == "C")
propgraz <- julgraz %>%
  filter(function. == "grazer")
nograz <- propgraz %>%
  filter(abundance != 0)
finalgraz <- nograz %>%
  group_by(site, species) %>%
  summarise(unique(species))
finalgraz$count <- 1
calcgraz <- finalgraz %>%
  group_by(site) %>%
  summarise(sum(count))
# change column name
names(calcgraz)[names(calcgraz)=="sum(count)"] <- "total"
mean(calcgraz$total)


# determine proportion of grazers within size classes

filtergraz <- names(epicomm_s[-c(1:3)])

newgraz <- data.s %>%
  filter(variable %in% filtergraz)
sizegraz <- newgraz %>%
  group_by(Sieve) %>%
  summarise(sum(value))

# determine how many spp range in all sites mid summer
test <- EpiC %>%
  select(-Sample) %>%
  group_by(site, Time.Code2) %>%
  summarise_all(sum)

# replace >1 with 1 
test[test > 0] <- 1
 
range(rowSums(test[3:36]))

#Average observed quadrat scale diversity
quaddiv <- EpiC
quaddiv %>% mutate_each(funs(replace(., . > 0, 1)))
quaddiv <- quaddiv %>%
  ungroup() %>%
  select(-c(Time.Code2, Sample))
quaddiv$spp <- rowSums(quaddiv[2:35])
quaddiv <- quaddiv %>%
  select(site, spp) %>%
  group_by(site) %>%
  summarise(meanrich = mean(spp))
range(quaddiv$meanrich)

# permutation test of multivariate homogeneity of group dispersions for all midusmmer sites
dispertest <- range01(vegdist(decostand(EpiC[,4:37], "hellinger"), "euclidean"))
groups <- factor(EpiC$site)
betamod <- betadisper(dispertest, groups, type = "centroid")
permutest(betamod)

# test across all time periods to determine if composition was stable.
gammatest <- epicomm_s
primsites <- c("DC", "WI", "RP", "NB", "CB")
primegamma <- gammatest %>%
  filter(site %in% primsites)
primeplot <- primegamma %>%
  select(-Sample) %>%
  group_by(site, Time.Code2) %>%
  summarise_all(sum)
adonismod <- adonis(primeplot[,3:36] ~ Time.Code2, data = primeplot, permutations = 999, method = "bray")
adonismod

# wizard jaccard spike?
Wizard <- epicomm_s %>%
  filter(site == "WI")
Wizard <- Wizard %>%
  select(-Sample) %>%
  group_by(site, Time.Code2) %>%
  summarise_all(sum)
