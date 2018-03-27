###################################################################################
#                                                                                ##
# Basic Diversity Metrics                                                        ##
# Data are current as of 2018-03-26                                              ##
# Data source: O'Connor Lab - UBC                                                ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2018-03-26                                                        ##
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

# 2018-03-26 Cleaned up for figures and added analyses from 2018-beta analysis.R
# 2018-03-24 Added MDS analysis
# 2018-03-05 Switched back to rawcomm dataset with treatment from Mary's script
# 2018-02-23 Added evenness code, and began Figure 3 panel.
# 2018-02-21 Created script. 

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

library(tidyverse) # manipulate data
library(reshape2) # manipulate data
library(vegan) # diversity analyses
library(viridis) # color palette
library(ggpubr) # combining plots
library(lubridate) # manipulate data
library(car) # normality tests
library(lme4)

# function to scale hellinger matrix between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# full community data
#epicomm_example <- read_csv("./data/epicomm_201802.csv")

data <- read.csv("./data/rawcomm.csv")
traits <- read.csv("./data/grazertraits3.csv")
sites <- read.csv("./data/site.info.csv")

# delete, add and redefine columns ----------------------------------------
traits <- traits[,-c(3,8:10)]

data$date1 <- mdy(data$date)

# data preparation ---------------------------------------------------
## Melt and recast so the datafile goes from long to wide (species as columns)
data.m <- melt(data, id = c(1,2,3,4,5,52, 53))

## Clean and correct species names. 
## in the data file, some names we initially used for taxa needed updating based on improvements in our ability to identify them. So these replacements reflect those updates. 
levels(data.m$variable)[levels(data.m$variable)== "Bittium.spp."] <- "Lirobittium.spp."
levels(data.m$variable)[levels(data.m$variable)== "Olivella.sp."] <- "Callianax.sp."
levels(data.m$variable)[levels(data.m$variable)== "Cypricercus."] <- "Cyprideis.beaconensis"
levels(data.m$variable)[levels(data.m$variable)== "Odontosyllis"] <- "Polychaete1"
#levels(data.m$variable)[levels(data.m$variable)== "Idotea.resecata"] <- "Pentidotea.resecata"

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
## Sum across size classes within plots (samples)
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
epicomm_s <- epicomm_s %>%
  select(-Date, -dfw, -order, -area, -salinity, -shoot.density, -fetch, -TimeID, -function., -taxon, -group, -eelgrss.epifauna)

# replace NAs with 0 
epicomm_s[is.na(epicomm_s)] <- 0

# summarize by sample
epicomm_s <- epicomm_s %>%
  group_by(site, Time.Code2, Sample) %>%
  summarize_all(funs(sum))

# extract July (midsummer samples) for separate analyses of alpha diversity (did not sure rarefied richness because requires subsampling of smallest species count. Per quadrat that was sometimes zero, so instead report actually observed richness)
epicomm_midsum <- epicomm_s %>%
  select(-Sample) %>%
  filter(Time.Code2 == "C")
epicomm_midsum <- epicomm_midsum[,-2]
epicomm_midsum[epicomm_midsum > 0] <- 1 
epicomm_midsum$richness <- richness <- rowSums(epicomm_midsum[,2:35])
epicomm_richness <- epicomm_midsum[, c(1,36)]
epicomm_richness$richness <- as.numeric(epicomm_richness$richness)
#reorder factors
epicomm_richness$site <- factor(epicomm_richness$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))




# create sitetime column for analyses
epicomm <- epicomm_s %>%
  select(-Sample) %>%
  unite(sitetime, site, Time.Code2, sep = "")


# remove redundant columns and summarise species occurrences
# epicomm_summ <- epicomm %>%
#  select(-date, -Sieve) %>%
#  unite(sitetime, site, 'Time Code', sep = "", remove = TRUE)
# fix CBC typo
# epicomm_summ$sitetime <- gsub("CBC ", "CBC", epicomm_summ$sitetime)
# finish summarization for every sample by summing
# epicomm_summ <- epicomm_summ %>%
#  group_by(sitetime, Sample) %>%
#  summarise_all(funs(sum))

# group all primary sites 
primary <- c("DCA", "DCC", "DCE", "WIA", "WIC", "WIE", "RPA", "RPC", "RPE", "NBA", "NBC", "NBE", "CBA", "CBC", "CBE")
epicomm_prim_full <- epicomm %>%
  subset(sitetime %in% primary) 

# separate each sitetime for analysis

DCA <- epicomm_prim_full %>%
  subset(sitetime == "DCA")
DCC <- epicomm_prim_full %>%
  subset(sitetime == "DCC")
DCE <- epicomm_prim_full %>%
  subset(sitetime == "DCE")

WIA <- epicomm_prim_full %>%
  subset(sitetime == "WIA")
WIC <- epicomm_prim_full %>%
  subset(sitetime == "WIC")
WIE <- epicomm_prim_full %>%
  subset(sitetime == "WIE")

RPA <- epicomm_prim_full %>%
  subset(sitetime == "RPA")
RPC <- epicomm_prim_full %>%
  subset(sitetime == "RPC")
RPE <- epicomm_prim_full %>%
  subset(sitetime == "RPE")

NBA <- epicomm_prim_full %>%
  subset(sitetime == "NBA")
NBC <- epicomm_prim_full %>%
  subset(sitetime == "NBC")
NBE <- epicomm_prim_full %>%
  subset(sitetime == "NBE")

CBA <- epicomm_prim_full %>%
  subset(sitetime == "CBA")
CBC <- epicomm_prim_full %>%
  subset(sitetime == "CBC")
CBE <- epicomm_prim_full %>%
  subset(sitetime == "CBE")

# group all secondary sites 
secondary <- c("BIC", "BEC", "CCC", "EIC")
epicomm_sec_full <- epicomm %>%
  subset(sitetime %in% secondary) 

# separate each sitetime for analysis and rename time period

BIC <- epicomm_sec_full %>%
  subset(sitetime == "BIC")
BEC <- epicomm_sec_full %>%
  subset(sitetime == "BEC")
EIC <- epicomm_sec_full %>%
  subset(sitetime == "EIC")
CCC <- epicomm_sec_full %>%
  subset(sitetime == "CCC")

# separate into time periods for shannon and ens
EpiA <- bind_rows(DCA, WIA, RPA, NBA, CBA)

EpiC <- bind_rows(DCC, WIC, BEC, EIC, RPC, NBC, CBC, BIC, CCC)
EpiC <- transform(EpiC, site = substr(sitetime, 1, 2))
EpiC$site <- as.factor(EpiC$site)

EpiE <- bind_rows(DCE, WIE, RPE, NBE, CBE)

###################################################################################
# COMMUNITY DESCRIPTION                                                           #
###################################################################################

# this code primarily corrects abundance values in Table 2.

# Overall rank abundance of taxa across all sites
total_abundance <- colSums(epicomm[,2:35], na.rm = TRUE) 
total_abundance <- sort(total_abundance, decreasing = TRUE)

# DCA rank abundance
DCA_abundance <- colSums(DCA[,2:35])
DCA_abundance <- sort(DCA_abundance, decreasing = TRUE)
# DCC rank abundance
DCC_abundance <- colSums(DCC[,2:35])
DCC_abundance <- sort(DCC_abundance, decreasing = TRUE)
# DCE rank abundance
DCE_abundance <- colSums(DCE[,2:35])
DCE_abundance <- sort(DCE_abundance, decreasing = TRUE)
# WIA rank abundance
WIA_abundance <- colSums(WIA[,2:35])
WIA_abundance <- sort(WIA_abundance, decreasing = TRUE)
# WIC rank abundance
WIC_abundance <- colSums(WIC[,2:35])
WIC_abundance <- sort(WIC_abundance, decreasing = TRUE)
# WIE rank abundance
WIE_abundance <- colSums(WIE[,2:35])
WIE_abundance <- sort(WIE_abundance, decreasing = TRUE)
# RPA rank abundance
RPA_abundance <- colSums(RPA[,2:35])
RPA_abundance <- sort(RPA_abundance, decreasing = TRUE)
# RPC rank abundance
RPC_abundance <- colSums(RPC[,2:35])
RPC_abundance <- sort(RPC_abundance, decreasing = TRUE)
# RPE rank abundance
RPE_abundance <- colSums(RPE[,2:35])
RPE_abundance <- sort(RPE_abundance, decreasing = TRUE)
# NBA rank abundance
NBA_abundance <- colSums(NBA[,2:35])
NBA_abundance <- sort(NBA_abundance, decreasing = TRUE)
# NBC rank abundance
NBC_abundance <- colSums(NBC[,2:35])
NBC_abundance <- sort(NBC_abundance, decreasing = TRUE)
# NBE rank abundance
NBE_abundance <- colSums(NBE[,2:35])
NBE_abundance <- sort(NBE_abundance, decreasing = TRUE)
# CBA rank abundance
CBA_abundance <- colSums(CBA[,2:35])
CBA_abundance <- sort(CBA_abundance, decreasing = TRUE)
# CBC rank abundance
CBC_abundance <- colSums(CBC[,2:35])
CBC_abundance <- sort(CBC_abundance, decreasing = TRUE)
# CBE rank abundance
CBE_abundance <- colSums(CBE[,2:35])
CBE_abundance <- sort(CBE_abundance, decreasing = TRUE)
# BEC rank abundance
BEC_abundance <- colSums(BEC[,2:35])
BEC_abundance <- sort(BEC_abundance, decreasing = TRUE)
# EIC rank abundance
EIC_abundance <- colSums(EIC[,2:35])
EIC_abundance <- sort(EIC_abundance, decreasing = TRUE)
# BIC rank abundance
BIC_abundance <- colSums(BIC[,2:35])
BIC_abundance <- sort(BIC_abundance, decreasing = TRUE)
# CCC rank abundance
CCC_abundance <- colSums(CCC[,2:35])
CCC_abundance <- sort(CCC_abundance, decreasing = TRUE)

###################################################################################
# DIVERSITY MEASURES                                                              #
###################################################################################

############### OBSERVED RICHNESS

DCA_rich <- DCA[,2:35] %>%
  specnumber()
DCC_rich <- DCC[,2:35] %>%
  specnumber()
DCE_rich <- DCE[,2:35] %>%
  specnumber()

WIA_rich <- WIA[,2:35] %>%
  specnumber()
WIC_rich <- WIC[,2:35] %>%
  specnumber()
WIE_rich <- WIE[,2:35] %>%
  specnumber()

RPA_rich <- RPA[,2:35] %>%
  specnumber()
RPC_rich <- RPC[,2:35] %>%
  specnumber()
RPE_rich <- RPE[,2:35] %>%
  specnumber()

NBA_rich <- NBA[,2:35] %>%
  specnumber()
NBC_rich <- NBC[,2:35] %>%
  specnumber()
NBE_rich <- NBE[,2:35] %>%
  specnumber()

CBA_rich <- CBA[,2:35] %>%
  specnumber()
CBC_rich <- CBC[,2:35] %>%
  specnumber()
CBE_rich <- CBE[,2:35] %>%
  specnumber()

# combine into single data frame
# make all same length
length(DCA_rich) <- 17                      
length(DCC_rich) <- 17  
length(DCE_rich) <- 17 
length(WIA_rich) <- 17 
length(WIC_rich) <- 17  
length(WIE_rich) <- 17  
length(RPA_rich) <- 17  
length(RPC_rich) <- 17  
length(RPE_rich) <- 17  
length(NBA_rich) <- 17  
length(NBC_rich) <- 17  
length(NBE_rich) <- 17  
length(CBA_rich) <- 17  
length(CBC_rich) <- 17  
length(CBE_rich) <- 17 
# combine
rich_prim <- melt(data.frame(DCA_rich, DCC_rich, DCE_rich, WIA_rich, WIC_rich, WIE_rich, RPA_rich, RPC_rich, RPE_rich, NBA_rich, NBC_rich, NBE_rich, CBA_rich, CBC_rich, CBE_rich))
# rename columns, reduce sitetime values, and split into site and time
colnames(rich_prim) <- c('sitetime', 'value') 
rich_prim$sitetime <- as.character(rich_prim$sitetime)
rich_prim$sitetime <- substr(rich_prim$sitetime,1,3)
rich_prim <- transform(rich_prim, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
rich_prim <- rich_prim %>%
  select(-sitetime) %>%
  na.omit()
# reorder factors
rich_prim$site <- factor(rich_prim$site, levels = c("DC", "WI", "RP", "NB", "CB"))



############### Richness ANOVA across all sites in midsummer

# test for homogeneity 
qqnorm(epicomm_richness$richness)
qqline(epicomm_richness$richness)
leveneTest(richness ~ site, data = epicomm_richness)

aovrich <- aov(richness ~ site, data = epicomm_richness)
summary(aovrich)
TukeyHSD(aovrich)



############### ENS

DCA_ens <- exp(diversity(DCA[,2:35], "shannon"))
DCC_ens <- exp(diversity(DCC[,2:35], "shannon"))
DCE_ens <- exp(diversity(DCE[,2:35], "shannon"))

WIA_ens <- exp(diversity(WIA[,2:35], "shannon"))
WIC_ens <- exp(diversity(WIC[,2:35], "shannon"))
WIE_ens <- exp(diversity(WIE[,2:35], "shannon"))

RPA_ens <- exp(diversity(RPA[,2:35], "shannon"))
RPC_ens <- exp(diversity(RPC[,2:35], "shannon"))
RPE_ens <- exp(diversity(RPE[,2:35], "shannon"))

NBA_ens <- exp(diversity(NBA[,2:35], "shannon"))
NBC_ens <- exp(diversity(NBC[,2:35], "shannon"))
NBE_ens <- exp(diversity(NBE[,2:35], "shannon"))

CBA_ens <- exp(diversity(CBA[,2:35], "shannon"))
CBC_ens <- exp(diversity(CBC[,2:35], "shannon"))
CBE_ens <- exp(diversity(CBE[,2:35], "shannon"))

############## SECONDARY SITES TO JOIN SEPARATELY FOR FIGURE 2

BEC_ens <- exp(diversity(BEC[,2:35], "shannon"))
EIC_ens <- exp(diversity(EIC[,2:35], "shannon"))
BIC_ens <- exp(diversity(BIC[,2:35], "shannon"))
CCC_ens <- exp(diversity(BIC[,2:35], "shannon"))

# combine into single data frame
# make all same length
length(DCA_ens) <- 17                      
length(DCC_ens) <- 17  
length(DCE_ens) <- 17 
length(WIA_ens) <- 17 
length(WIC_ens) <- 17  
length(WIE_ens) <- 17  
length(RPA_ens) <- 17  
length(RPC_ens) <- 17  
length(RPE_ens) <- 17  
length(NBA_ens) <- 17  
length(NBC_ens) <- 17  
length(NBE_ens) <- 17  
length(CBA_ens) <- 17  
length(CBC_ens) <- 17  
length(CBE_ens) <- 17 
#secondaries
length(BEC_ens) <- 17
length(EIC_ens) <- 17
length(BIC_ens) <- 17
length(CCC_ens) <- 17
# combine primaries
ens_prim <- melt(data.frame(DCA_ens, DCC_ens, DCE_ens, WIA_ens, WIC_ens, WIE_ens, RPA_ens, RPC_ens, RPE_ens, NBA_ens, NBC_ens, NBE_ens, CBA_ens, CBC_ens, CBE_ens))
# combine July all sites
ens_midsum <-  melt(data.frame(DCC_ens, WIC_ens, BEC_ens, EIC_ens, RPC_ens, NBC_ens, CBC_ens, BIC_ens, CCC_ens))
# rename columns, reduce sitetime values, and split into site and time PRIMARY
colnames(ens_prim) <- c('sitetime', 'value') 
ens_prim$sitetime <- as.character(ens_prim$sitetime)
ens_prim$sitetime <- substr(ens_prim$sitetime,1,3)
ens_prim <- transform(ens_prim, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
ens_prim <- ens_prim %>%
  select(-sitetime) %>%
  na.omit()
ens_prim$site <- factor(ens_prim$site, levels = c("DC", "WI", "RP", "NB", "CB"))
# rename columns, reduce sitetime values, and split into site and time MIDSUMMER
colnames(ens_midsum) <- c('sitetime', 'value') 
ens_midsum$sitetime <- as.character(ens_midsum$sitetime)
ens_midsum$sitetime <- substr(ens_midsum$sitetime,1,3)
ens_midsum <- transform(ens_midsum, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
ens_midsum <- ens_midsum %>%
  select(-sitetime, -time) %>%
  na.omit()
# reorder factors
ens_midsum$site <- factor(ens_midsum$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))

# FULL TIME PERIOD SHANNON


Shannon <- diversity(EpiC[,-c(1,36,37)], index ="shannon")
# Compile indices into one dataframe to make FIGURE 2
div.summary <- cbind(EpiC[36], Shannon)
#div.summaryE <- merge(div.summary, RR.data, by.x = c("site", "Date","Sample", "alpha.p", "N"), by.y = c("site", "Date","alpha.p","N", "Sample"))
div.summary2 <- merge(div.summary, sites, by.x = c("site"), by.y = c("site"))


############### ENS ANOVA across all sites in midsummer

# test for homogeneity 
qqnorm(ens_midsum$value)
qqline(ens_midsum$value)
leveneTest(value ~ site, data = ens_midsum)

aovens <- aov(value ~ site, data = ens_midsum)
summary(aovens)
TukeyHSD(aovens)


############### SHANNON DIVERSITY

DCA_shannon <- DCA[,2:35] %>%
  diversity("shannon")
DCC_shannon <- DCC[,2:35] %>%
  diversity("shannon")
DCE_shannon <- DCE[,2:35] %>%
  diversity("shannon")

WIA_shannon <- WIA[,2:35] %>%
  diversity("shannon")
WIC_shannon <- WIC[,2:35] %>%
  diversity("shannon")
WIE_shannon <- WIE[,2:35] %>%
  diversity("shannon")

RPA_shannon <- RPA[,2:35] %>%
  diversity("shannon")
RPC_shannon <- RPC[,2:35] %>%
  diversity("shannon")
RPE_shannon <- RPE[,2:35] %>%
  diversity("shannon")

NBA_shannon <- NBA[,2:35] %>%
  diversity("shannon")
NBC_shannon <- NBC[,2:35] %>%
  diversity("shannon")
NBE_shannon <- NBE[,2:35] %>%
  diversity("shannon")

CBA_shannon <- CBA[,2:35] %>%
  diversity("shannon")
CBC_shannon <- CBC[,2:35] %>%
  diversity("shannon")
CBE_shannon <- CBE[,2:35] %>%
  diversity("shannon")

############## SECONDARY SITES TO JOIN SEPARATELY FOR FIGURE 2

BEC_shannon <- BEC[,2:35] %>%
  diversity("shannon")
EIC_shannon <- EIC[,2:35] %>%
  diversity("shannon")
BIC_shannon <- BIC[,2:35] %>%
  diversity("shannon")
CCC_shannon <- CCC[,2:35] %>%
  diversity("shannon")



# combine primary into single data frame
# make all same length
length(DCA_shannon) <- 17                      
length(DCC_shannon) <- 17  
length(DCE_shannon) <- 17 
length(WIA_shannon) <- 17 
length(WIC_shannon) <- 17  
length(WIE_shannon) <- 17  
length(RPA_shannon) <- 17  
length(RPC_shannon) <- 17  
length(RPE_shannon) <- 17  
length(NBA_shannon) <- 17  
length(NBC_shannon) <- 17  
length(NBE_shannon) <- 17  
length(CBA_shannon) <- 17  
length(CBC_shannon) <- 17  
length(CBE_shannon) <- 17 
#secondaries
length(BEC_shannon) <- 17
length(EIC_shannon) <- 17
length(BIC_shannon) <- 17
length(CCC_shannon) <- 17
# combine primaries
shannon_prim <- melt(data.frame(DCA_shannon, DCC_shannon, DCE_shannon, WIA_shannon, WIC_shannon, WIE_shannon, RPA_shannon, RPC_shannon, RPE_shannon, NBA_shannon, NBC_shannon, NBE_shannon, CBA_shannon, CBC_shannon, CBE_shannon))
# combine July all sites
shannon_midsum <-  melt(data.frame(DCC_shannon, WIC_shannon, BEC_shannon, EIC_shannon, RPC_shannon, NBC_shannon, CBC_shannon, BIC_shannon, CCC_shannon))
# rename columns, reduce sitetime values, and split into site and time PRIMARY
colnames(shannon_prim) <- c('sitetime', 'value') 
shannon_prim$sitetime <- as.character(shannon_prim$sitetime)
shannon_prim$sitetime <- substr(shannon_prim$sitetime,1,3)
shannon_prim <- transform(shannon_prim, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
shannon_prim <- shannon_prim %>%
  select(-sitetime) %>%
  na.omit()
shannon_prim$site <- factor(shannon_prim$site, levels = c("DC", "WI", "RP", "NB", "CB"))
# rename columns, reduce sitetime values, and split into site and time MIDSUMMER
colnames(shannon_midsum) <- c('sitetime', 'value') 
shannon_midsum$sitetime <- as.character(shannon_midsum$sitetime)
shannon_midsum$sitetime <- substr(shannon_midsum$sitetime,1,3)
shannon_midsum <- transform(shannon_midsum, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
shannon_midsum <- shannon_midsum %>%
  select(-sitetime, -time) %>%
  na.omit()
# reorder factors
shannon_midsum$site <- factor(shannon_midsum$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))




############### Shannon diversity ANOVA across all sites in midsummer

# test for homogeneity 
qqnorm(shannon_midsum$value)
qqline(shannon_midsum$value)
leveneTest(value ~ site, data = shannon_midsum)

aovrich <- aov(value ~ site, data = shannon_midsum)
summary(aovrich)
TukeyHSD(aovrich)


############### HELLINGER DISTANCE

DCA_hell <- range01(vegdist(decostand(DCA[,2:35], "hellinger"), "euclidean"))
DCC_hell <- range01(vegdist(decostand(DCC[,2:35], "hellinger"), "euclidean"))
DCE_hell <- range01(vegdist(decostand(DCE[,2:35], "hellinger"), "euclidean"))

WIA_hell <- range01(vegdist(decostand(WIA[,2:35], "hellinger"), "euclidean"))
WIC_hell <- range01(vegdist(decostand(WIC[,2:35], "hellinger"), "euclidean"))
WIE_hell <- range01(vegdist(decostand(WIE[,2:35], "hellinger"), "euclidean"))

RPA_hell <- range01(vegdist(decostand(RPA[,2:35], "hellinger"), "euclidean"))
RPC_hell <- range01(vegdist(decostand(RPC[,2:35], "hellinger"), "euclidean"))
RPE_hell <- range01(vegdist(decostand(RPE[,2:35], "hellinger"), "euclidean"))

NBA_hell <- range01(vegdist(decostand(NBA[,2:35], "hellinger"), "euclidean"))
NBC_hell <- range01(vegdist(decostand(NBC[,2:35], "hellinger"), "euclidean"))
NBE_hell <- range01(vegdist(decostand(NBE[,2:35], "hellinger"), "euclidean"))

CBA_hell <- range01(vegdist(decostand(CBA[,2:35], "hellinger"), "euclidean"))
CBC_hell <- range01(vegdist(decostand(CBC[,2:35], "hellinger"), "euclidean"))
CBE_hell <- range01(vegdist(decostand(CBE[,2:35], "hellinger"), "euclidean"))

############## SECONDARY SITES TO JOIN SEPARATELY FOR FIGURE 2

BEC_hell <- range01(vegdist(decostand(BEC[,2:35], "hellinger"), "euclidean"))
EIC_hell <- range01(vegdist(decostand(EIC[,2:35], "hellinger"), "euclidean"))
BIC_hell <- range01(vegdist(decostand(BIC[,2:35], "hellinger"), "euclidean"))
CCC_hell <- range01(vegdist(decostand(CCC[,2:35], "hellinger"), "euclidean"))



# combine primary into single data frame
# make all same length
length(DCA_hell) <- 17                      
length(DCC_hell) <- 17  
length(DCE_hell) <- 17 
length(WIA_hell) <- 17 
length(WIC_hell) <- 17  
length(WIE_hell) <- 17  
length(RPA_hell) <- 17  
length(RPC_hell) <- 17  
length(RPE_hell) <- 17  
length(NBA_hell) <- 17  
length(NBC_hell) <- 17  
length(NBE_hell) <- 17  
length(CBA_hell) <- 17  
length(CBC_hell) <- 17  
length(CBE_hell) <- 17 
#secondaries
length(BEC_hell) <- 17
length(EIC_hell) <- 17
length(BIC_hell) <- 17
length(CCC_hell) <- 17
# combine primaries
hell_prim <- melt(data.frame(DCA_hell, DCC_hell, DCE_hell, WIA_hell, WIC_hell, WIE_hell, RPA_hell, RPC_hell, RPE_hell, NBA_hell, NBC_hell, NBE_hell, CBA_hell, CBC_hell, CBE_hell))
# combine July all sites
hell_midsum <-  melt(data.frame(DCC_hell, WIC_hell, BEC_hell, EIC_hell, RPC_hell, NBC_hell, CBC_hell, BIC_hell, CCC_hell))
# rename columns, reduce sitetime values, and split into site and time PRIMARY
colnames(hell_prim) <- c('sitetime', 'value') 
hell_prim$sitetime <- as.character(hell_prim$sitetime)
hell_prim$sitetime <- substr(hell_prim$sitetime,1,3)
hell_prim <- transform(hell_prim, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
hell_prim <- hell_prim %>%
  select(-sitetime) %>%
  na.omit()
hell_prim$site <- factor(hell_prim$site, levels = c("DC", "WI", "RP", "NB", "CB"))
# rename columns, reduce sitetime values, and split into site and time MIDSUMMER
colnames(hell_midsum) <- c('sitetime', 'value') 
hell_midsum$sitetime <- as.character(hell_midsum$sitetime)
hell_midsum$sitetime <- substr(hell_midsum$sitetime,1,3)
hell_midsum <- transform(hell_midsum, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
hell_midsum <- hell_midsum %>%
  select(-sitetime, -time) %>%
  na.omit()
# reorder factors
hell_midsum$site <- factor(hell_midsum$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))




############### Hellinger ANOVA across all sites in midsummer

# test for homogeneity 
qqnorm(hell_midsum$value)
qqline(hell_midsum$value)
leveneTest(value ~ site, data = hell_midsum)

aovhell <- aov(value ~ site, data = hell_midsum)
summary(aovhell)
TukeyHSD(aovhell)



############### EVENNESS

DCA_even <- diversity(DCA[,2:35])/log(specnumber(DCA[,2:35]))
DCC_even <- diversity(DCC[,2:35])/log(specnumber(DCC[,2:35]))
DCE_even <- diversity(DCE[,2:35])/log(specnumber(DCE[,2:35]))

WIA_even <- diversity(WIA[,2:35])/log(specnumber(WIA[,2:35]))
WIC_even <- diversity(WIC[,2:35])/log(specnumber(WIC[,2:35]))
WIE_even <- diversity(WIE[,2:35])/log(specnumber(WIE[,2:35]))

RPA_even <- diversity(RPA[,2:35])/log(specnumber(RPA[,2:35]))
RPC_even <- diversity(RPC[,2:35])/log(specnumber(RPC[,2:35]))
RPE_even <- diversity(RPE[,2:35])/log(specnumber(RPE[,2:35]))

NBA_even <- diversity(NBA[,2:35])/log(specnumber(NBA[,2:35]))
NBC_even <- diversity(NBC[,2:35])/log(specnumber(NBC[,2:35]))
NBE_even <- diversity(NBE[,2:35])/log(specnumber(NBE[,2:35]))

CBA_even <- diversity(CBA[,2:35])/log(specnumber(CBA[,2:35]))
CBC_even <- diversity(CBC[,2:35])/log(specnumber(CBC[,2:35]))
CBE_even <- diversity(CBE[,2:35])/log(specnumber(CBE[,2:35]))

# combine into single data frame
# make all same length
length(DCA_even) <- 17                      
length(DCC_even) <- 17  
length(DCE_even) <- 17 
length(WIA_even) <- 17 
length(WIC_even) <- 17  
length(WIE_even) <- 17  
length(RPA_even) <- 17  
length(RPC_even) <- 17  
length(RPE_even) <- 17  
length(NBA_even) <- 17  
length(NBC_even) <- 17  
length(NBE_even) <- 17  
length(CBA_even) <- 17  
length(CBC_even) <- 17  
length(CBE_even) <- 17 
# combine
even_prim <- melt(data.frame(DCA_even, DCC_even, DCE_even, WIA_even, WIC_even, WIE_even, RPA_even, RPC_even, RPE_even, NBA_even, NBC_even, NBE_even, CBA_even, CBC_even, CBE_even))
# rename columns, reduce sitetime values, and split into site and time
colnames(even_prim) <- c('sitetime', 'value') 
even_prim$sitetime <- as.character(even_prim$sitetime)
even_prim$sitetime <- substr(even_prim$sitetime,1,3)
even_prim <- transform(even_prim, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
even_prim <- even_prim %>%
  filter(!is.na(value)) %>%
  select(-sitetime)
# reorder factors
even_prim$site <- factor(even_prim$site, levels = c("DC", "WI", "RP", "NB", "CB"))




############### RARIFIED RICHNESS


DCA_rare <- DCA[,2:35] %>%
  rarefy(5)

DCA_rare <- DCA[,2:35] %>%
  rarefy(5)
DCC_rare <- DCC[,2:35] %>%
  rarefy(5)
DCE_rare <- DCE[,2:35] %>%
  rarefy(5)

WIA_rare <- WIA[,2:35] %>%
  rarefy(5)
WIC_rare <- WIC[,2:35] %>%
  rarefy(5)
WIE_rare <- WIE[,2:35] %>%
  rarefy(5)

RPA_rare <- RPA[,2:35] %>%
  rarefy(5)
RPC_rare <- RPC[,2:35] %>%
  rarefy(5)
RPE_rare <- RPE[,2:35] %>%
  rarefy(5)

NBA_rare <- NBA[,2:35] %>%
  rarefy(5)
NBC_rare <- NBC[,2:35] %>%
  rarefy(5)
NBE_rare <- NBE[,2:35] %>%
  rarefy(5)

CBA_rare <- CBA[,2:35] %>%
  rarefy(5)
CBC_rare <- CBC[,2:35] %>%
  rarefy(5)
CBE_rare <- CBE[,2:35] %>%
  rarefy(5)

# combine into single data frame
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


############### MORISITAS INDEX

# Calculates the Morisita index of dispersion, standardized index values, and the 
# so called clumpedness and uniform indices. crit = two-sided p-value used to 
# calculate critical Chi-squared values. Used for bold values in Table 2.

I.DCA <- dispindmorisita(DCA[,-1], unique.rm = TRUE, na.rm = TRUE)
I.DCC <- dispindmorisita(DCC[,-1], unique.rm = TRUE, na.rm = TRUE)
I.DCE <- dispindmorisita(DCE[,-1], unique.rm = TRUE, na.rm = TRUE)
I.WIA <- dispindmorisita(WIA[,-1], unique.rm = TRUE, na.rm = TRUE)
I.WIC <- dispindmorisita(WIC[,-1], unique.rm = TRUE, na.rm = TRUE)





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

epi_mat <- epicomm_mds[,3:36]

epiMDS <- metaMDS(epi_mat)
epicomm_mds_points <- epiMDS$points
epicomm_mds_points <- data.frame(epicomm_mds_points)
plot_data_tax <- data.frame(epicomm_mds_points, epicomm_mds[,1:2])
library(plyr)
chulls_tax <- ddply(plot_data_tax, .(site), function(df) df[chull(df$MDS1, df$MDS2), ])
detach(package:plyr)

###################################################################################
# TIME SERIES BETA FIGURES                                                        #
###################################################################################

############### JACCARD DIVERSITY

DCA_jaccard <- DCA[,2:35] %>%
  vegdist(method = "jaccard")
DCC_jaccard <- DCC[,2:35] %>%
  vegdist(method = "jaccard")
DCE_jaccard <- DCE[,2:35] %>%
  vegdist(method = "jaccard")

WIA_jaccard <- WIA[,2:35] %>%
  vegdist(method = "jaccard")
WIC_jaccard <- WIC[,2:35] %>%
  vegdist(method = "jaccard")
WIE_jaccard <- WIE[,2:35] %>%
  vegdist(method = "jaccard")

RPA_jaccard <- RPA[,2:35] %>%
  vegdist(method = "jaccard")
RPC_jaccard <- RPC[,2:35] %>%
  vegdist(method = "jaccard")
RPE_jaccard <- RPE[,2:35] %>%
  vegdist(method = "jaccard")

NBA_jaccard <- NBA[,2:35] %>%
  vegdist(method = "jaccard")
NBC_jaccard <- NBC[,2:35] %>%
  vegdist(method = "jaccard")
NBE_jaccard <- NBE[,2:35] %>%
  vegdist(method = "jaccard")

CBA_jaccard <- CBA[,2:35] %>%
  vegdist(method = "jaccard")
CBC_jaccard <- CBC[,2:35] %>%
  vegdist(method = "jaccard")
CBE_jaccard <- CBE[,2:35] %>%
  vegdist(method = "jaccard")

BEC_jaccard <- BEC[,2:35] %>%
  vegdist(method = "jaccard")
BIC_jaccard <- BIC[,2:35] %>%
  vegdist(method = "jaccard")
CCC_jaccard <- CCC[,2:35] %>%
  vegdist(method = "jaccard")
EIC_jaccard <- EIC[,2:35] %>%
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

############### BRAY DIVERSITY

DCA_bray <- DCA[,2:35] %>%
  vegdist(method = "bray")
DCC_bray <- DCC[,2:35] %>%
  vegdist(method = "bray")
DCE_bray <- DCE[,2:35] %>%
  vegdist(method = "bray")

WIA_bray <- WIA[,2:35] %>%
  vegdist(method = "bray")
WIC_bray <- WIC[,2:35] %>%
  vegdist(method = "bray")
WIE_bray <- WIE[,2:35] %>%
  vegdist(method = "bray")

RPA_bray <- RPA[,2:35] %>%
  vegdist(method = "bray")
RPC_bray <- RPC[,2:35] %>%
  vegdist(method = "bray")
RPE_bray <- RPE[,2:35] %>%
  vegdist(method = "bray")

NBA_bray <- NBA[,2:35] %>%
  vegdist(method = "bray")
NBC_bray <- NBC[,2:35] %>%
  vegdist(method = "bray")
NBE_bray <- NBE[,2:35] %>%
  vegdist(method = "bray")

CBA_bray <- CBA[,2:35] %>%
  vegdist(method = "bray")
CBC_bray <- CBC[,2:35] %>%
  vegdist(method = "bray")
CBE_bray <- CBE[,2:35] %>%
  vegdist(method = "bray")

BEC_bray <- BEC[,2:35] %>%
  vegdist(method = "bray")
BIC_bray <- BIC[,2:35] %>%
  vegdist(method = "bray")
CCC_bray <- CCC[,2:35] %>%
  vegdist(method = "bray")
EIC_bray <- EIC[,2:35] %>%
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




# Beta as raw alpha/gamma

DCA_rawbeta <- ncol(DCA[,2:35])/mean(specnumber(DCA[,2:35])) - 1
DCC_rawbeta <- ncol(DCC[,2:35])/mean(specnumber(DCC[,2:35])) - 1
DCE_rawbeta <- ncol(DCE[,2:35])/mean(specnumber(DCE[,2:35])) - 1

WIA_rawbeta <- ncol(WIA[,2:35])/mean(specnumber(WIA[,2:35])) - 1
WIC_rawbeta <- ncol(WIC[,2:35])/mean(specnumber(WIC[,2:35])) - 1
WIE_rawbeta <- ncol(WIE[,2:35])/mean(specnumber(WIE[,2:35])) - 1

RPA_rawbeta <- ncol(RPA[,2:35])/mean(specnumber(RPA[,2:35])) - 1
RPC_rawbeta <- ncol(RPC[,2:35])/mean(specnumber(RPC[,2:35])) - 1
RPE_rawbeta <- ncol(RPE[,2:35])/mean(specnumber(RPE[,2:35])) - 1

NBA_rawbeta <- ncol(NBA[,2:35])/mean(specnumber(NBA[,2:35])) - 1
NBC_rawbeta <- ncol(NBC[,2:35])/mean(specnumber(NBC[,2:35])) - 1
NBE_rawbeta <- ncol(NBE[,2:35])/mean(specnumber(NBE[,2:35])) - 1

CBA_rawbeta <- ncol(CBA[,2:35])/mean(specnumber(CBA[,2:35])) - 1
CBC_rawbeta <- ncol(CBC[,2:35])/mean(specnumber(CBC[,2:35])) - 1
CBE_rawbeta <- ncol(CBE[,2:35])/mean(specnumber(CBE[,2:35])) - 1

BIC_rawbeta <- ncol(BIC[,2:35])/mean(specnumber(BIC[,2:35])) - 1

BEC_rawbeta <- ncol(BEC[,2:35])/mean(specnumber(BEC[,2:35])) - 1

EIC_rawbeta <- ncol(EIC[,2:35])/mean(specnumber(EIC[,2:35])) - 1

CCC_rawbeta <- ncol(CCC[,2:35])/mean(specnumber(CCC[,2:35])) - 1

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
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Lirobittium.spp."] <- "Lirobittium"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Leptochelia.dubia"] <- "Leptochelia dubia"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Lacuna.spp."] <- "Lacuna spp."
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Janua.pagastecheri"] <- "Janua pagastecheri"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Idotea.resecata"] <- "Pentidotea resecata"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Haminoea.sp."] <- "Haminoea spp."
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Halacarid.mite"] <- "Halacarid mite"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Eogammarus.confervicolus"] <- "Eogammarus confervicolus"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Cyprideis.beaconensis"] <- "Cyprideis beaconensis"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Copepod"] <- "Harpacticoid copepod"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Cirolana.harfordi"] <- "Cirolana harfordi"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Caprella.spp."] <- "Caprella spp."
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Balanus.spp."] <- "Balanus spp."
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Aoroides.columbiae"] <- "Aoroides columbiae"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Amphithoe.spp."] <- "Amphithoe spp."
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Amph.E..dorsal.teeth."] <- "Amphipoda sp."
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Alia.carinata"] <- "Alia carinata"
levels(checkerboard_epicomm$species)[levels(checkerboard_epicomm$species)== "Nematode"] <- "Nematode sp. 1"

###################################################################################
# FIGURES                                                                         #
###################################################################################


########### FIGURE 2 - Diversity metrics for 9 sites in midsummer




####### MARY'S CODE FOR TABLE 2





H.plot <- ggplot(data = div.summary2, aes(reorder(site, -order.dfw), y = H, ymin = 0, ymax = 2)) + 
  theme_bw() +
  geom_boxplot() +
  geom_point(x = 5, y = 1.9, pch = '*', size = 8, colour = "gray50") +
  xlab("Site") +
  ylab("Shannon Diversity")

H.plot
ggsave("Jan2017Hplot.png", device = "png", width = 4, height = 2.5)

S.plot <- ggplot(data = div.summary2, aes(reorder(site, -order.dfw), y = S, ymin = 0, ymax = 1)) + 
  theme_bw() +
  geom_boxplot() +
  geom_point(x = 5, y = 1, pch = '*', size = 8, colour = "gray50") +
  geom_point(x = 1, y = 1, pch = '*', size = 8, colour = "gray50") +
  xlab("Site") +
  ylab("Simpson Evenness")

S.plot
ggsave("Jan2017Splot.png", device = "png", width = 4, height = 2.5)


RR.plot <- ggplot(data = div.summary2, aes(reorder(site, -order.dfw), y = RR, ymin = 0, ymax = 4)) + 
  theme_bw() +
  geom_boxplot() +
  geom_point(x = 6, y = 3.5, pch = '*', size = 8, colour = "gray50") +
  geom_point(x = 4, y = 3.5, pch = '*', size = 8, colour = "gray50") +
  geom_point(x = 1, y = 3.5, pch = '*', size = 8, colour = "gray50") +
  xlab("Site") +
  ylab("RR") 

RR.plot
ggsave("Jan2017RRplot.png", device = "png", width = 4, height = 2.5)


ggsave("histogram", device = "png", width = 4, height = 4)

# factor levels
div.summary2$site <- factor(div.summary2$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))

div.plot <- ggplot(data = div.summary2, aes(x= site, y = alpha.p, ymin = 0, ymax = 12, fill = site)) +
  geom_boxplot() +
  fill_palette(viridis(9, option = "viridis")) +
  theme_minimal() +
  annotate("text", x = 4, y = 18.5, label = "*", size = 8) +
  annotate("text", x = 5, y = 18.5, label = "*") +
  annotate("text", x = 7, y = 18.5, label = "*") +
  theme(axis.text.x=element_blank()) +
  labs(x="", y = "Richness")

# Richness (ANOVA and TUKEY stats are above)

#rich_midsum_plot <- ggplot(epicomm_richness, aes(x = site, y = richness, fill = site)) +
#  geom_boxplot() +
#  fill_palette(viridis(9, option = "viridis")) +
#  theme_minimal() +
#  theme(axis.text.x=element_blank()) +
#  labs(x="", y = "Richness") +
#  annotate("text", x = 1:9, y = 11.9, label = c("ac", "ac", "ac", "a", "b", "ac", "ac", "a", "c"))



# ENS (ANOVA and TUKEY stats are above)

ens_midsum_plot <- ggplot(ens_midsum, aes(x = site, y = value, fill = site)) + 
  geom_boxplot() + 
  fill_palette(viridis(9, option = "viridis")) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="ENS") +
  annotate("text", x = 1:9, y = 6.3, label = c("a", "ab", "a", "ab", "b", "a", "a", "ab", "ab"))

# SHANNON DIVERSITY (ANOVA and TUKEY stats are above)

shannon_midsum_plot <- ggplot(shannon_midsum, aes(x = site, y = value, fill = site)) + 
  geom_boxplot() + 
  fill_palette(viridis(9, option = "viridis")) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x ="", y="Shannon Diversity") +
  annotate("text", x = 1:9, y = 1.9, label = c("ac", "abc", "ac", "ab", "b", "abc", "ac", "ab", "c"))



ggplot(div.summary, aes(x = sitetime, y = Shannon, fill = sitetime)) + 
  geom_boxplot() + 
  fill_palette(viridis(9, option = "viridis")) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x ="", y="Shannon Diversity") +
  annotate("text", x = 1:9, y = 1.9, label = c("ac", "abc", "ac", "ab", "b", "abc", "ac", "ab", "c"))

# HELLINGER DISTANCE (ANOVA and TUKEY stats are above)

hell_midsum_plot <- ggplot(hell_midsum, aes(x = site, y = value, fill = site)) + 
  geom_boxplot() + 
  fill_palette(viridis(9, option = "viridis")) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x ="", y="Hellinger Distance") +
  annotate("text", x = 1:9, y = 1.1, label = c("a", "a", "a", "b", "a", "ab", "a", "a", "a"))

checkerboard_plot <- ggplot(checkerboard_epicomm, aes(site, species, fill = abun)) +
  geom_tile(colour = "gray30", stat = "identity") +
  scale_fill_viridis(option = "B") +
  theme_minimal() +
  guides(fill = guide_colorbar(label = TRUE, ticks = FALSE, title = "abundance", fill = c("high", "low")))

# FULL FIGURE 2

Figure2 <- ggarrange(ggarrange(rich_midsum_plot, ens_midsum_plot, shannon_midsum_plot, hell_midsum_plot,
                     labels = c("A", "B", "C", "D"),
                     ncol = 2, nrow = 2,
                     common.legend = TRUE, legend = "right"),
                     ggarrange(checkerboard_plot,
                               labels = "E",
                               ncol = 1, nrow = 1,
                               legend = "right"),
                     ncol = 1, nrow = 2)
# annotate_figure(Figure2, bottom = text_grob("Figure 2: Measures of A) observed richness, B) shannon diversity, C) effective number of species (ENS), and  \n D) Hellinger distance across nine seagrass habitats types sampled in midsummer.", size = 10))

# best size: ~700x1020





########### FIGURE 4

########### JACCARD DISTANCE THROUGH TIME

jacc_plot <- ggplot(all_jaccard, aes(x = time, y = value, group = site)) + 
  geom_point(size=4, aes(colour = site)) +
  geom_line(aes(color = site)) +
  scale_color_viridis(discrete=TRUE) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Mean Jaccard Distance")

# best size: ~600x400

########### BRAY DISTANCE THROUGH TIME

bray_plot <- ggplot(all_bray, aes(x = time, y = value, group = site)) + 
  geom_point(size=4, aes(colour = site)) +
  geom_line(aes(color = site)) +
  scale_color_viridis(discrete=TRUE, begin = 0.3) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Mean Bray-Curtis Dist.")

# best size: ~600x400

# raw beta plot
rawbeta_plot <- ggplot(na.omit(Raw_beta), aes(x = Times, y = rawbeta, group = Sites)) +
  geom_point(size=4, aes(color = Sites)) + 
  geom_line(aes(color = Sites)) +
  scale_color_viridis(discrete = TRUE, begin = 0.3) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Gamma/Mean(Alpha)") 

# OBSERVED RICHNESS

rich_plot <- ggplot(rich_prim, aes(x = time, y = value, fill = site)) + 
  geom_boxplot() + 
  fill_palette(viridis(5, begin = 0.3)) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Richness")

# RAREFIED RICHNESS

#rare_plot <- ggplot(rare_prim, aes(x = time, y = value, fill = site)) + 
#  geom_boxplot() + 
#  fill_palette(viridis(5, begin = 0.3)) +
#  theme_minimal() +
#  theme(axis.text.x=element_blank()) +
#  labs(x="", y="Rarefied Richness")

# ENS

ens_plot <- ggplot(ens_prim, aes(x = time, y = value, fill = site)) + 
  geom_boxplot() + 
  fill_palette(viridis(5, begin = 0.3)) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="ENS")

# SHANNON DIVERSITY

shannon_plot <- ggplot(shannon_prim, aes(x = time, y = value, fill = site)) + 
  geom_boxplot() + 
  fill_palette(viridis(5, begin = 0.3)) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x=" May                  June/July                 August", y="Shannon Diversity")

# EVENNESS

even_plot <- ggplot(even_prim, aes(x = time, y = value, fill = site)) + 
  geom_boxplot() + 
  fill_palette(viridis(5, begin = 0.3)) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x=" May                 June/July                 August", y="Evenness")

# nMDS
mds_plot <- ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, pch = month, color = site)) +
  scale_color_viridis(discrete = TRUE, begin = 0.3) +
  theme_minimal() +
  geom_point(size = 4) + 
  geom_polygon(data=chulls_tax, aes(x=MDS1, y=MDS2, group=site), fill=NA) 


# FULL FIGURE 4


Figure4 <- ggarrange(ggarrange(rawbeta_plot, bray_plot,
                               labels = c("A", "B"),
                               ncol = 2, nrow = 1,
                               legend = FALSE), 
                     ggarrange(rich_plot, ens_plot, shannon_plot, even_plot, 
                               labels = c("C", "D", "E", "F"), 
                               ncol = 2, nrow = 2,
                               legend = FALSE),
                     ggarrange(mds_plot, 
                               labels = "G",
                               ncol = 1, nrow = 1),                        
                     ncol = 1, nrow = 3,
                     common.legend = TRUE, legend = "right",
                     heights = c(1,2,2))
                  
annotate_figure(Figure3, bottom = text_grob("Figure 4: Measures of A) observed richness, B) shannon diversity, and C) effective number of species (ENS) across five seagrass habitats types \n sampled in May, June/July, and August", size = 10))

# best size: ~800x1100




#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#


##### SCRATCH PAD


test1 <- div.summary2 %>%
  group_by(site) %>%
  summarise(mean(H))
test2 <- div.summary %>%
  group_by(site) %>%
  summarise(mean(Shannon))

div.data$site <- factor(epicomm_richness$site, levels = c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"))
div.data <- div.data %>%
  group_by(site)
levels(div.data$site)
levels(EpiC$sitetime)
EpiC$sitetime <- as.factor(EpiC$sitetime)
EpiC <- EpiC %>%
  group_by(site)



Captest <- epicomm_z %>%
  filter(TimeID == "DC.C")



filtertest <- data.tr %>%
  group_by(species, eelgrss.epifauna) %>%
  filter(species == "Caprella.spp.")


data <- read.csv("./data/rawcomm.csv")
traits <- read.csv("./data/grazertraits3.csv")
sites <- read.csv("./data/site.info.csv")

# delete, add and redefine columns ----------------------------------------
traits <- traits[,-c(3,8:10)]

data$date1 <- mdy(data$date)

# data preparation ---------------------------------------------------
## Melt and recast so the datafile goes from long to wide (species as columns)
data.m <- melt(data, id = c(1,2,3,4,5,52, 53))

## Clean and correct species names. 
## in the data file, some names we initially used for taxa needed updating based on improvements in our ability to identify them. So these replacements reflect those updates. 
levels(data.m$variable)[levels(data.m$variable)== "Bittium.spp."] <- "Lirobittium.spp."
levels(data.m$variable)[levels(data.m$variable)== "Olivella.sp."] <- "Callianax.sp."
levels(data.m$variable)[levels(data.m$variable)== "Cypricercus."] <- "Cyprideis.beaconensis"
levels(data.m$variable)[levels(data.m$variable)== "Odontosyllis"] <- "Polychaete1"
#levels(data.m$variable)[levels(data.m$variable)== "Idotea.resecata"] <- "Pentidotea.resecata"

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
## Sum across size classes within plots (samples)
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
epicomm_z <- data.tr %>% filter(eelgrss.epifauna == c("yes", "sometimes"))
epicomm_s <- epicomm_z %>%
  spread(species, abundance)

# remove unused columns
epicomm_s <- epicomm_s %>%
  select(-Date, -dfw, -order, -area, -salinity, -shoot.density, -fetch, -TimeID, -function., -taxon, -group, -eelgrss.epifauna)

# replace NAs with 0 
epicomm_s[is.na(epicomm_s)] <- 0
