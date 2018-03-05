###################################################################################
#                                                                                ##
# Beta Analysis                                                                  ##
# Data are current as of 2018-02-20                                              ##
# Data source: O'Connor Lab, UBC                                                 ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2018-02-20                                                        ##
#                                                                                ##
###################################################################################

# SUMMARY:


# Required Files (check that script is loading latest version):
# epicomm_201802.csv

# Associated Scripts:
# raup_crick.R

# TO DO

# need to remove empty rows from individual sitetime DF for distance matrices

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# RECENT CHANGES TO SCRIPT                                                        #
# LOAD PACKAGES                                                                   #
# READ IN DATA                                                                    #
# PREPARE DATA FOR ANALYSES                                                       #
# RAUP-CRICK METRICS                                                              #
# TIME SERIES BETA FIGURES                                                        #
#                                                                                 #
###################################################################################

###################################################################################
# RECENT CHANGES TO SCRIPT                                                        #
###################################################################################

# 2018-03-05 using rawcomm data with treatment from Mary
# 2018-02-20 Script tidied from previous analysis and new community (epicomm_201802.csv) used instead of old community (rawcomm.csv)

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

library(plyr)
library(tidyverse)
library(reshape2)
library(lme4)
library(vegan)
library(lubridate) # data manipulation

# YOU MUST RUN RAUP_CRICK.R FOR FUNCTIONALITY 

###################################################################################
# READ IN DATA                                                                    #
###################################################################################

# full community data
#epicomm <- read.csv("./data/epicomm_201802.csv")

# remove redundant columns and summarise species occurrences
#epicomm_summ <- epicomm %>%
 # select(-date, -Sieve) %>%
#  unite(sitetime, site, Time.Code, sep = "", remove = TRUE)
# fix CBC typo
#epicomm_summ$sitetime <- gsub("CBC ", "CBC", epicomm_summ$sitetime)
# finish summarization for every sample by summing
#epicomm_summ <- epicomm_summ %>%
#  group_by(sitetime, Sample) %>%
#  summarise_all(funs(sum))
# summarise entire sites with total spp abund per site
#epicomm_site <- epicomm_summ %>%
#  select(-Sample) %>%
#  group_by(sitetime) %>%
#  summarise_all(funs(sum))


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
levels(data.m$variable)[levels(data.m$variable)== "Idotea.resecata"] <- "Pentidotea.resecata"

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

## Sum across size classes within plots (samples)
data.p <- ddply(data.s, .(site, date1, Sample, Time.Code2, variable, dfw,order.dfw,area,salinity, shoot.density, fetch.jc), summarise, sum(value))

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

# summarize by sample
epicomm_s <- epicomm_s %>%
  group_by(site, Time.Code2, Sample) %>%
  summarize_all(funs(sum))

# separate out middle time for Raup-Crick test
epicomm_C <- epicomm_s %>%
  filter(Time.Code2 == "C") %>%
  select(-Sample) %>%
  unite(sitetime, site, Time.Code2, sep = "")

# create sitetime column for analyses
epicomm <- epicomm_s %>%
  select(-Sample) %>%
  unite(sitetime, site, Time.Code2, sep = "") 

# pull out full unsummed sites for site by site analysis
epicomm_full <- epicomm_s %>%
  unite(sitetime, site, Time.Code2, sep = "") %>%
  select(-Sample)
  
# sum all sites
epicomm <- epicomm %>%
  group_by(sitetime) %>%
  summarize_all(funs(sum))

###################################################################################
# PREPARE DATA FOR ANALYSES                                                       #
###################################################################################

beta_sim_data <- epicomm

###################################################################################
# RAUP-CRICK METRICS                                                              #
###################################################################################

# Divided into three time periods, all sites

# Time A 
beta_A <- beta_sim_data[c(3, 7, 11, 14, 17),]
# Time A raup_crick analysis
T_rc_1 <- raup_crick(beta_A, plot_names_in_col1 = TRUE)
T_mat_1 <- as.matrix(T_rc_1)
t2 <- melt(T_mat_1)[melt(upper.tri(T_mat_1))$value,]
names(t2) <- c("c1", "c2", "distance")
t2$time <- rep("May")
t2$pair <- seq(from = 1, to=10, by =1)

#Time B, C, D
beta_C <- beta_sim_data[c(1,2,4,6,8,10,12,15,18),]
T_rc_2 <- raup_crick(beta_C, plot_names_in_col1 = TRUE)
T_mat_2 <- as.matrix(T_rc_2)
t3 <- melt(T_mat_2)[melt(upper.tri(T_mat_2))$value,]
names(t3) <- c("c1", "c2", "distance")
t3$time <- rep("June/July")
t3$pair <- seq(from = 1, to=36, by =1)

# Time E
beta_E <- beta_sim_data[c(5,9,13,16,19),]
T_rc_3 <- raup_crick(beta_E, plot_names_in_col1 = TRUE)
T_mat_3 <- as.matrix(T_rc_3)
t4 <- melt(T_mat_3)[melt(upper.tri(T_mat_3))$value,]
names(t4) <- c("c1", "c2", "distance")
t4$time <- rep("August")
t4$pair <- seq(from = 1, to=10, by =1)

#t1_summary <- as.table(summary(T_rc_1))
#t2_summary <- as.table(summary(T_rc_2))
#t3_summary <- as.table(summary(T_rc_3))


############### PLOT ACROSS THREE TIME PERIODS

#join data

t1 <- bind_rows(t2, t3, t4)

# Summary of raup-crick values between sites
inter_raup <- t1 %>%
  group_by(time) %>%
  summarise(mean(distance))

# raup-crick variation
inter_raup_var <- t1 %>%
  group_by(time) %>%
  summarise(sd(distance))

#m1$site <- factor(m1$site, levels = c("DC", "WI", "RP", "NB", "CB"))

# Cross site comparison by time
ggplot(t1, aes(time, distance)) + 
  geom_violin(trim = TRUE, fill = "#9ecae1") +
  geom_boxplot(width = 0.1) +
  #facet_grid(site~.) +
  scale_y_continuous(breaks=seq(-1,1, 1)) +
  geom_hline(yintercept=0, color = "red", linetype = "dashed") +
  ylab("Raup-Crick Distance") +
  xlab("Time Period") +
  theme(legend.position = "none")


############### SUBSECTION HERE


###################################################################################
# TIME SERIES BETA FIGURES                                                        #
###################################################################################

# separate each sitetime for analysis

DCA <- epicomm_full %>%
  subset(sitetime == "DCA")
DCC <- epicomm_full %>%
  subset(sitetime == "DCC")
DCE <- epicomm_full %>%
  subset(sitetime == "DCE")

WIA <- epicomm_full %>%
  subset(sitetime == "WIA")
WIC <- epicomm_full %>%
  subset(sitetime == "WIC")
WIE <- epicomm_full %>%
  subset(sitetime == "WIE")

RPA <- epicomm_full %>%
  subset(sitetime == "RPA")
RPC <- epicomm_full %>%
  subset(sitetime == "RPC")
RPE <- epicomm_full %>%
  subset(sitetime == "RPE")

NBA <- epicomm_full %>%
  subset(sitetime == "NBA")
NBC <- epicomm_full %>%
  subset(sitetime == "NBC")
NBE <- epicomm_full %>%
  subset(sitetime == "NBE")

CBA <- epicomm_full %>%
  subset(sitetime == "CBA")
CBC <- epicomm_full %>%
  subset(sitetime == "CBC")
CBE <- epicomm_full %>%
  subset(sitetime == "CBE")

BEB <- epicomm_full %>%
  subset(sitetime == "BEB")
BIB <- epicomm_full %>%
  subset(sitetime == "BIB")
CCD <- epicomm_full %>%
  subset(sitetime == "CCD")
EID <- epicomm_full %>%
  subset(sitetime == "EID")

############### JACCARD DIVERSITY

DCA_jaccard <- DCA[,2:34] %>%
  vegdist(method = "jaccard")
DCC_jaccard <- DCC[,2:34] %>%
  vegdist(method = "jaccard")
DCE_jaccard <- DCE[,2:34] %>%
  vegdist(method = "jaccard")

WIA_jaccard <- WIA[,2:34] %>%
  vegdist(method = "jaccard")
WIC_jaccard <- WIC[,2:34] %>%
  vegdist(method = "jaccard")
WIE_jaccard <- WIE[,2:34] %>%
  vegdist(method = "jaccard")

RPA_jaccard <- RPA[,2:34] %>%
  vegdist(method = "jaccard")
RPC_jaccard <- RPC[,2:34] %>%
  vegdist(method = "jaccard")
RPE_jaccard <- RPE[,2:34] %>%
  vegdist(method = "jaccard")

NBA_jaccard <- NBA[,2:34] %>%
  vegdist(method = "jaccard")
NBC_jaccard <- NBC[,2:34] %>%
  vegdist(method = "jaccard")
NBE_jaccard <- NBE[,2:34] %>%
  vegdist(method = "jaccard")

CBA_jaccard <- CBA[,2:34] %>%
  vegdist(method = "jaccard")
CBC_jaccard <- CBC[,2:34] %>%
  vegdist(method = "jaccard")
CBE_jaccard <- CBE[,2:34] %>%
  vegdist(method = "jaccard")

BEB_jaccard <- BEB[,2:34] %>%
  vegdist(method = "jaccard")
BIB_jaccard <- BIB[,2:34] %>%
  vegdist(method = "jaccard")
CCD_jaccard <- CCD[,2:34] %>%
  vegdist(method = "jaccard")
EID_jaccard <- EID[,2:34] %>%
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
BIC_jmean <- mean(BIB_jaccard)
BEC_jmean <- mean(BEB_jaccard)
CCC_jmean <- mean(CCD_jaccard)
EIC_jmean <- mean(EID_jaccard)

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

DCA_bray <- DCA[,2:34] %>%
  vegdist(method = "bray")
DCC_bray <- DCC[,2:34] %>%
  vegdist(method = "bray")
DCE_bray <- DCE[,2:34] %>%
  vegdist(method = "bray")

WIA_bray <- WIA[,2:34] %>%
  vegdist(method = "bray")
WIC_bray <- WIC[,2:34] %>%
  vegdist(method = "bray")
WIE_bray <- WIE[,2:34] %>%
  vegdist(method = "bray")

RPA_bray <- RPA[,2:34] %>%
  vegdist(method = "bray")
RPC_bray <- RPC[,2:34] %>%
  vegdist(method = "bray")
RPE_bray <- RPE[,2:34] %>%
  vegdist(method = "bray")

NBA_bray <- NBA[,2:34] %>%
  vegdist(method = "bray")
NBC_bray <- NBC[,2:34] %>%
  vegdist(method = "bray")
NBE_bray <- NBE[,2:34] %>%
  vegdist(method = "bray")

CBA_bray <- CBA[,2:34] %>%
  vegdist(method = "bray")
CBC_bray <- CBC[,2:34] %>%
  vegdist(method = "bray")
CBE_bray <- CBE[,2:34] %>%
  vegdist(method = "bray")

BEB_bray <- BEB[,2:34] %>%
  vegdist(method = "bray")
BIB_bray <- BIB[,2:34] %>%
  vegdist(method = "bray")
CCD_bray <- CCD[,2:34] %>%
  vegdist(method = "bray")
EID_bray <- EID[,2:34] %>%
  vegdist(method = "bray")


# take mean of bray distance

DCA_bmean <- mean(DCA_bray)
DCC_bmean <- mean(DCC_bray)
DCE_bmean <- mean(DCE_bray)

WIA_bmean <- mean(WIA_bray)
WIC_bmean <- mean(WIC_bray)
WIE_bmean <- mean(WIE_bray)

RPA_bmean <- mean(RPA_bray)
RPC_bmean <- mean(RPC_bray)
RPE_bmean <- mean(RPE_bray)

NBA_bmean <- mean(NBA_bray)
NBC_bmean <- mean(NBC_bray)
NBE_bmean <- mean(NBE_bray)

CBA_bmean <- mean(CBA_bray)
CBC_bmean <- mean(CBC_bray)
CBE_bmean <- mean(CBE_bray)

# renamed time code for secondaries
BIC_bmean <- mean(BIB_bray)
BEC_bmean <- mean(BEB_bray)
CCC_bmean <- mean(CCD_bray)
EIC_bmean <- mean(EID_bray)

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

###################################################################################
# FIGURES                                                                         #
###################################################################################

########### JACCARD DISTANCE THROUGH TIME (fig 3+)

jacc_plot <- ggplot(all_jaccard, aes(x = time, y = value, group = site)) + 
  geom_point(size=4, aes(colour = site)) +
  geom_line(aes(color = site)) +
  scale_color_viridis(discrete=TRUE) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Mean Jaccard Distance")

# best size: ~600x400

########### BRAY DISTANCE THROUGH TIME (fig 3+)

bray_plot <- ggplot(all_bray, aes(x = time, y = value, group = site)) + 
  geom_point(size=4, aes(colour = site)) +
  geom_line(aes(color = site)) +
  scale_color_viridis(discrete=TRUE) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Mean Bray-Curtis Distance")

# best size: ~600x400


######## TEST OF RAUP-CRICK MIDDLE TIME ALL SAMPLES SEPARATELY
ALL_rc_C <- raup_crick(epicomm_C, plot_names_in_col1 = TRUE)





# Single sites through time (5 core sites)

# DCA - shared absences (makes little difference)
DCA_rc  <- epicomm_full[96:111,-1]
DCA_rc <- data.frame(DCA_rc)
DCA_rc_matrix <- raup_crick(DCA_rc, plot_names_in_col1 = FALSE)
#beanplot(DCA_rc_matrix)
D_mat <- as.matrix(DCA_rc_matrix)
d2 <- melt(D_mat)[melt(upper.tri(D_mat))$value,]
names(d2) <- c("c1", "c2", "distance")
d2$time <- rep("May")
d2$pair <- seq(from = 1, to=120, by =1)

# DCC 
DCC_rc <- epicomm_full[112:127, -1]
DCC_rc <- data.frame(DCC_rc)
DCC_rc_matrix <- raup_crick(DCC_rc, plot_names_in_col1 = FALSE)
#beanplot(DCC_rc_matrix)
D_mat <- as.matrix(DCC_rc_matrix)
d3 <- melt(D_mat)[melt(upper.tri(D_mat))$value,]
names(d3) <- c("c1", "c2", "distance")
d3$time <- rep("June/July")
d3$pair <- seq(from = 1, to=120, by =1)

# DCE
DCE_rc <- epicomm_full[128:143,-1]
DCE_rc <- data.frame(DCE_rc)
DCE_rc_matrix <- raup_crick(DCE_rc, plot_names_in_col1 = FALSE)
#beanplot(DCE_rc_matrix)
D_mat <- as.matrix(DCE_rc_matrix)
d4 <- melt(D_mat)[melt(upper.tri(D_mat))$value,]
names(d4) <- c("c1", "c2", "distance")
d4$time <- rep("August")
d4$pair <- seq(from = 1, to=120, by =1)

d5 <- bind_rows(d2, d3)
d5 <- bind_rows(d5, d4)

#ggplot(d5, aes(time, distance)) + 
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1) 

#DCA_summary <- as.table(summary(DCA_rc_matrix))
#DCC_summary <- as.table(summary(DCC_rc_matrix))
#DCE_summary <- as.table(summary(DCE_rc_matrix))

#

# WIA - shared absences (makes little difference)
WIA_rc  <- epicomm_full[255:270,-1]
WIA_rc <- data.frame(WIA_rc)
WIA_rc_matrix <- raup_crick(WIA_rc, plot_names_in_col1 = FALSE)
#beanplot(WIA_rc_matrix)
W_mat <- as.matrix(WIA_rc_matrix)
w2 <- melt(W_mat)[melt(upper.tri(W_mat))$value,]
names(w2) <- c("c1", "c2", "distance")
w2$time <- rep("May")
w2$pair <- seq(from = 1, to=120, by =1)

# WIC 
WIC_rc <- epicomm_full[271:287, -1]
WIC_rc <- data.frame(WIC_rc)
WIC_rc_matrix <- raup_crick(WIC_rc, plot_names_in_col1 = FALSE)
#beanplot(WIC_rc_matrix)
W_mat <- as.matrix(WIC_rc_matrix)
w3 <- melt(W_mat)[melt(upper.tri(W_mat))$value,]
names(w3) <- c("c1", "c2", "distance")
w3$time <- rep("June/July")
w3$pair <- seq(from = 1, to=136, by =1)

# WIE
WIE_rc <- epicomm_full[288:303,-1]
WIE_rc <- data.frame(WIE_rc)
WIE_rc_matrix <- raup_crick(WIE_rc, plot_names_in_col1 = FALSE)
#beanplot(WIE_rc_matrix)
W_mat <- as.matrix(WIE_rc_matrix)
w4 <- melt(W_mat)[melt(upper.tri(W_mat))$value,]
names(w4) <- c("c1", "c2", "distance")
w4$time <- rep("August")
w4$pair <- seq(from = 1, to=120, by =1)

w5 <- bind_rows(w2, w3)
w5 <- bind_rows(w5, w4)

#ggplot(w5, aes(time, distance)) + 
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1) 

#WIA_summary <- as.table(summary(WIA_rc_matrix))
#WIC_summary <- as.table(summary(WIC_rc_matrix))
#WIE_summary <- as.table(summary(WIE_rc_matrix))

#

# RPA - shared absences (makes little difference)
RPA_rc  <- epicomm_full[208:222,-1]
RPA_rc <- data.frame(RPA_rc)
RPA_rc_matrix <- raup_crick(RPA_rc, plot_names_in_col1 = FALSE)
#beanplot(RPA_rc_matrix)
R_mat <- as.matrix(RPA_rc_matrix)
r2 <- melt(R_mat)[melt(upper.tri(R_mat))$value,]
names(r2) <- c("c1", "c2", "distance")
r2$time <- rep("May")
r2$pair <- seq(from = 1, to=105, by =1)

# RPC 
RPC_rc <- epicomm_full[223:238, -1]
RPC_rc <- data.frame(RPC_rc)
RPC_rc_matrix <- raup_crick(RPC_rc, plot_names_in_col1 = FALSE)
#beanplot(RPC_rc_matrix)
R_mat <- as.matrix(RPC_rc_matrix)
r3 <- melt(R_mat)[melt(upper.tri(R_mat))$value,]
names(r3) <- c("c1", "c2", "distance")
r3$time <- rep("June/July")
r3$pair <- seq(from = 1, to=120, by =1)

# RPE
RPE_rc <- epicomm_full[239:254,-1]
RPE_rc <- data.frame(RPE_rc)
RPE_rc_matrix <- raup_crick(RPE_rc, plot_names_in_col1 = FALSE)
#beanplot(RPE_rc_matrix)
R_mat <- as.matrix(RPE_rc_matrix)
r4 <- melt(R_mat)[melt(upper.tri(R_mat))$value,]
names(r4) <- c("c1", "c2", "distance")
r4$time <- rep("August")
r4$pair <- seq(from = 1, to=120, by =1)

r5 <- bind_rows(r2, r3)
r5 <- bind_rows(r5, r4)

#ggplot(r5, aes(time, distance)) + 
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1) 

#RPA_summary <- as.table(summary(RPA_rc_matrix))
#RPC_summary <- as.table(summary(RPC_rc_matrix))
#RPE_summary <- as.table(summary(RPE_rc_matrix))

#

# NBA - shared absences (makes little difference)
NBA_rc  <- epicomm_full[160:175,-1]
NBA_rc <- data.frame(NBA_rc)
NBA_rc_matrix <- raup_crick(NBA_rc, plot_names_in_col1 = FALSE)
#beanplot(NBA_rc_matrix)
N_mat <- as.matrix(NBA_rc_matrix)
n2 <- melt(N_mat)[melt(upper.tri(N_mat))$value,]
names(n2) <- c("c1", "c2", "distance")
n2$time <- rep("May")
n2$pair <- seq(from = 1, to=120, by =1)

# NBC 
NBC_rc <- epicomm_full[176:191, -1]
NBC_rc <- data.frame(NBC_rc)
NBC_rc_matrix <- raup_crick(NBC_rc, plot_names_in_col1 = FALSE)
#beanplot(NBC_rc_matrix)
N_mat <- as.matrix(NBC_rc_matrix)
n3 <- melt(N_mat)[melt(upper.tri(N_mat))$value,]
names(n3) <- c("c1", "c2", "distance")
n3$time <- rep("June/July")
n3$pair <- seq(from = 1, to=120, by =1)

# NBE
NBE_rc <- epicomm_full[192:207,-1]
NBE_rc <- data.frame(NBE_rc)
NBE_rc_matrix <- raup_crick(NBE_rc, plot_names_in_col1 = FALSE)
#beanplot(NBE_rc_matrix)
N_mat <- as.matrix(NBE_rc_matrix)
n4 <- melt(N_mat)[melt(upper.tri(N_mat))$value,]
names(n4) <- c("c1", "c2", "distance")
n4$time <- rep("August")
n4$pair <- seq(from = 1, to=120, by =1)

n5 <- bind_rows(n2, n3)
n5 <- bind_rows(n5, n4)

#ggplot(n5, aes(time, distance)) + 
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1) 

#NBA_summary <- as.table(summary(NBA_rc_matrix))
#NBC_summary <- as.table(summary(NBC_rc_matrix))
#NBE_summary <- as.table(summary(NBE_rc_matrix))

#

# CBA - shared absences (makes little difference)
CBA_rc  <- epicomm_full[33:48,-1]
CBA_rc <- data.frame(CBA_rc)
CBA_rc_matrix <- raup_crick(CBA_rc, plot_names_in_col1 = FALSE)
#beanplot(CBA_rc_matrix)
C_mat <- as.matrix(CBA_rc_matrix)
c2 <- melt(C_mat)[melt(upper.tri(C_mat))$value,]
names(c2) <- c("c1", "c2", "distance")
c2$time <- rep("May")
c2$pair <- seq(from = 1, to=120, by =1)

# CBC 
CBC_rc <- epicomm_full[49:64, -1]
CBC_rc <- data.frame(CBC_rc)
CBC_rc_matrix <- raup_crick(CBC_rc, plot_names_in_col1 = FALSE)
#beanplot(CBC_rc_matrix)
C_mat <- as.matrix(CBC_rc_matrix)
c3 <- melt(C_mat)[melt(upper.tri(C_mat))$value,]
names(c3) <- c("c1", "c2", "distance")
c3$time <- rep("June/July")
c3$pair <- seq(from = 1, to=120, by =1)

# CBE
CBE_rc <- epicomm_full[65:79,-1]
CBE_rc <- data.frame(CBE_rc)
CBE_rc_matrix <- raup_crick(CBE_rc, plot_names_in_col1 = FALSE)
#beanplot(CBE_rc_matrix)
C_mat <- as.matrix(CBE_rc_matrix)
c4 <- melt(C_mat)[melt(upper.tri(C_mat))$value,]
names(c4) <- c("c1", "c2", "distance")
c4$time <- rep("August")
c4$pair <- seq(from = 1, to=105, by =1)

c5 <- bind_rows(c2, c3)
c5 <- bind_rows(c5, c4)

#ggplot(c5, aes(time, distance)) + 
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1) 

#CBA_summary <- as.table(summary(CBA_rc_matrix))
#CBC_summary <- as.table(summary(CBC_rc_matrix))
#CBE_summary <- as.table(summary(CBE_rc_matrix))

####PLOT ALL CORE SITES IN SINGLE GRAPH

#join data

d5$site <- rep("DC")
w5$site <- rep("WI")
r5$site <- rep("RP")
n5$site <- rep("NB")
c5$site <- rep("CB")

m1 <- bind_rows(d5, w5)
m1 <- bind_rows(m1, r5)
m1 <- bind_rows(m1, n5)
m1 <- bind_rows(m1, c5)

m1$site <- factor(m1$site, levels = c("DC", "WI", "RP", "NB", "CB"))

#Plot core sites only
#ggplot(m1, aes(time, distance)) + 
#  geom_violin(trim = TRUE, aes(fill = factor(site))) +
#  scale_fill_brewer(palette = "Blues", direction = -1) +
#  geom_boxplot(width = 0.1) +
#  facet_grid(site~.) +
#  scale_y_continuous(breaks=seq(-1,1, 1)) +
#  geom_hline(yintercept=0, color = "red", linetype = "dashed") +
#  ylab("Raup-Crick Distance") +
#  xlab("Time Period")

# Secondary Sites

# BEB - shared absences (makes little difference)
BEB_rc  <- epicomm_full[1:16,-1]
BEB_rc <- data.frame(BEB_rc)
BEB_rc_matrix <- raup_crick(BEB_rc, plot_names_in_col1 = FALSE)
B_mat <- as.matrix(BEB_rc_matrix)
b2 <- melt(B_mat)[melt(upper.tri(B_mat))$value,]
names(b2) <- c("c1", "c2", "distance")
b2$time <- rep("June/July")
b2$pair <- seq(from = 1, to=120, by =1)

#BEB_summary <- as.table(summary(BEB_rc_matrix))

# BIB - shared absences (makes little difference)
BIB_rc  <- epicomm_full[17:32,-1]
BIB_rc <- data.frame(BIB_rc)
BIB_rc_matrix <- raup_crick(BIB_rc, plot_names_in_col1 = FALSE)
B_mat <- as.matrix(BIB_rc_matrix)
b3 <- melt(B_mat)[melt(upper.tri(B_mat))$value,]
names(b3) <- c("c1", "c2", "distance")
b3$time <- rep("June/July")
b3$pair <- seq(from = 1, to=120, by =1)

#BIB_summary <- as.table(summary(BIB_rc_matrix))

# CCD - shared absences (makes little difference)
CCD_rc  <- epicomm_full[80:95,-1]
CCD_rc <- data.frame(CCD_rc)
CCD_rc_matrix <- raup_crick(CCD_rc, plot_names_in_col1 = FALSE)
B_mat <- as.matrix(CCD_rc_matrix)
b4 <- melt(B_mat)[melt(upper.tri(B_mat))$value,]
names(b4) <- c("c1", "c2", "distance")
b4$time <- rep("June/July")
b4$pair <- seq(from = 1, to=120, by =1)

#CCD_summary <- as.table(summary(CCD_rc_matrix))

# EID - shared absences (makes little difference)
EID_rc  <- epicomm_full[144:159,-1]
EID_rc <- data.frame(EID_rc)
EID_rc_matrix <- raup_crick(EID_rc, plot_names_in_col1 = FALSE)
B_mat <- as.matrix(EID_rc_matrix)
b5 <- melt(B_mat)[melt(upper.tri(B_mat))$value,]
names(b5) <- c("c1", "c2", "distance")
b5$time <- rep("June/July")
b5$pair <- seq(from = 1, to=120, by =1)

#EID_summary <- as.table(summary(EID_rc_matrix))

# Export all summary values of analysis
#RC_all <- rbind(DCA_summary, DCC_summary, DCE_summary, WIA_summary, WIC_summary, WIE_summary, RPA_summary, RPC_summary, RPE_summary, NBA_summary, NBC_summary, NBE_summary, CBA_summary, CBC_summary, CBE_summary, BEB_summary, BIB_summary, EID_summary, CCD_summary)
#write.csv(RC_all, file = "RC_all.csv")

#RC_time_all <- rbind(t1_summary, t2_summary, t3_summary)

#join data

b2$site <- rep("BE")
b3$site <- rep("BI")
b4$site <- rep("CC")
b5$site <- rep("EI")



b1 <- bind_rows(b2,b3)
b1 <- bind_rows(b1,b4)
b1 <- bind_rows(b1,b5)


#b1$site <- factor(b1$site, levels = c("BE", "EI", "BI", "CC"))

#ggplot(b1, aes(site, distance)) + 
#  geom_violin(trim = TRUE, aes(fill = factor(site))) +
# scale_fill_brewer(palette = "Blues", direction = -1) +
#  geom_boxplot(width = 0.1) +
#  #facet_grid(site~.) +
#  scale_y_continuous(breaks=seq(-1,1, 1)) +
#  geom_hline(yintercept=0, color = "red", linetype = "dashed") +
#  ylab("Raup-Crick Distance") +
#  xlab("Site")

# Join primary with secondary sites

x1 <- bind_rows(m1, b1)
x1$site <- factor(x1$site, levels = c("DC","WI","BE", "EI","RP","NB","CB", "BI", "CC"))
x1$time <- as.factor(x1$time)

# intra-site raup-crick summary
intra_raup <- x1 %>%
  group_by(time) %>%
  summarise(mean(distance))

# intra-site raup-crick variation
intra_raup_var <- x1 %>%
  group_by(time) %>%
  summarise(sd(distance))


# Create palette for plot
bluepal <- c("#08519c", "#2271b5", "#4292c6", "#6baed6", "#9ecae1", "#c6dbef", "#deebf7", "#f7fbff", "#ffffff")


# Plot all times and sites
ggplot(x1, aes(site, distance)) + 
  geom_violin(trim = TRUE, aes(fill = factor(site))) +
  scale_fill_manual(values = bluepal) +
  geom_boxplot(width = 0.1) +
  facet_grid(time~.) +
  scale_y_continuous(breaks=seq(-1,1, 1)) +
  geom_hline(yintercept=0, color = "red", linetype = "dashed") +
  ylab("Raup-Crick Probability") +
  xlab("Site") +
  theme(legend.position = "none") + 
  theme_classic() +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))

#1100 x 500 best size

# Compare within meadow variation to between meadow variation through time:

x0 <- x1[,3:5]
x0$loc <- rep("within")

t0 <- t1
t0$loc <- rep("across")

y1 <- bind_rows(x0,t0)
y1$time <- as.factor(y1$time)
y1$loc <- as.factor(y1$loc)  

# summary of raup-crick within/between

raup_crick_summary <- y1 %>%
  group_by(time, loc) %>%
  summarise(mean(distance))

# anova of within/between raup crick data
raup_crick.aov <- aov(distance ~ loc * time, data = y1)
summary(raup_crick.aov)

#Df Sum Sq Mean Sq F value Pr(>F)    
#loc            1    0.0   0.034   0.144  0.704    
#time           2   21.1  10.556  45.222 <2e-16 ***
#  loc:time       2    0.6   0.309   1.322  0.267    
#Residuals   2316  540.6   0.233                   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#library(car)

# levene test of variance
raup_crick.var <- leveneTest(distance ~ loc * time, data = y1)
#raup_crick.var

#Levene's Test for Homogeneity of Variance (center = median)
#Df F value    Pr(>F)    
#group    5  9.4375 6.409e-09 ***
#  2316                      
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

greens <- c("#e5f5e0", "#a1d99b", "#e5f5e0")

# Across/within site comparison by time
ggplot(y1, aes(loc, distance)) + 
  geom_violin(trim = TRUE, aes(fill = factor(time))) +
  scale_fill_manual(values = greens) +
  geom_boxplot(width = 0.1) +
  facet_grid(time~.) +
  scale_y_continuous(breaks=seq(-1,1, 1)) +
  geom_hline(yintercept=0, color = "red", linetype = "dashed") +
  ylab("Raup-Crick Distance") +
  xlab("Scale of Comparison") +
  theme(legend.position = "none") + 
  theme_classic() +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))

# 500 x 600

#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#


#######

#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#


#######



x_med <- y1 %>%
  group_by(time, loc) %>%
  summarise(median(distance))

x_min <- y1 %>%
  group_by(time, loc) %>%
  summarise(min(distance))

x_max <- y1 %>%
  group_by(time, loc) %>%
  summarise(max(distance))

x_mean <- y1 %>%
  group_by(time,loc) %>%
  summarise(mean(distance))


#### MDS for all data

newMDScomm <- new.comm[c(1:32,49:64,80:95,112:127,144:159,176:191,223:238,271:287),]
newMDScomm <- newMDScomm[rowSums(newMDScomm[,3:48])!=0,]
commMDS <- new.comm[c(1:32,49:64,80:95,112:127,144:159,176:191,223:238,271:287),3:48]
commMDS <- commMDS[rowSums(commMDS)!=0,]

palette(c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999"))
newMDS <- metaMDS(commMDS)
MDSpoints <- as.data.frame(newMDS$points)
plot_data <- bind_cols(newMDScomm, MDSpoints)
sitetime <- as.factor(plot_data$sitetime)
plot(plot_data$MDS1, plot_data$MDS2, pch = 19, cex = 2 ,col=sitetime, lwd = 2, xlab = "MDS1", ylab = "MDS2")
legend(-1.7, -0.2, c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"), bty = "n", pch = 19, col = c("#ff7f00", "#999999", "#e41a1c", "#ffff33", "#f781bf", "#a65628", "#4daf4a", "#377eb8", "#984ea3"))
col = c("#ff7f00", "#999999", "#e41a1c", "#ffff33", "#f781bf", "#a65628", "#4daf4a", "#377eb8", "#984ea3")


library(RColorBrewer)

####################### BRAY-CURTIS BETA FOR ALL SITES/TIMES


#make sitetime column and sum data per sample w/in site
full.comm <- read.csv("rawcomm_20170103.csv")

beta.new.comm <- full.comm %>%
  select(-date, -Sieve) %>%
  unite(sitetime, site, Time.Code, sep = "", remove = FALSE)
#fix CBC typo
beta.new.comm$sitetime <- gsub("CBC ", "CBC", beta.new.comm$sitetime)
# remove rows that are empty
beta.new.comm <- beta.new.comm[rowSums(beta.new.comm[,5:50])!=0,]
#finish summarization
beta.new.comm <- beta.new.comm %>%
  select(-Time.Code2) %>%
  group_by(sitetime, site, Time.Code, Sample) %>%
  summarise_each(funs(sum))

#WITHIN SITES

DCA <- beta.new.comm %>%
  subset(sitetime == "DCA") 
DCC <- beta.new.comm %>%
  subset(sitetime == "DCC")
DCE <- beta.new.comm %>%
  subset(sitetime == "DCE")
WIA <- beta.new.comm %>%
  subset(sitetime == "WIA")
WIC <- beta.new.comm %>%
  subset(sitetime == "WIC2")
WIE <- beta.new.comm %>%
  subset(sitetime == "WIE")
RPA <- beta.new.comm %>%
  subset(sitetime == "RPA")
RPC <- beta.new.comm %>%
  subset(sitetime == "RPC")
RPE <- beta.new.comm %>%
  subset(sitetime == "RPE")
NBA <- beta.new.comm %>%
  subset(sitetime == "NBA")
NBC <- beta.new.comm %>%
  subset(sitetime == "NBC")
NBE <- beta.new.comm %>%
  subset(sitetime == "NBE")
CBA <- beta.new.comm %>%
  subset(sitetime == "CBA")
CBC <- beta.new.comm %>%
  subset(sitetime == "CBC")
CBE <- beta.new.comm %>%
  subset(sitetime == "CBE")
BEB <- beta.new.comm %>%
  subset(sitetime == "BEB")
EID <- beta.new.comm %>%
  subset(sitetime == "EID")
BIB <- beta.new.comm %>%
  subset(sitetime == "BIB")
CCD <- beta.new.comm %>%
  subset(sitetime == "CCD")

# Beta as raw alpha/gamma

DCA_rawbeta <- ncol(DCA[,5:50])/mean(specnumber(DCA[,5:50])) - 1
DCC_rawbeta <- ncol(DCC[,5:50])/mean(specnumber(DCC[,5:50])) - 1
DCE_rawbeta <- ncol(DCE[,5:50])/mean(specnumber(DCE[,5:50])) - 1

WIA_rawbeta <- ncol(WIA[,5:50])/mean(specnumber(WIA[,5:50])) - 1
WIC_rawbeta <- ncol(WIC[,5:50])/mean(specnumber(WIC[,5:50])) - 1
WIE_rawbeta <- ncol(WIE[,5:50])/mean(specnumber(WIE[,5:50])) - 1

RPA_rawbeta <- ncol(RPA[,5:50])/mean(specnumber(RPA[,5:50])) - 1
RPC_rawbeta <- ncol(RPC[,5:50])/mean(specnumber(RPC[,5:50])) - 1
RPE_rawbeta <- ncol(RPE[,5:50])/mean(specnumber(RPE[,5:50])) - 1

NBA_rawbeta <- ncol(NBA[,5:50])/mean(specnumber(NBA[,5:50])) - 1
NBC_rawbeta <- ncol(NBC[,5:50])/mean(specnumber(NBC[,5:50])) - 1
NBE_rawbeta <- ncol(NBE[,5:50])/mean(specnumber(NBE[,5:50])) - 1

CBA_rawbeta <- ncol(CBA[,5:50])/mean(specnumber(CBA[,5:50])) - 1
CBC_rawbeta <- ncol(CBC[,5:50])/mean(specnumber(CBC[,5:50])) - 1
CBE_rawbeta <- ncol(CBE[,5:50])/mean(specnumber(CBE[,5:50])) - 1

BIB_rawbeta <- ncol(BIB[,5:50])/mean(specnumber(BIB[,5:50])) - 1

BEB_rawbeta <- ncol(BEB[,5:50])/mean(specnumber(BEB[,5:50])) - 1

EID_rawbeta <- ncol(EID[,5:50])/mean(specnumber(EID[,5:50])) - 1

CCD_rawbeta <- ncol(CCD[,5:50])/mean(specnumber(CCD[,5:50])) - 1

# Beta as Bray Curtis

DCA.mat <- DCA[,5:50] %>%
  vegdist(method = "bray")
DCA_BC <- mean(DCA.mat)
DCC.mat <- DCC[,5:50] %>%
  vegdist(method = "bray")
DCC_BC <- mean(DCC.mat)
DCE.mat <- DCE[,5:50] %>%
  vegdist(method = "bray")
DCE_BC <- mean(DCE.mat)

WIA.mat <- WIA[,5:50] %>%
  vegdist(method = "bray")
WIA_BC <- mean(WIA.mat)
WIC.mat <- WIC[,5:50] %>%
  vegdist(method = "bray")
WIC_BC <- mean(WIC.mat)
WIE.mat <- WIE[,5:50] %>%
  vegdist(method = "bray")
WIE_BC <- mean(WIE.mat)

RPA.mat <- RPA[,5:50] %>%
  vegdist(method = "bray")
RPA_BC <- mean(RPA.mat)
RPC.mat <- RPC[,5:50] %>%
  vegdist(method = "bray")
RPC_BC <- mean(RPC.mat)
RPE.mat <- RPE[,5:50] %>%
  vegdist(method = "bray")
RPE_BC <- mean(RPE.mat)

NBA.mat <- NBA[,5:50] %>%
  vegdist(method = "bray")
NBA_BC <- mean(NBA.mat)
NBC.mat <- NBC[,5:50] %>%
  vegdist(method = "bray")
NBC_BC <- mean(NBC.mat)
NBE.mat <- NBE[,5:50] %>%
  vegdist(method = "bray")
NBE_BC <- mean(NBE.mat)

CBA.mat <- CBA[,5:50] %>%
  vegdist(method = "bray")
CBA_BC <- mean(CBA.mat)
CBC.mat <- CBC[,5:50] %>%
  vegdist(method = "bray")
CBC_BC <- mean(CBC.mat)
CBE.mat <- CBE[,5:50] %>%
  vegdist(method = "bray")
CBE_BC <- mean(CBE.mat)

BEB.mat <- BEB[,5:50] %>%
  vegdist(method = "bray")
BEB_BC <- mean(BEB.mat)

EID.mat <- EID[,5:50] %>%
  vegdist(method = "bray")
EID_BC <- mean(EID.mat)

BIB.mat <- BIB[,5:50] %>%
  vegdist(method = "bray")
BIB_BC <- mean(BIB.mat)

CCD.mat <- CCD[,5:50] %>%
  vegdist(method = "bray")
CCD_BC <- mean(CCD.mat)


Sites <- rep(c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC"),3)
Times <- c(rep("May", 9), rep("June/July",9), rep("August",9))
rawbeta <- c(DCA_rawbeta, WIA_rawbeta, "NA", "NA", RPA_rawbeta, NBA_rawbeta, CBA_rawbeta, "NA", "NA", DCC_rawbeta, WIC_rawbeta, BEB_rawbeta, EID_rawbeta, RPC_rawbeta, NBC_rawbeta, CBC_rawbeta, BIB_rawbeta, CCD_rawbeta, DCE_rawbeta, WIE_rawbeta, "NA", "NA", RPE_rawbeta, NBE_rawbeta, CBE_rawbeta, "NA", "NA" )
BCbeta <- c(DCA_BC, WIA_BC, "NA", "NA", RPA_BC, NBA_BC, CBA_BC, "NA", "NA", DCC_BC, WIC_BC, BEB_BC, EID_BC, RPC_BC, NBC_BC, CBC_BC, BIB_BC, CCD_BC, DCE_BC, WIE_BC, "NA", "NA", RPE_BC, NBE_BC, CBE_BC, "NA", "NA")
  
  
Raw_beta <- data.frame(Sites, Times, rawbeta, BCbeta)
Raw_beta$rawbeta <- as.numeric(as.character(Raw_beta$rawbeta))
Raw_beta$BCbeta <- as.numeric(as.character(Raw_beta$BCbeta))
Raw_beta$Sites <- factor(Raw_beta$Sites, levels = c("DC","WI","BE", "EI","RP","NB","CB", "BI", "CC"))


# export beta values as csv
#write.csv(Raw_beta, file = "beta_all.csv")

# raw beta plot
ggplot(na.omit(Raw_beta), aes(x = Times, y = rawbeta, color = Sites, group = Sites)) +
  geom_point(size=4) + geom_line() +
  xlab("Time Period") + ylab("beta = gamma/mean(alpha)") + theme_classic() + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))

# Bray Curtis beta plot
ggplot(na.omit(Raw_beta), aes(x = Times, y = BCbeta, color = Sites, group = Sites)) +
  geom_point(size=4) + geom_line() +
  xlab("Time Period") + ylab("beta = mean dissimilarity") + theme_classic() + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))
 

# Other Bray Curtis Plot PCoA

target <- c("B","C","C2","D")
BrayCommB <- beta.new.comm %>%
  filter(Time.Code %in% target)
braycomm_mat <- vegdist(BrayCommB[!names(BrayCommB)%in%c("sitetime","site","Time.Code","Sample")], method = "bray")
locB <- BrayCommB$site
brayB <- with(BrayCommB,betadisper(braycomm_mat,locB))
#scores(brayE)
#TukeyHSD(brayE)
plot(brayB, main = "",ylim = c(-0.4,0.5), xlim = c(-0.5,0.5))
text(-0.45, 0.3, "DC")
text(-0.35, -0.2, "RP")
text(-0.2, -0.4, "WI")
text(0.425, -0.2, "CB")
text(0.3, 0.35, "NB")
#eig <- brayE$eig
#eig/sum(eig)
#PCoA1 = 0.36, PCoA2 = 0.15



# Make presence absence matrix from abundance
comm.alphas <- new.comm[,3:48]
comm.alphas <- comm.alphas %>%
  decostand(method = "pa")

# Add ID columns back to PA matrix
comm.alphas$sitetime <- new.comm$sitetime
comm.alphas$Sample <- new.comm$Sample

# Add rowsums for total alpha per sample
comm.alphas$spp_num <- rowSums(comm.alphas[1:46])

# remove pa info
comm.alpha <- comm.alphas[,47:49]

mean_alpha <- comm.alpha %>%
  group_by(sitetime) %>%
  summarize_each(funs(mean))
  
mean_alpha$beta_1 <- 46 / mean_alpha$spp_num

##################################################

#create euclidean distance matrix for samples
euc.mat <- read.csv("dist1.csv")
euc.mat <- as.matrix(euc.mat)
euc.dist <- dist(euc.mat)

#euclidean distance matices for sites w/ missing plots

RPA.euc.mat <- read.csv("RPAdist.csv")
euc.mat <- as.matrix(RPA.euc.mat)
RPA.euc.dist <- dist(RPA.euc.mat)

NBA.euc.mat <- read.csv("NBAdist.csv")
euc.mat <- as.matrix(NBA.euc.mat)
NBA.euc.dist <- dist(RPA.euc.mat)

CBA.euc.mat <- read.csv("CBAdist.csv")
euc.mat <- as.matrix(CBA.euc.mat)
CBA.euc.dist <- dist(CBA.euc.mat)

CBE.euc.mat <- read.csv("CBEdist.csv")
euc.mat <- as.matrix(CBE.euc.mat)
CBE.euc.dist <- dist(CBE.euc.mat)

BEC.euc.mat <- read.csv("BECdist.csv")
euc.mat <- as.matrix(BEC.euc.mat)
BEC.euc.dist <- dist(BEC.euc.mat)

WIC.euc.mat <- read.csv("WICdist.csv")
euc.mat <- as.matrix(WIC.euc.mat)
WIC.euc.dist <- dist(WIC.euc.mat)

#create bray curtis matrices for each sitetime

DCA <- new.comm %>%
  subset(sitetime == "DCA") 
DCC <- new.comm %>%
  subset(sitetime == "DCC")
DCE <- new.comm %>%
  subset(sitetime == "DCE")
WIA <- new.comm %>%
  subset(sitetime == "WIA")
WIC <- new.comm %>%
  subset(sitetime == "WIC")
WIE <- new.comm %>%
  subset(sitetime == "WIE")
RPA <- new.comm %>%
  subset(sitetime == "RPA")
RPC <- new.comm %>%
  subset(sitetime == "RPC")
RPE <- new.comm %>%
  subset(sitetime == "RPE")
NBA <- new.comm %>%
  subset(sitetime == "NBA")
NBC <- new.comm %>%
  subset(sitetime == "NBC")
NBE <- new.comm %>%
  subset(sitetime == "NBE")
CBA <- new.comm %>%
  subset(sitetime == "CBA")
CBC <- new.comm %>%
  subset(sitetime == "CBC")
CBE <- new.comm %>%
  subset(sitetime == "CBE")
BEC <- new.comm %>%
  subset(sitetime == "BEC")
EIC <- new.comm %>%
  subset(sitetime == "EIC")
BIC <- new.comm %>%
  subset(sitetime == "BIC")
CCC <- new.comm %>%
  subset(sitetime == "CCC")


#create distance matrices for each

DCA.mat <- DCA[,3:48] %>%
  vegdist(method = "bray")
DCC.mat <- DCC[,3:48] %>%
  vegdist(method = "bray")
DCE.mat <- DCE[,3:48] %>%
  vegdist(method = "bray")
WIA.mat <- WIA[,3:48] %>%
  vegdist(method = "bray")
WIC.mat <- WIC[,3:48] %>%
  vegdist(method = "bray")
WIE.mat <- WIE[,3:48] %>%
  vegdist(method = "bray")
RPA.mat <- RPA[,3:48] %>%
  vegdist(method = "bray")
RPC.mat <- RPC[,3:48] %>%
  vegdist(method = "bray")
RPE.mat <- RPE[,3:48] %>%
  vegdist(method = "bray")
NBA.mat <- NBA[,3:48] %>%
  vegdist(method = "bray")
NBC.mat <- NBC[,3:48] %>%
  vegdist(method = "bray")
NBE.mat <- NBE[,3:48] %>%
  vegdist(method = "bray")
CBA.mat <- CBA[,3:48] %>%
  vegdist(method = "bray")
CBC.mat <- CBC[,3:48] %>%
  vegdist(method = "bray")
CBE.mat <- CBE[,3:48] %>%
  vegdist(method = "bray")
BEC.mat <- BEC[,3:48] %>%
  vegdist(method = "bray")
EIC.mat <- EIC[,3:48] %>%
  vegdist(method = "bray")
BIC.mat <- BIC[,3:48] %>%
  vegdist(method = "bray")
CCC.mat <- CCC[,3:48] %>%
  vegdist(method = "bray")

#turn distance matrices into tables for analysis
euc.tab <- as.matrix(euc.dist)
euc.tab <- melt(euc.tab)[melt(lower.tri(euc.tab))$value,]
names(euc.tab) <- c("c1","c2","euc.dist")
#euc.tab

RPA.euc.tab <- as.matrix(RPA.euc.dist)
RPA.euc.tab <- melt(RPA.euc.tab)[melt(lower.tri(RPA.euc.tab))$value,]
names(RPA.euc.tab) <- c("c1","c2","RPA.euc.dist")
#RPA.euc.tab

WIC.euc.tab <- as.matrix(WIC.euc.dist)
WIC.euc.tab <- melt(WIC.euc.tab)[melt(lower.tri(WIC.euc.tab))$value,]
names(WIC.euc.tab) <- c("c1","c2","WIC.euc.dist")
#WIC.euc.tab

BEC.euc.tab <- as.matrix(BEC.euc.dist)
BEC.euc.tab <- melt(BEC.euc.tab)[melt(lower.tri(BEC.euc.tab))$value,]
names(BEC.euc.tab) <- c("c1","c2","BEC.euc.dist")
#BEC.euc.tab

NBA.euc.tab <- as.matrix(NBA.euc.dist)
NBA.euc.tab <- melt(NBA.euc.tab)[melt(lower.tri(NBA.euc.tab))$value,]
names(NBA.euc.tab) <- c("c1","c2","NBA.euc.dist")
#NBA.euc.tab

CBA.euc.tab <- as.matrix(CBA.euc.dist)
CBA.euc.tab <- melt(CBA.euc.tab)[melt(lower.tri(CBA.euc.tab))$value,]
names(CBA.euc.tab) <- c("c1","c2","CBA.euc.dist")
#CBA.euc.tab

CBE.euc.tab <- as.matrix(CBE.euc.dist)
CBE.euc.tab <- melt(CBE.euc.tab)[melt(lower.tri(CBE.euc.tab))$value,]
names(CBE.euc.tab) <- c("c1","c2","CBE.euc.dist")
#CBE.euc.tab



DCA.tab <- as.matrix(DCA.mat)
DCA.tab <- melt(DCA.tab)[melt(lower.tri(DCA.tab))$value,]
names(DCA.tab) <- c("c1","c2","bray.dist")
#DCA.tab
DCA.com <- bind_cols(DCA.tab, euc.tab)
DCA.com <- DCA.com[,c(1:3,6)]
#DCA.com

DCC.tab <- as.matrix(DCC.mat)
DCC.tab <- melt(DCC.tab)[melt(lower.tri(DCC.tab))$value,]
names(DCC.tab) <- c("c1","c2","bray.dist")
#DCC.tab
DCC.com <- bind_cols(DCC.tab, euc.tab)
DCC.com <- DCC.com[,c(1:3,6)]
#DCC.com

DCE.tab <- as.matrix(DCE.mat)
DCE.tab <- melt(DCE.tab)[melt(lower.tri(DCE.tab))$value,]
names(DCE.tab) <- c("c1","c2","bray.dist")
#DCE.tab
DCE.com <- bind_cols(DCE.tab, euc.tab)
DCE.com <- DCE.com[,c(1:3,6)]
#DCE.com

WIA.tab <- as.matrix(WIA.mat)
WIA.tab <- melt(WIA.tab)[melt(lower.tri(WIA.tab))$value,]
names(WIA.tab) <- c("c1","c2","bray.dist")
#WIA.tab
WIA.com <- bind_cols(WIA.tab, euc.tab)
WIA.com <- WIA.com[,c(1:3,6)]
#WIA.com

WIC.tab <- as.matrix(WIC.mat)
WIC.tab <- melt(WIC.tab)[melt(lower.tri(WIC.tab))$value,]
names(WIC.tab) <- c("c1","c2","bray.dist")
#WIC.tab
WIC.com <- bind_cols(WIC.tab, WIC.euc.tab)
WIC.com <- WIC.com[,c(1:3,6)]
#WIC.com

WIE.tab <- as.matrix(WIE.mat)
WIE.tab <- melt(WIE.tab)[melt(lower.tri(WIE.tab))$value,]
names(WIE.tab) <- c("c1","c2","bray.dist")
#WIE.tab
WIE.com <- bind_cols(WIE.tab, euc.tab)
WIE.com <- WIE.com[,c(1:3,6)]
#WIE.com

#Turn RPA euclidean distances into table
RPA.euc.tab <- as.matrix(RPA.euc.dist)
RPA.euc.tab <- melt(RPA.euc.tab)[melt(lower.tri(RPA.euc.tab))$value,]
names(RPA.euc.tab) <- c("c1","c2","euc.dist")
#RPA.euc.tab

RPA.tab <- as.matrix(RPA.mat)
RPA.tab <- melt(RPA.tab)[melt(lower.tri(RPA.tab))$value,]
names(RPA.tab) <- c("c1","c2","bray.dist")
#RPA.tab
RPA.com <- bind_cols(RPA.tab, RPA.euc.tab)
RPA.com <- RPA.com[,c(1:3,6)]
#RPA.com

RPC.tab <- as.matrix(RPC.mat)
RPC.tab <- melt(RPC.tab)[melt(lower.tri(RPC.tab))$value,]
names(RPC.tab) <- c("c1","c2","bray.dist")
#RPC.tab
RPC.com <- bind_cols(RPC.tab, euc.tab)
RPC.com <- RPC.com[,c(1:3,6)]
#RPC.com

RPE.tab <- as.matrix(RPE.mat)
RPE.tab <- melt(RPE.tab)[melt(lower.tri(RPE.tab))$value,]
names(RPE.tab) <- c("c1","c2","bray.dist")
#RPE.tab
RPE.com <- bind_cols(RPE.tab, euc.tab)
RPE.com <- RPE.com[,c(1:3,6)]
#RPE.com

NBA.tab <- as.matrix(NBA.mat)
NBA.tab <- melt(NBA.tab)[melt(lower.tri(NBA.tab))$value,]
names(NBA.tab) <- c("c1","c2","bray.dist")
#NBA.tab
NBA.com <- bind_cols(NBA.tab, NBA.euc.tab)
NBA.com <- NBA.com[,c(1:3,6)]
#NBA.com

NBC.tab <- as.matrix(NBC.mat)
NBC.tab <- melt(NBC.tab)[melt(lower.tri(NBC.tab))$value,]
names(NBC.tab) <- c("c1","c2","bray.dist")
#NBC.tab
NBC.com <- bind_cols(NBC.tab, euc.tab)
NBC.com <- NBC.com[,c(1:3,6)]
#NBC.com

NBE.tab <- as.matrix(NBE.mat)
NBE.tab <- melt(NBE.tab)[melt(lower.tri(NBE.tab))$value,]
names(NBE.tab) <- c("c1","c2","bray.dist")
#NBE.tab
NBE.com <- bind_cols(NBE.tab, euc.tab)
NBE.com <- NBE.com[,c(1:3,6)]
#NBE.com

CBA.tab <- as.matrix(CBA.mat)
CBA.tab <- melt(CBA.tab)[melt(lower.tri(CBA.tab))$value,]
names(CBA.tab) <- c("c1","c2","bray.dist")
#CBA.tab
CBA.com <- bind_cols(CBA.tab, CBA.euc.tab)
CBA.com <- CBA.com[,c(1:3,6)]
#CBA.com

CBC.tab <- as.matrix(CBC.mat)
CBC.tab <- melt(CBC.tab)[melt(lower.tri(CBC.tab))$value,]
names(CBC.tab) <- c("c1","c2","bray.dist")
#CBC.tab
CBC.com <- bind_cols(CBC.tab, euc.tab)
CBC.com <- CBC.com[,c(1:3,6)]
#CBC.com

CBE.tab <- as.matrix(CBE.mat)
CBE.tab <- melt(CBE.tab)[melt(lower.tri(CBE.tab))$value,]
names(CBE.tab) <- c("c1","c2","bray.dist")
#CBE.tab
CBE.com <- bind_cols(CBE.tab, CBE.euc.tab)
CBE.com <- CBE.com[,c(1:3,6)]
#CBE.com

CCC.tab <- as.matrix(CCC.mat)
CCC.tab <- melt(CCC.tab)[melt(lower.tri(CCC.tab))$value,]
names(CCC.tab) <- c("c1","c2","bray.dist")
#CCC.tab
CCC.com <- bind_cols(CCC.tab, euc.tab)
CCC.com <- CCC.com[,c(1:3,6)]
#CCC.com

BIC.tab <- as.matrix(BIC.mat)
BIC.tab <- melt(BIC.tab)[melt(lower.tri(BIC.tab))$value,]
names(BIC.tab) <- c("c1","c2","bray.dist")
#BIC.tab
BIC.com <- bind_cols(BIC.tab, euc.tab)
BIC.com <- BIC.com[,c(1:3,6)]
#BIC.com

EIC.tab <- as.matrix(EIC.mat)
EIC.tab <- melt(EIC.tab)[melt(lower.tri(EIC.tab))$value,]
names(EIC.tab) <- c("c1","c2","bray.dist")
#EIC.tab
EIC.com <- bind_cols(EIC.tab, euc.tab)
EIC.com <- EIC.com[,c(1:3,6)]
#EIC.com

BEC.tab <- as.matrix(BEC.mat)
BEC.tab <- melt(BEC.tab)[melt(lower.tri(BEC.tab))$value,]
names(BEC.tab) <- c("c1","c2","bray.dist")
#BEC.tab
BEC.com <- bind_cols(BEC.tab, BEC.euc.tab)
BEC.com <- BEC.com[,c(1:3,6)]
#BEC.com

#lm and plots

plot(euc.dist, DCA.mat, main = "DCA", ylim = c(0,1))
DCA.lm <- lm(bray.dist ~ euc.dist, data = DCA.com)
abline(DCA.lm)
text(4, 0.95, "p = 0.09")
text(3.8, 0.9, "adj R^2 = 0.02")

plot(euc.dist, DCC.mat, main = "DCC", ylim = c(0,1))
DCC.lm <- lm(bray.dist ~ euc.dist, data = DCC.com)
abline(DCC.lm)
text(4, 0.95, "p < 0.0005")
text(3.8, 0.9, "adj R^2 = 0.13")

plot(euc.dist, DCE.mat, main = "DCE", ylim = c(0,1))
DCE.lm <- lm(bray.dist ~ euc.dist, data = DCE.com)
abline(DCE.lm)
text(4, 0.95, "p = 0.58")
text(3.8, 0.9, "adj R^2 < 0.01")

plot(euc.dist, WIA.mat, main = "WIA", ylim = c(0,1))
WIA.lm <- lm(bray.dist ~ euc.dist, data = WIA.com)
abline(WIA.lm)
text(4, 0.95, "p = 0.052")
text(3.8, 0.9, "adj R^2 = 0.02")


plot(WIC.euc.dist, WIC.mat, main = "WIC", ylim = c(0,1))
WIC.lm <- lm(bray.dist ~ WIC.euc.dist, data = WIC.com)
abline(WIC.lm)
text(4, 0.95, "p = 0.06")
text(3.8, 0.9, "adj R^2 = 0.02")


plot(euc.dist, WIE.mat, main = "WIE", ylim = c(0,1))
WIE.lm <- lm(bray.dist ~ euc.dist, data = WIE.com)
abline(WIE.lm)
text(4, 0.95, "p = 0.07")
text(3.8, 0.9, "adj R^2 = 0.02")

plot(RPA.euc.dist, RPA.mat, main = "RPA", ylim = c(0,1))
RPA.lm <- lm(bray.dist ~ RPA.euc.dist, data = RPA.com)
abline(RPA.lm)
text(4, 0.95, "p = 0.13")
text(3.8, 0.9, "adj R^2 = 0.01")

plot(euc.dist, RPC.mat, main = "RPC", ylim = c(0,1))
RPC.lm <- lm(bray.dist ~ euc.dist, data = RPC.com)
abline(RPC.lm)
text(4, 0.95, "p = 0.23")
text(3.8, 0.9, "adj R^2 < 0.01")


plot(euc.dist, RPE.mat, main = "RPE", ylim = c(0,1))
RPE.lm <- lm(bray.dist ~ euc.dist, data = RPE.com)
abline(RPE.lm)
text(4, 0.95, "p = 0.36")
text(3.8, 0.9, "adj R^2 < 0.01")


plot(euc.dist, BIC.mat, main = "BIC", ylim = c(0,1))
BIC.lm <- lm(bray.dist ~ euc.dist, data = BIC.com)
abline(BIC.lm)
text(4, 0.95, "p = 0.8")
text(3.8, 0.9, "adj R^2 < 0.01")


plot(euc.dist, CCC.mat, main = "CCC", ylim = c(0,1))
CCC.lm <- lm(bray.dist ~ euc.dist, data = CCC.com)
abline(CCC.lm)
text(4, 0.95, "p = 0.35")
text(3.8, 0.9, "adj R^2 < 0.01")


plot(euc.dist, EIC.mat, main = "EIC", ylim = c(0,1))
EIC.lm <- lm(bray.dist ~ euc.dist, data = EIC.com)
abline(EIC.lm)
text(4, 0.95, "p = 0.002")
text(3.8, 0.9, "adj R^2 = 0.07")

plot(BEC.euc.dist, BEC.mat, main = "BEC", ylim = c(0,1))
BEC.lm <- lm(bray.dist ~ BEC.euc.dist, data = BEC.com)
abline(BEC.lm)
text(4, 0.95, "p = 0.10")
text(3.8, 0.9, "adj R^2 = 0.02")

plot(NBA.euc.dist, NBA.mat, main = "NBA", ylim = c(0,1))
NBA.lm <- lm(bray.dist ~ NBA.euc.dist, data = NBA.com)
abline(NBA.lm)
text(4, 0.95, "p = 0.78")
text(3.8, 0.9, "adj R^2 < 0.01")


plot(euc.dist, NBC.mat, main = "NBC", ylim = c(0,1))
NBC.lm <- lm(bray.dist ~ euc.dist, data = NBC.com)
abline(NBC.lm)
text(4, 0.95, "p = 0.41")
text(3.8, 0.9, "adj R^2 < 0.01")


plot(euc.dist, NBE.mat, main = "NBE", ylim = c(0,1))
NBE.lm <- lm(bray.dist ~ euc.dist, data = NBE.com)
abline(NBE.lm)
text(4, 0.95, "p = 0.32")
text(3.8, 0.9, "adj R^2 < 0.01")


plot(CBA.euc.dist, CBA.mat, main = "CBA", ylim = c(0,1), xlim = c(1,4.3))
CBA.lm <- lm(bray.dist ~ CBA.euc.dist, data = CBA.com)
abline(CBA.lm)
text(4, 0.95, "p = 0.28")
text(3.8, 0.9, "adj R^2 < 0.01")


plot(euc.dist, CBC.mat, main = "CBC", ylim = c(0,1))
CBC.lm <- lm(bray.dist ~ euc.dist, data = CBC.com)
abline(CBC.lm)
text(4, 0.9, "p = 0.81")
text(3.8, 0.85, "adj R^2 < 0.01")


plot(CBE.euc.dist, CBE.mat, main = "CBE", ylim = c(0,1))
CBE.lm <- lm(bray.dist ~ CBE.euc.dist, data = CBE.com)
abline(CBE.lm)
text(4, 0.95, "p = 0.91")
text(3.8, 0.9, "adj R^2 < 0.01")













######DISTANCE TO CENTROID FOR ALL TIME PERIODS INDIVIDUALLY

setwd("~/Desktop/Add To Personal/Beta_analysis_2015-4")
library(dplyr)
library(tidyr)
library(reshape2)
library(lme4)
library(ggplot2)
library(vegan)

#make sitetime column and sum data per sample w/in site
full.comm <- read.csv("rawcomm-correctedsamplenumbers.csv")
new.comm <- full.comm %>%
  select(-date, -Sieve) %>%
  unite(sitetime, site, Time.Code, sep = "", remove = TRUE)
new.comm <- new.comm %>%
  group_by(sitetime, Sample) %>%
  summarise_each(funs(sum))

#distance matrix of all sites and times, remove empty rows
new.comm$sums <- rowSums(new.comm[,3:48]) #%>%
new.comm <- filter(new.comm, sums > 0) #%>%
new.comm <- select(new.comm, -sums)
full.mat <- new.comm[,3:48] #%>%
  
full.dist <- vegdist(full.mat, method = "bray")
all.loc <- new.comm$sitetime
bray.all <- with(new.comm, betadisper(full.dist,all.loc))



#ordination for all sitetimes
plot(bray.all)

axes <- c(1,2)
g <- scores (bray.all, choices = axes)
centroids <- points(g$centroids)
ch <- chull(g$sites)
ch <- c(ch, ch[1])
polygons <- lines(bray.all$vectors[, axes][ch, ])

group.color <- as.data.frame(levels(bray.all$group))
group.color$color <- 
  #rainbow(19, alpha = .15)
  #c("gray","gray","lightgray","gray","darkgray","gray","lightgray","gray","darkgray","gray","lightgray","gray","darkgray","lightgray","gray","darkgray","lightgray","gray","darkgray")
colnames(group.color) <- c("sitetime", "color")
plot(g$centroids, , pch = 16, cex = 2, col = "white", xlim = c(-.5,.5), ylim = c(-.5,.5))
text(g$centroids, label = group.color$sitetime, col = "black", cex = 1.5)


for (i in levels(bray.all$group)) {
  j <- which(levels(bray.all$group) == i)
  
    ch <- chull(g$sites[bray.all$group == i, ])
    ch <- c(ch, ch[1])
    polygon(bray.all$vectors[bray.all$group == i, axes][ch, ], 
            col = group.color[i==levels(bray.all$group),2], lty =0)
}



#distance matrix of May############
A.comm <- new.comm[grep("May", new.comm$sitetime),] 

#reorder factors
A.comm$sitetime <- factor(A.comm$sitetime, levels = c("DCA", "WIA", "RPA", "NBA", "CBA")) 

#Bray Curtis matrix
A.bray <- vegdist(A.comm[,3:48])

#Match sample ID with matrix
A.loc <- A.comm$sitetime
bray.A <- with(A.comm, betadisper(A.bray,A.loc))

#ordination
plot(bray.A)

#distance to centroid plot
boxplot(bray.A, ylab = "mean multivariate dispersion", main = "May", xaxt = 'n',
        ylim = c(0:1))
asites <- c("DC","WI","RP","NB","CB")
axis(1, at=1:5, labels=asites)

#permutest (better than TUKEY)
permutest(bray.A, pairwise = TRUE)

#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
#Groups     4 0.67601 0.169002 11.017    999  0.001 ***
#  Residuals 71 1.08919 0.015341                         
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#       DCA        WIA        RPA        NBA   CBA
#DCA            1.0000e-03 6.2000e-01 5.9000e-02 0.761
#WIA 4.2235e-05            1.0000e-03 1.6000e-02 0.001
#RPA 6.2424e-01 1.9848e-06            1.7000e-02 0.803
#NBA 6.9445e-02 1.8891e-02 1.7008e-02            0.029
#CBA 7.8694e-01 7.4334e-06 8.0636e-01 3.2291e-02 

text(1, .7, "AB")
text(2, .85, "C")
text(3, .55, "A")
text(4, .83, "B")
text(5, .58, "A")





#distance matrix of June/July##########

#create subset and (filter out non C sites)
C.comm <- new.comm[grep("C", new.comm$sitetime),] 
C.comm <- filter(C.comm, !grepl('DCA|DCE|CBA|CBE', sitetime))

#reorder factors
C.comm$sitetime <- factor(C.comm$sitetime, levels = c("DCC", "WIC", "BEC", "EIC", "RPC",
                                                      "NBC","CBC","BIC","CCC")) 

#Bray Curtis matrix
C.bray <- vegdist(C.comm[,3:48])

#Match sample ID with matrix
C.loc <- C.comm$sitetime
bray.C <- with(C.comm, betadisper(C.bray,C.loc))

#ordination
plot(bray.C)

#distance to centroid plot
boxplot(bray.C, ylab = "mean multivariate dispersion", main = "June/July", xaxt = 'n',
        ylim = c(0:1))
csites <- c("DC","WI","BE","EI","RP","NB","CB","BI","CC")
axis(1, at=1:9, labels=csites)

#permutest (better than TUKEY)
permutest(bray.C, pairwise = TRUE)

#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
#Groups      8 0.46314 0.057893 2.4734    999  0.013 *
#  Residuals 133 3.11307 0.023407                       
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#         DCC       WIC       BEC       EIC       RPC       NBC       CBC       BIC   CCC
#DCC           0.0240000 0.3930000 0.9650000 0.2930000 0.1000000 0.4190000 0.0630000 0.055
#WIC 0.0204886           0.0040000 0.0190000 0.2530000 0.4390000 0.1660000 0.5540000 0.505
#BEC 0.3745045 0.0014604           0.2950000 0.0670000 0.0120000 0.1170000 0.0050000 0.004
#EIC 0.9667559 0.0108386 0.2886707           0.3000000 0.0630000 0.4080000 0.0310000 0.033
#RPC 0.3046187 0.2449232 0.0673714 0.2772511           0.5830000 0.8220000 0.4760000 0.478
#NBC 0.0917767 0.4372296 0.0090554 0.0644647 0.6133337           0.4340000 0.7930000 0.856
#CBC 0.4407222 0.1640547 0.1169214 0.4202317 0.8156339 0.4496836           0.3280000 0.358
#BIC 0.0508497 0.5606710 0.0035834 0.0305234 0.4695852 0.8141242 0.3287350           0.951
#CCC 0.0478288 0.5011401 0.0026881 0.0264998 0.4876437 0.8558280 0.3386071 0.9457608       




#distance matrix of August##########

#create subset and (filter out non E sites)
E.comm <- new.comm[grep("C", new.comm$sitetime),] 
E.comm <- filter(E.comm, !grepl('BEC|EIC', sitetime))

#reorder factors
E.comm$sitetime <- factor(E.comm$sitetime, levels = c("DCE", "WIE", "RPE", "NBE", "CBE")) 

#Bray Curtis matrix
E.bray <- vegdist(E.comm[,3:48])

#Match sample ID with matrix
E.loc <- E.comm$sitetime
bray.E <- with(E.comm, betadisper(E.bray,E.loc))

#ordination
plot(bray.E)

#distance to centroid plot
boxplot(bray.E, ylab = "mean multivariate dispersion", main = "August", xaxt = 'n',
        ylim = c(0:1))
esites <- c("DC","WI","RP","NB","CB")
axis(1, at=1:5, labels=esites)

#permutest (better than TUKEY)
permutest(bray.E, pairwise = TRUE)

#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
#Groups     4 0.44428 0.111069 6.4316    999  0.002 **
#  Residuals 74 1.27793 0.017269                        
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#        DCE        WIE        RPE        NBE   CBE
#DCE            1.7000e-01 5.4800e-01 1.5000e-02 0.001
#WIE 1.5202e-01            3.3600e-01 3.9600e-01 0.015
#RPE 5.2008e-01 3.1482e-01            3.4000e-02 0.001
#NBE 9.9434e-03 3.7732e-01 2.8242e-02            0.083
#CBE 1.6841e-05 1.4976e-02 4.3381e-05 7.4228e-02    

text(1, .55, "A")
text(2, .7, "AB")
text(3, .57, "A")
text(4, .7, "BC")
text(5, .78, "C")


#par(mfrow = c(1,3))


par(mfrow = c(1,1))

#comparing null distance to centroid expectations to observations, each site/time
#using mean rarefied richness and mean abundance for each sample, log normal dist
par(mfrow = c(5,2))
par(mar = c(4,4.1,2,0.5))
boxplot(DCA.mat, ylim = c(0,1), xlab = "DCA")
sad_beta(9,70, Nsite = 1)
boxplot(WIA.mat, ylim = c(0,1), xlab = "WIA")
sad_beta(12,25, Nsite = 1)
boxplot(RPA.mat, ylim = c(0,1), ylab = "dist to centroid", xlab = "RPA")
sad_beta(8,52, Nsite = 1)
boxplot(NBA.mat, ylim = c(0,1), xlab = "NBA")
sad_beta(9, 40, Nsite = 1)
boxplot(CBA.mat, ylim = c(0,1), xlab = "CBA")
sad_beta(7, 38, Nsite = 1)

par(mfrow = c(5,4))
boxplot(DCC.mat, ylim = c(0,1), xlab = "DCC")
sad_beta(7, 400, Nsite = 1)
boxplot(WIC.mat, ylim = c(0,1), xlab = "WIC")
sad_beta(10,100, Nsite = 1)
boxplot(BEC.mat, ylim = c(0,1), xlab = "BEC")
sad_beta(9,100, Nsite = 1)
boxplot(EIC.mat, ylim = c(0,1), xlab = "EIC")
sad_beta(10, 500, Nsite = 1)
boxplot(RPC.mat, ylim = c(0,1), xlab = "RPC", ylab = "dist to centroid")
sad_beta(16, 300, Nsite = 1)
boxplot(NBC.mat, ylim = c(0,1), xlab = "NBC")
sad_beta(8, 100, Nsite = 1)
boxplot(CBC.mat, ylim = c(0,1), xlab = "CBC")
sad_beta(9, 200, Nsite = 1)
boxplot(BIC.mat, ylim = c(0,1), xlab = "BIC")
sad_beta(12, 100, Nsite = 1)
boxplot(CCC.mat, ylim = c(0,1), xlab = "CCC")
sad_beta(7, 50, Nsite = 1)

par(mfrow = c(5,2))
par(mar = c(4,4.1,2,0.5))
boxplot(DCE.mat, ylim = c(0,1), xlab = "DCE")
sad_beta(8,1300, Nsite = 1)
boxplot(WIE.mat, ylim = c(0,1), xlab = "WIE")
sad_beta(12,500, Nsite = 1)
boxplot(RPE.mat, ylim = c(0,1), ylab = "dist to centroid", xlab = "RPE")
sad_beta(8,750, Nsite = 1)
boxplot(NBE.mat, ylim = c(0,1), xlab = "NBE")
sad_beta(7, 100, Nsite = 1)
boxplot(CBE.mat, ylim = c(0,1), xlab = "CBE")
sad_beta(9, 100, Nsite = 1)


#comparing null distance to centroid expectations to observations, each site/time
#using mean rarefied richness and mean abundance for each sample, geometric dist
par(mfrow = c(5,2))
par(mar = c(4,4.1,2,0.5))
boxplot(DCA.mat, ylim = c(0,1), xlab = "DCA")
sad_beta(9,70, Nsite = 1, dist = "geom")
boxplot(WIA.mat, ylim = c(0,1), xlab = "WIA")
sad_beta(12,25, Nsite = 1, dist = "geom")
boxplot(RPA.mat, ylim = c(0,1), ylab = "dist to centroid", xlab = "RPA")
sad_beta(8,52, Nsite = 1, dist = "geom")
boxplot(NBA.mat, ylim = c(0,1), xlab = "NBA")
sad_beta(9, 40, Nsite = 1, dist = "geom")
boxplot(CBA.mat, ylim = c(0,1), xlab = "CBA")
sad_beta(7, 38, Nsite = 1, dist = "geom")

par(mfrow = c(5,4))
boxplot(DCC.mat, ylim = c(0,1), xlab = "DCC")
sad_beta(7, 400, Nsite = 1, dist = "geom")
boxplot(WIC.mat, ylim = c(0,1), xlab = "WIC")
sad_beta(10,100, Nsite = 1, dist = "geom")
boxplot(BEC.mat, ylim = c(0,1), xlab = "BEC")
sad_beta(9,100, Nsite = 1, dist = "geom")
boxplot(EIC.mat, ylim = c(0,1), xlab = "EIC")
sad_beta(10, 500, Nsite = 1, dist = "geom")
boxplot(RPC.mat, ylim = c(0,1), xlab = "RPC", ylab = "dist to centroid")
sad_beta(16, 300, Nsite = 1, dist = "geom")
boxplot(NBC.mat, ylim = c(0,1), xlab = "NBC")
sad_beta(8, 100, Nsite = 1, dist = "geom")
boxplot(CBC.mat, ylim = c(0,1), xlab = "CBC")
sad_beta(9, 200, Nsite = 1, dist = "geom")
boxplot(BIC.mat, ylim = c(0,1), xlab = "BIC")
sad_beta(12, 100, Nsite = 1, dist = "geom")
boxplot(CCC.mat, ylim = c(0,1), xlab = "CCC")
sad_beta(7, 50, Nsite = 1, dist = "geom")

par(mfrow = c(5,2))
par(mar = c(4,4.1,2,0.5))
boxplot(DCE.mat, ylim = c(0,1), xlab = "DCE")
sad_beta(8,1300, Nsite = 1, dist = "geom")
boxplot(WIE.mat, ylim = c(0,1), xlab = "WIE")
sad_beta(12,500, Nsite = 1, dist = "geom")
boxplot(RPE.mat, ylim = c(0,1), ylab = "dist to centroid", xlab = "RPE")
sad_beta(8,750, Nsite = 1, dist = "geom")
boxplot(NBE.mat, ylim = c(0,1), xlab = "NBE")
sad_beta(7, 100, Nsite = 1, dist = "geom")
boxplot(CBE.mat, ylim = c(0,1), xlab = "CBE")
sad_beta(9, 100, Nsite = 1, dist = "geom")

#comparing null distance to centroid expectations to observations, each site/time
#using mean rarefied richness and mean abundance for each sample, zipf dist
par(mfrow = c(5,2))
par(mar = c(4,4.1,2,0.5))
boxplot(DCA.mat, ylim = c(0,1), xlab = "DCA")
sad_beta(9,70, Nsite = 1, dist = "zipf")
boxplot(WIA.mat, ylim = c(0,1), xlab = "WIA")
sad_beta(12,25, Nsite = 1, dist = "zipf")
boxplot(RPA.mat, ylim = c(0,1), ylab = "dist to centroid", xlab = "RPA")
sad_beta(8,52, Nsite = 1, dist = "zipf")
boxplot(NBA.mat, ylim = c(0,1), xlab = "NBA")
sad_beta(9, 40, Nsite = 1, dist = "zipf")
boxplot(CBA.mat, ylim = c(0,1), xlab = "CBA")
sad_beta(7, 38, Nsite = 1, dist = "zipf")

par(mfrow = c(5,4))
boxplot(DCC.mat, ylim = c(0,1), xlab = "DCC")
sad_beta(7, 400, Nsite = 1, dist = "zipf")
boxplot(WIC.mat, ylim = c(0,1), xlab = "WIC")
sad_beta(10,100, Nsite = 1, dist = "zipf")
boxplot(BEC.mat, ylim = c(0,1), xlab = "BEC")
sad_beta(9,100, Nsite = 1, dist = "zipf")
boxplot(EIC.mat, ylim = c(0,1), xlab = "EIC")
sad_beta(10, 500, Nsite = 1, dist = "zipf")
boxplot(RPC.mat, ylim = c(0,1), xlab = "RPC", ylab = "dist to centroid")
sad_beta(16, 300, Nsite = 1, dist = "zipf")
boxplot(NBC.mat, ylim = c(0,1), xlab = "NBC")
sad_beta(8, 100, Nsite = 1, dist = "zipf")
boxplot(CBC.mat, ylim = c(0,1), xlab = "CBC")
sad_beta(9, 200, Nsite = 1, dist = "zipf")
boxplot(BIC.mat, ylim = c(0,1), xlab = "BIC")
sad_beta(12, 100, Nsite = 1, dist = "zipf")
boxplot(CCC.mat, ylim = c(0,1), xlab = "CCC")
sad_beta(7, 50, Nsite = 1, dist = "zipf")

par(mfrow = c(5,2))
par(mar = c(4,4.1,2,0.5))
boxplot(DCE.mat, ylim = c(0,1), xlab = "DCE")
sad_beta(8,1300, Nsite = 1, dist = "zipf")
boxplot(WIE.mat, ylim = c(0,1), xlab = "WIE")
sad_beta(12,500, Nsite = 1, dist = "zipf")
boxplot(RPE.mat, ylim = c(0,1), ylab = "dist to centroid", xlab = "RPE")
sad_beta(8,750, Nsite = 1, dist = "zipf")
boxplot(NBE.mat, ylim = c(0,1), xlab = "NBE")
sad_beta(7, 100, Nsite = 1, dist = "zipf")
boxplot(CBE.mat, ylim = c(0,1), xlab = "CBE")
sad_beta(9, 100, Nsite = 1, dist = "zipf")


#full community ordination for beta 

comm.code <- full.comm$Code
full.dis <- with(full.comm,betadisper(full.mat,comm.code))
full.plot <- plot(full.dis, xlim = c(-0.6,0.6))
text(full.plot$centroids, row.names(full.plot$centroids), cex=2, pos=c(3,4,1,4,1,4,2,2,1,4,4,3,2,4,1,2,2,2,4), col="black")
glimpse(full.comm)

#community ordination for time C only (look at clumping of secondaries)

comm.code <- full.comm$Code
full.dis <- with(full.comm,betadisper(full.mat,comm.code))
full.plot <- plot(full.dis, xlim = c(-0.6,0.6))
text(full.plot$centroids, row.names(full.plot$centroids), cex=1, pos=c(3,4,1,4,1,4,2,2,1,4,4,3,2,4,1,2,2,2,4), col="red")
glimpse(full.comm)


#below is source code from thesis

locA <- commsmallA$site
braysmallA <- with(commsmallA,betadisper(braycommsmallA,locA))
#scores(braysmallA)
#TukeyHSD(braysmallA)
par(mar = c(3,2,0,0.5))
plot(braysmallA, main = "", ylim = c(-0.4,0.5), xlim = c(-0.5,0.5))
text(0.4, 0.8, "A", cex = 2)
text(-0.4, -0.1, "DC")
text(-0.25, 0.25, "RP")
text(-0.2, -0.4, "WI")
text(0.47, -0.25, "NB")
text(0.3, 0.45, "CB")
#eigA <- braysmallA$eig
#eigA/sum(eigA)
#PCoA1 = 0.37, PCoA2 = 0.15

braycommsmallE <- vegdist(commsmallE[!names(commsmall)%in%c("site","time","Sample")], method = "bray")
locE <- commsmallE$site
braysmallE <- with(commsmallE,betadisper(braycommsmallE,locE))
#scores(braysmallE)
#TukeyHSD(braysmallE)
plot(braysmallE, main = "", ylim = c(-0.4,0.5), xlim = c(-0.5,0.5), yaxt = 'n')
text(0.4, 0.8, "B", cex = 2)
text(-0.45, 0.3, "DC")
text(-0.25, -0.15, "RP")
text(-0.18, -0.45, "WI")
text(0.45, -0.2, "CB")
text(0.35, 0.35, "NB")
#eig <- braysmallE$eig
#eig/sum(eig)
#PCoA1 = 0.42, PCoA2 = 0.19

@
  \end{center}
\caption{Beta diversity as Bray-Curtis distance for each site in May (A) and August (B). 
         Red points represent a multidimensional centroid for each site, blue vectors are distance of 
         sampled 0.28 m\squared plots within site to the centroid. Axes are dimensionless.}



