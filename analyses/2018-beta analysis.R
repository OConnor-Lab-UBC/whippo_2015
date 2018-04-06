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

library(tidyverse)
library(reshape2)
library(lme4)
library(vegan)
library(lubridate) # data manipulation
library(viridis) # plots

# YOU MUST RUN RAUP_CRICK.R FOR FUNCTIONALITY 

###################################################################################
# READ IN DATA                                                                    #
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
ggplot(t1, aes(time, distance, fill = time)) + 
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1) +
  #fill_palette(viridis(3, option = "D", begin = 0.3)) +
  #facet_grid(site~.) +
  scale_y_continuous(breaks=seq(-1,1, 1)) +
  geom_hline(yintercept=0, color = "red", linetype = "dashed") +
  ylab("Raup-Crick Distance") +
  theme_minimal() +
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

BEC <- epicomm_full %>%
  subset(sitetime == "BEC")
BIC <- epicomm_full %>%
  subset(sitetime == "BIC")
CCC <- epicomm_full %>%
  subset(sitetime == "CCC")
EIC <- epicomm_full %>%
  subset(sitetime == "EIC")

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

BEC_jaccard <- BEC[,2:34] %>%
  vegdist(method = "jaccard")
BIC_jaccard <- BIC[,2:34] %>%
  vegdist(method = "jaccard")
CCC_jaccard <- CCC[,2:34] %>%
  vegdist(method = "jaccard")
EIC_jaccard <- EIC[,2:34] %>%
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

BEC_bray <- BEC[,2:34] %>%
  vegdist(method = "bray")
BIC_bray <- BIC[,2:34] %>%
  vegdist(method = "bray")
CCC_bray <- CCC[,2:34] %>%
  vegdist(method = "bray")
EIC_bray <- EIC[,2:34] %>%
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
# ALL_rc_C <- raup_crick(epicomm_C, plot_names_in_col1 = TRUE)





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

ggplot(r5, aes(time, distance)) + 
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

ggplot(n5, aes(time, distance)) + 
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

ggplot(c5, aes(time, distance)) + 
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

# BEC - shared absences (makes little difference)
BEC_rc  <- epicomm_full[1:16,-1]
BEC_rc <- data.frame(BEC_rc)
BEC_rc_matrix <- raup_crick(BEC_rc, plot_names_in_col1 = FALSE)
B_mat <- as.matrix(BEC_rc_matrix)
b2 <- melt(B_mat)[melt(upper.tri(B_mat))$value,]
names(b2) <- c("c1", "c2", "distance")
b2$time <- rep("June/July")
b2$pair <- seq(from = 1, to=120, by =1)

#BEC_summary <- as.table(summary(BEC_rc_matrix))

# BIC - shared absences (makes little difference)
BIC_rc  <- epicomm_full[17:32,-1]
BIC_rc <- data.frame(BIC_rc)
BIC_rc_matrix <- raup_crick(BIC_rc, plot_names_in_col1 = FALSE)
B_mat <- as.matrix(BIC_rc_matrix)
b3 <- melt(B_mat)[melt(upper.tri(B_mat))$value,]
names(b3) <- c("c1", "c2", "distance")
b3$time <- rep("June/July")
b3$pair <- seq(from = 1, to=120, by =1)

#BIC_summary <- as.table(summary(BIC_rc_matrix))

# CCC - shared absences (makes little difference)
CCC_rc  <- epicomm_full[80:95,-1]
CCC_rc <- data.frame(CCC_rc)
CCC_rc_matrix <- raup_crick(CCC_rc, plot_names_in_col1 = FALSE)
B_mat <- as.matrix(CCC_rc_matrix)
b4 <- melt(B_mat)[melt(upper.tri(B_mat))$value,]
names(b4) <- c("c1", "c2", "distance")
b4$time <- rep("June/July")
b4$pair <- seq(from = 1, to=120, by =1)

#CCC_summary <- as.table(summary(CCC_rc_matrix))

# EIC - shared absences (makes little difference)
EIC_rc  <- epicomm_full[144:159,-1]
EIC_rc <- data.frame(EIC_rc)
EIC_rc_matrix <- raup_crick(EIC_rc, plot_names_in_col1 = FALSE)
B_mat <- as.matrix(EIC_rc_matrix)
b5 <- melt(B_mat)[melt(upper.tri(B_mat))$value,]
names(b5) <- c("c1", "c2", "distance")
b5$time <- rep("June/July")
b5$pair <- seq(from = 1, to=120, by =1)

#EIC_summary <- as.table(summary(EIC_rc_matrix))

# Export all summary values of analysis
#RC_all <- rbind(DCA_summary, DCC_summary, DCE_summary, WIA_summary, WIC_summary, WIE_summary, RPA_summary, RPC_summary, RPE_summary, NBA_summary, NBC_summary, NBE_summary, CBA_summary, CBC_summary, CBE_summary, BEC_summary, BIC_summary, EIC_summary, CCC_summary)
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
x1$time <- factor(x1$time, levels = c("May", "June/July", "August"))

# Plot all times and sites
ggplot(x1, aes(site, distance)) + 
  geom_violin(trim = TRUE, aes(fill = factor(site))) +
  scale_fill_manual(values = bluepal, guide = FALSE) +
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

greens <- c("#BBDF27FF", "#2FB47CFF", "#BBDF27FF")
y1$time <- factor(y1$time, levels = c("May", "June/July", "August"))

# Across/within site comparison by time
ggplot(y1, aes(loc, distance)) + 
  geom_violin(trim = TRUE, aes(fill = factor(time))) +
  scale_fill_manual(values = greens, guide = FALSE) +
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


##### SCRATCH PAD


# actual mean values for Raup-Crick figure
test <- y1
test <- test %>%
  select(time, distance, loc) %>%
  group_by(time, loc) %>%
  summarise(mean(distance))
