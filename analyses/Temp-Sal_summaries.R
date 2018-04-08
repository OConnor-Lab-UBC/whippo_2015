###################################################################################
#                                                                                ##
# Title of Script                                                                ##
# Data are current as of 2018-02-14                                              ##
# Data source: UBC                                                               ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2018-02-14                                                        ##
#                                                                                ##
###################################################################################

# SUMMARY:


# Required Files (check that script is loading latest version):
# TempSalData.csv

# Associated Scripts:
# none

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

# 2018-02-14 Script created 

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

library(tidyverse)

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

tempsal <- read_csv("./data/TempSalData.csv")

str(tempsal)

# filter out unused sites and make data long
tempsal_used <- tempsal %>%
  filter(Site %in% c("DC", "WI", "BE", "EI", "RP", "NB", "CB", "BI", "CC")) %>%
  select(Date, Time, Site, Station, Temp.s, Temp.2, Temp.b, Sal.s, Sal.2, Sal.b) %>%
  arrange(Site, Date, Station) %>%
  gather(Measurement, Value, -Date, -Time, -Site, -Station) 

# create site_id to extract mean of measurements
tempsal_used <- tempsal_used %>%
  unite(site.id, Date, Time, Site, remove = FALSE)

# remove rows with no data
tempsal_used <- tempsal_used %>%
  subset(Value != "-")

# make dates into date format
tempsal_used$Date <- as.Date(tempsal_used$Date, "%m/%d/%y")

# relabel temp/sal

tempsal_used <- transform(tempsal_used, Unit = substr(Measurement, 1, 3))

# make measurements numeric 
tempsal_used$Value <- as.numeric(tempsal_used$Value)

# separate into appx 30 day sample periods for each site

# DCC = 2012-06-27
DCC_tempsal <- tempsal_used %>%
  filter(Site == "DC") %>%
  filter(Date > "2012-06-12" & Date < "2012-07-12")

# WIC = 2012-06-30
WIC_tempsal <- tempsal_used %>%
  filter(Site == "WI") %>%
  filter(Date > "2012-06-15" & Date < "2012-07-15")

# BEC = 2012-06-15
BEC_tempsal <- tempsal_used %>%
  filter(Site == "BE") %>%
  filter(Date > "2012-05-31" & Date < "2012-06-30")

# EIC = 2012-06-23
EIC_tempsal <- tempsal_used %>%
  filter(Site == "EI") %>%
  filter(Date > "2012-06-08" & Date < "2012-07-08")

# RPC = 2012-07-03
RPC_tempsal <- tempsal_used %>%
  filter(Site == "RP") %>%
  filter(Date > "2012-06-18" & Date < "2012-07-18")

# NB = 2012-07-10
NBC_tempsal <- tempsal_used %>%
  filter(Site == "NB") %>%
  filter(Date > "2012-06-25" & Date < "2012-07-25")

# CBC = 2012-07-05
CBC_tempsal <- tempsal_used %>%
  filter(Site == "CB") %>%
  filter(Date > "2012-06-20" & Date < "2012-07-20")

# BIC = 2012-06-20
BIC_tempsal <- tempsal_used %>%
  filter(Site == "BI") %>%
  filter(Date > "2012-06-05" & Date < "2012-07-05")

# CCC = 2012-07-25
CCC_tempsal <- tempsal_used %>%
  filter(Site == "CC") %>%
  filter(Date > "2012-07-10" & Date < "2012-08-10")

# one big datafile:
temps <- rbind(DCC_tempsal, WIC_tempsal, BEC_tempsal, EIC_tempsal, RPC_tempsal, NBC_tempsal, CBC_tempsal, BIC_tempsal, CCC_tempsal)

#summarize
temps1 <- temps %>%
  group_by(Site) %>%
  filter(Unit == "Tem") %>%
  summarise(., mean(Value))

## are temps and salinities correlated?
mod1 <- lm(temps[(temps$Unit == "Tem"),]$Value ~ temps[(temps$Unit == "Sal"),]$Value)
mod2 <- lm(temps[(temps$Unit == "Tem"),]$Value ~ temps[(temps$Unit == "Sal"),]$Value*temps[(temps$Unit == "Sal"),]$Site)

# TEMP MEANS
BEC_tempsal %>%
  filter(Unit == "Tem") %>%
  summarise(mean(Value))
BIC_tempsal %>%
  filter(Unit == "Tem") %>%
  summarise(mean(Value))
CBC_tempsal %>%
  filter(Unit == "Tem") %>%
  summarise(mean(Value))
CCC_tempsal %>%
  filter(Unit == "Tem") %>%
  summarise(mean(Value))
DCC_tempsal %>%
  filter(Unit == "Tem") %>%
  summarise(mean(Value))
EIC_tempsal %>%
  filter(Unit == "Tem") %>%
  summarise(mean(Value))
NBC_tempsal %>%
  filter(Unit == "Tem") %>%
  summarise(mean(Value))
RPC_tempsal %>%
  filter(Unit == "Tem") %>%
  summarise(mean(Value))
WIC_tempsal %>%
  filter(Unit == "Tem") %>%
  summarise(mean(Value))




# SAL MEANS
BEC_tempsal %>%
  filter(Unit == "Sal") %>%
  summarise(mean(Value))
BIC_tempsal %>%
  filter(Unit == "Sal") %>%
  summarise(mean(Value))
CBC_tempsal %>%
  filter(Unit == "Sal") %>%
  summarise(mean(Value))
CCC_tempsal %>%
  filter(Unit == "Sal") %>%
  summarise(mean(Value))
DCC_tempsal %>%
  filter(Unit == "Sal") %>%
  summarise(mean(Value))
EIC_tempsal %>%
  filter(Unit == "Sal") %>%
  summarise(mean(Value))
NBC_tempsal %>%
  filter(Unit == "Sal") %>%
  summarise(mean(Value))
RPC_tempsal %>%
  filter(Unit == "Sal") %>%
  summarise(mean(Value))
WIC_tempsal %>%
  filter(Unit == "Sal") %>%
  summarise(mean(Value))

###################################################################################
# MANIPULATE DATA                                                                 #
###################################################################################

############### SUBSECTION HERE

#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#


# overall salinity means

overallsal <- tempsal_used %>%
  group_by(Site, Unit) %>%
  summarise(mean(Value))



