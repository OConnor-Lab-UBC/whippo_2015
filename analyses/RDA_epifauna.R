###################################################################################
#                                                                                ##
# RDA of epifauna in Barkley Sound Seagrass for Ecosphere re-submission          ##
# Data are current as of 2018-01-21                                              ##
# Data source: Whippo Thesis Data                                                ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2018-01-21                                                        ##
#                                                                                ##
###################################################################################

# SUMMARY:
# This analysis uses constrained ordination (RDA) analysis to determine 
# relationships between epifaunal community and environmental variables as part of
# a manuscript describing beta diversity of epifauna on multiple scales in Barkley
# Sound, BC seagrass communities. Analysis uses vegan functions and is based on 
# similar analyses conducted by Bostrom et al. 2006.

# Required Files (check that script is loading latest version):
# RawComm_201801.csv (identical to rawcomm-correctedsamplenumbers.csv)
# site.info_201801.csv (updated 2018-01-22 with correct values)

# TO DO

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# RECENT CHANGES TO SCRIPT                                                        #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# SUMMARY STATS                                                                   #    
#                                                                                 #
###################################################################################

###################################################################################
# RECENT CHANGES TO SCRIPT                                                        #
###################################################################################

# 2018-01-21  Script created.

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

# Load packages here
library(vegan)
library(tidyverse)
library(MVN)

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# full community data
epicomm <- read.csv("RawComm_201801.csv")

# remove redundant columns and summarise species occurrences for middle time period
epicomm_summ <- epicomm %>%
  select(-date, -Sieve) %>%
  filter(Time.Code %in% c("B", "C", "D")) %>%
  unite(sitetime, site, Time.Code, sep = "", remove = TRUE)
# finish summarization for every sample by summing
epicomm_summ <- epicomm_summ %>%
  group_by(sitetime, Sample) %>%
  summarise_all(funs(sum))
# summarise entire sites with total spp abund per site
epicomm_site <- epicomm_summ %>%
  select(-Sample) %>%
  group_by(sitetime) %>%
  summarise_all(funs(sum))
# remove time code from site name to map to environ data
epicomm_site$sitetime <- substr(epicomm_site$sitetime, 1, nchar(epicomm_site$sitetime)-1)

# full environmental data
environ <- read.csv("site.info_201801.csv")


###################################################################################
# SUMMARY STATS                                                                   #
###################################################################################

# convert community to matrix for hellinger transformation as suggested by 
# Legendre & Gallagher 2001
comm_mat <- as.matrix(epicomm_site[,2:47])
# hellinger transform community data for use in RDA
comm_hell <- decostand(comm_mat, method = "hellinger")



############### SUBSECTION HERE

#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#