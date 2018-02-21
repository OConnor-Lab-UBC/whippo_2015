###################################################################################
#                                                                                ##
# Basic Diversity Metrics                                                        ##
# Data are current as of 2018-02-21                                              ##
# Data source: O'Connor Lab - UBC                                                ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2018-02-21                                                        ##
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
# NONE

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

# 2018-02-21 Created script. 

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

library(tidyverse)
library(vegan)

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# full community data
epicomm <- read_csv("./data/epicomm_201802.csv")

# remove redundant columns and summarise species occurrences
epicomm_summ <- epicomm %>%
  select(-date, -Sieve) %>%
  unite(sitetime, site, 'Time Code', sep = "", remove = TRUE)
# fix CBC typo
epicomm_summ$sitetime <- gsub("CBC ", "CBC", epicomm_summ$sitetime)
# finish summarization for every sample by summing
epicomm_summ <- epicomm_summ %>%
  group_by(sitetime, Sample) %>%
  summarise_all(funs(sum))
# summarise entire sites with total spp abund per site
epicomm_site <- epicomm_summ %>%
  select(-Sample) %>%
  group_by(sitetime) %>%
  summarise_all(funs(sum))

###################################################################################
# DIVERSITY MEASURES                                                              #
###################################################################################

############### RARIFIED RICHNESS




#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#