###################################################################################
#                                                                                ##
# nMDS of epifauna in Barkley Sound Seagrass for Ecosphere re-submission          ##
# Data are current as of 2018-01-24                                              ##
# Data source: Whippo Thesis Data                                                ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2018-01-24                                                        ##
#                                                                                ##
###################################################################################

# SUMMARY:
# This analysis visualizes epifaunal communities in seagrass of Barkley Sound from 
# a single time period in mid-summer. The goal is to determine the amount of 
# community overlap (within vs. across meadows) for similarities.

# Required Files (check that script is loading latest version):
# RawComm_201801.csv (identical to rawcomm-correctedsamplenumbers.csv)

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

# 2018-01-24  Script created.

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

# Load packages here
library(vegan)
library(tidyverse)
library(viridis)
library(ggpubr) # multi-plot visualizations

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
# remove time code from site name to map to environ data
epicomm_summ$sitetime <- substr(epicomm_summ$sitetime, 1, nchar(epicomm_summ$sitetime)-1)




###################################################################################
# nMDS                                                                            #
###################################################################################


# run matrix 
epicomm_mat <- epicomm_summ[,3:48]

# MDS of community
epicomm_mds <- metaMDS(epicomm_mat, distance = "altGower")
epicomm_mds_points <- epicomm_mds$points
epicomm_mds_points <- data.frame(epicomm_mds_points)
plot_data_epi <- data.frame(epicomm_mds_points, epicomm_summ[,1])
vec_sp_epi <- envfit(epicomm_mds$points, epicomm_mat, perm = 1000)
vec_sp_epi_df <- as.data.frame(vec_sp_epi$vectors$arrows*sqrt(vec_sp_epi$vectors$r))
vec_sp_epi_df$species <- rownames(vec_sp_epi_df)
library(plyr)
chulls_epi <- ddply(plot_data_epi, .(sitetime), function(df) df[chull(df$MDS1, df$MDS2), ])
detach(package:plyr)



############### PLOT nMDS

# 2C Taxonomic nMDS
epi_MDS <- ggplot(plot_data_epi, aes(x=MDS1, y=MDS2, color = sitetime)) +  
  theme_minimal() +
  geom_point(size = 4) +
  geom_polygon(data=chulls_epi, aes(x=MDS1, y=MDS2, group=sitetime), fill=NA) 
  #habScale + 
  #rremove("legend") 
epi_MDS

#MDS never reaches a solution and it's an overlapping mess. 

#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#