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

tempsal <- read_csv("TempSalData.csv")

str(tempsal)

# filter out unused sites and make data long
tempsal_used <- tempsal %>%
  filter(Site %in% c("DC", "WI", "BC", "EI", "RP", "NB", "CB", "BI", "CC")) %>%
  select(Date, Time, Site, Station, Temp.s, Temp.2, Temp.b, Sal.s, Sal.2, Sal.b) %>%
  arrange(Site, Date, Station) %>%
  gather(Measurement, Value, -Date, -Time, -Site, -Station) 

# create site_id to extract mean of measurements
tempsal_used <- tempsal_used %>%
  unite(site.id, Date, Time, Site, remove = FALSE)

# remove rows with no data
tempsal_used <- tempsal_used %>%
  subset(Value != "-")

###################################################################################
# MANIPULATE DATA                                                                 #
###################################################################################

############### SUBSECTION HERE

#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#






stocks <- tibble(
  time = as.Date('2009-01-01') + 0:9,
  X = rnorm(10, 0, 1),
  Y = rnorm(10, 0, 2),
  Z = rnorm(10, 0, 4)
)

gather(stocks, stock, price, -time)
stocks %>% gather(stock, price, -time)