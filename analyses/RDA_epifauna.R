###################################################################################
#                                                                                ##
# RDA of epifauna in Barkley Sound Seagrass for Ecosphere re-submission          ##
# Data are current as of 2018-01-21                                              ##
# Data source: Whippo Thesis Data                                                ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2018-02-02                                                        ##
#                                                                                ##
###################################################################################

# SUMMARY:
# This analysis uses constrained ordination (RDA) analysis to determine 
# relationships between epifaunal community and environmental variables as part of
# a manuscript describing beta diversity of epifauna on multiple scales in Barkley
# Sound, BC seagrass communities. Analysis uses vegan functions and is based on 
# similar analyses conducted by Bostrom et al. 2006.

# Required Files (check that script is loading latest version):
# rawcomm.csv (identified as most-accurate dataset)
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

# 2018-02-02 Updata data source to 'rawcomm.csv'. Started analysis derived from
# Numerical Ecology with R (Legendre)
# 2018-01-21  Script created.

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

# Load packages here
library(ade4)
library(vegan) # distance matrices, model selection
library(MASS)
library(ellipse)
library(FactoMineR)
library(adespatial) # forward selection
library(tidyverse)

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# full community data
epicomm <- read.csv("rawcomm.csv")

# remonve non-epifaunal organisms
epicomm <- epicomm %>%
  select(-Amphipholis.pugetana, -Cockle, -Dinophilus.sp., -Lyonsia.californica, -Nephtys.sp., -Nereis.sp., -Pisaster.ochraceus, -Solaster.sp., -Strongylocentrotus.sp., -Telmessus.cheiragonus, -Nemertea)

# remove redundant columns and summarise species occurrences for middle time period
epicomm_summ <- epicomm %>%
  select(-date, -Sieve, -Time.Code2) %>%
  filter(Time.Code %in% c("B", "C", "C2", "D")) %>%
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
# fix WIC
epicomm_site["9", "sitetime"] <- "WI"

#remove missing species
epicomm_site <- epicomm_site %>%
  select(-Odontosyllis, -Nebalia.sp., -Olivella.sp., -Margarites.helicinus)

# full environmental data
environ <- read.csv("site.info_201801.csv")


###################################################################################
# TRANSFORM COMMUNITY                                                             #
###################################################################################

# convert community to matrix for hellinger transformation as suggested by 
# Legendre & Gallagher 2001
comm_mat <- as.matrix(epicomm_site[,2:32])
# hellinger transform community data for use in RDA
comm_hell <- decostand(comm_mat, method = "hellinger")

###################################################################################
# RDA                                                                             #
###################################################################################

spe.rda <- rda(comm_hell ~., environ)
# '.' calls all variables from env2, default scale = FALSE, scaling = 2
summary(spe.rda)

# how to obtain canonical coefficients from an rda() object
coef(spe.rda)

# retrieval of unadjusted R^2
R2 <- RsquareAdj(spe.rda)$r.squared

# retrieval of adjusted R^2
R2adj <- RsquareAdj(spe.rda)$adj.r.squared







# scaling 2 (default): correlation triplot
plot(spe.rda, main = "Triplot RDA comm_hell ~ environ - scaling 2 - wa scores")
spe2.sc <- scores(spe.rda, choices = 1:2, display = "sp")
arrows(0, 0, spe2.sc[, 1], spe2.sc[, 2], length = 0, lty = 1, col = "red")



############### SUBSECTION HERE

#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#