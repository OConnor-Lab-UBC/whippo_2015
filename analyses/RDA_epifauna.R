###################################################################################
#                                                                                ##
# RDA of epifauna in Barkley Sound Seagrass for Ecosphere re-submission          ##
# Data are current as of 2018-02-20                                              ##
# Data source: Whippo Thesis Data                                                ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2018-02-20                                                        ##
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

# 2018-02-20 Added epicomm_201802.csv to replace rawcomm.csv for all community analysis and altered specifics of scripts as needed.
# 2018-02-02 Updated data source to 'rawcomm.csv'. Started analysis derived from
# Numerical Ecology with R (Legendre)
# 2018-01-21  Script created.

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

# Load packages here
library(ade4)
library(vegan) # distance matrices, model selection
library(ellipse)
library(FactoMineR)
library(adespatial) # forward selection
library(tidyverse)
library(ggvegan) # vegan plots in ggplot2 framework
library(viridis) # plotting palette

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# full community data
epicomm <- read.csv("./data/epicomm_201802.csv")

# remove redundant columns and summarise species occurrences for middle time period
epicomm_summ <- epicomm %>%
  select(-date, -Sieve) %>%
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


#remove missing species
epicomm_site <- epicomm_site %>%
  select(-Nebalia.sp., -Margarites.helicinus)

# full environmental data
environ_full <- read.csv("./data/site.info_201801.csv")


###################################################################################
# TRANSFORM COMMUNITY & ENVIRONMENTAL VARIABLES                                   #
###################################################################################

# convert community to matrix for hellinger transformation as suggested by 
# Legendre & Gallagher 2001
comm_mat <- as.matrix(epicomm_site[,2:29])
# hellinger transform community data for use in RDA
comm_hell <- decostand(comm_mat, method = "hellinger")
comm_hell <- as.data.frame(comm_hell)
row.names(comm_hell) <- epicomm_site$sitetime

# reduce variable list to ecologically relevant factors
environ <- environ_full %>%
  select(-site, -dfAI, -area, -dfw)
# center and scale environmental variables
environ_scaled <- transform(environ, area.ha = scale(area.ha), salinity = scale(salinity), shoot.density = scale(shoot.density), epiphytes = scale(epiphytes), fetch.meters = scale(fetch.meters))
row.names(environ_scaled) <- epicomm_site$sitetime



###################################################################################
# RDA                                                                             #
###################################################################################

spe_rda <- rda(comm_hell ~., environ_scaled)
# '.' calls all variables from env2, default scale = FALSE, scaling = 2
summary(spe_rda)

# how to obtain canonical coefficients from an rda() object
coef(spe_rda)

# retrieval of unadjusted R^2
R2 <- RsquareAdj(spe_rda)$r.squared

# retrieval of adjusted R^2
R2adj <- RsquareAdj(spe_rda)$adj.r.squared




###################################################################################
# TRIPLOT VISUALIZATION                                                           #
###################################################################################


# scaling 2 (default): correlation triplot
# plot(spe_rda, main = "Triplot RDA Epifaunal Community ~ Environmental Variables")
#  scaling 2 - wa scores
# 700 x 650
# spe2_sc <- scores(spe_rda, choices = 1:2, display = "sp")
# arrows(0, 0, spe2_sc[, 1], spe2_sc[, 2], length = 0, lty = 1, col = "red")


# attempt to use ggvegan to plot as ggplot

# transform vegan RDA output into ggplot readable format
gg_rda <- fortify(spe_rda)

# remove unused factors
gg_rda_min <- gg_rda %>%
  filter(Score != 'constraints')

# add size vector
txt_size <-c(rep(4.99,28), rep(5,9), rep(5.01,5))

# run plot
ggplot(data = gg_rda_min, aes(RDA1, RDA2)) +
  geom_text(aes(label=Label, color = Score), position = position_jitter(width = 0.1, height = 0.1) ) +
  scale_color_viridis(end = 0.85, discrete = TRUE, option = "viridis") +
  theme_minimal() + 
  theme(legend.position='none') +
  xlim(-1, 1) +
  ylim(-1, 1) +
  geom_segment(x = 0, y = 0, xend = gg_rda_min$RDA1, yend = gg_rda_min$RDA2, aes(color = Score))

###################################################################################
# PERMUTATION TESTS OF RDA                                                        #
###################################################################################

# global test of the RDA result
anova.cca(spe_rda, step = 1000)

# tests of all canonical axes
anova.cca(spe_rda, by = "axis", step = 1000)

# apply Kaiser-Guttman criterion to residual axes
spe_rda$CA$eig[spe_rda$CA$eig > mean(spe_rda$CA$eig)]

###################################################################################
# VARIATION INFLATION FACTORS (VIF)                                               #
###################################################################################

# test collinearity of variables
vif.cca(spe_rda)
# values over 10 should be examined, values over 20 are strongly collinear

###################################################################################
# FORWARD SELECTION                                                               #
###################################################################################

# to reduce or eliminate variables with high collinearity, forward selection can be
# used to select an appropriate model.

# using a double stopping criterion (Blanchet et al. 2008)
# RDA with all explanatory variables
spe_rda_all <- rda(comm_hell ~., environ_scaled)

# global adjusted R^2
(R2a_all <- RsquareAdj(spe_rda_all)$adj.r.squared)

# forward selection using adespatial's forward.sel()
library(adespatial)
forward.sel(comm_hell, environ_scaled, adjR2thresh = R2a_all)
# analysis stops with variable that exceeds R2 threshold


############### SUBSECTION HERE

#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#