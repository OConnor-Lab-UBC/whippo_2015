###########################################
#physical and community analysis - 2015-05#
###########################################

library(dplyr)
library(tidyr)
library(reshape2)
library(lme4)
library(ggplot2)
library(vegan)

setwd("~/Dropbox (Personal)/R_Scripts/Thesis_Manuscript/Beta_Diversity")

#make sitetime column and sum data per sample w/in site
full.comm <- read.csv("rawcomm_20170103.csv")
new.comm <- full.comm %>%
  select(-date, -Sieve) %>%
  unite(sitetime, site, Time.Code, sep = "", remove = TRUE)
#fix CBC typo
new.comm$sitetime <- gsub("CBC ", "CBC", new.comm$sitetime)
#finish summarization
new.comm <- new.comm %>%
  select(-Time.Code2) %>%
  group_by(sitetime, Sample) %>%
  summarise_each(funs(sum))

# sum total number of orgs per sample
new.comm$plot_tot <- rowSums(new.comm[,3:48])
mean_abun <- new.comm %>%
  group_by(sitetime) %>%
  summarise(mean(plot_tot), sd(plot_tot))

# rename columns
names(mean_abun)[names(mean_abun)=="mean(plot_tot)"] <- "plot_mean"
names(mean_abun)[names(mean_abun)=="sd(plot_tot)"] <- "plot_sd"

# vector of number of samples to calculate SE
samples <- table(new.comm$sitetime)
mean_abun$samples <- samples

# column of SE
mean_abun$SE <- mean_abun$plot_sd/sqrt(mean_abun$samples)

# export as excel
write.csv(mean_abun, "20170408_abundance.csv")

