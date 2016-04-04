## Metacom for Ross's data

library(metacom)
library(dplyr)
library(broom)

data <- read.csv("rawcomm.csv")
traits <- read.csv("grazertraits.csv")


env.data2 <- as.data.frame(c('DC', 'RP', 'WI', 'NB', 'CB'))
colnames(env.data2) <- 'site'

## step 1: create a presence absence matrix for the 9 meadows (or for all plots?)

timeC <- data[(data$Time.Code2 == 'C'),]
dim(data)
dim(timeC)

timeC <- timeC[,-52]
timeC <- timeC[,-(2:5)]

## get site totals for each species and size class: 
totals <- colSums(timeC[timeC$site == 'DC',2:47])
totals <- rbind(totals, colSums(timeC[timeC$site == 'RP',2:47]))
totals <- rbind(totals, colSums(timeC[timeC$site == 'WI',2:47]))
totals <- rbind(totals, colSums(timeC[timeC$site == 'NB',2:47]))
totals <- rbind(totals, colSums(timeC[timeC$site == 'CB',2:47]))
totals <- rbind(totals, colSums(timeC[timeC$site == 'EI',2:47]))
totals <- rbind(totals, colSums(timeC[timeC$site == 'CC',2:47]))
totals <- rbind(totals, colSums(timeC[timeC$site == 'BI',2:47]))
totals <- rbind(totals, colSums(timeC[timeC$site == 'BE',2:47]))
totals <- as.data.frame(totals)

totals <- rbind(totals, colSums(timeC[,2:47]))
rownames(totals)<- c('DC', 'RP', 'WI', 'NB', 'CB', 'EI', 'CC', 'BI', 'BE', 'Tot')

totals1 <- totals[,c(order(-totals[10,]))]
totals1 <- totals1[1:9,1:35]

## convert to presence absence: 

totals1[totals1 > 0] <- 1
totals1

test1 <- Coherence(totals1)

