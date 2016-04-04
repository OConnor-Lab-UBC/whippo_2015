## Metacom for Ross's data

library(metacom)
library(dplyr)
library(broom)

data <- read.csv("rawcomm.csv")
traits <- read.csv("grazertraits.csv")

## step 1: create a presence absence matrix for the 9 meadows (or for all plots?)

timeC <- data[(data$Time.Code2 == 'C'),]
dim(data)
dim(timeC)

timeC <- timeC[,-52]
timeC <- timeC[,-(2:5)]

## Analysis #1: collapse samples within sites
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

Metacommunity(totals1)

## Analysis #1: don't collapse samples

timeC <- data[(data$Time.Code2 == 'C'),]
dim(data)
dim(timeC)

timeC$site1 <- paste(timeC$site, timeC$Sample)
timeC <- timeC[,-52]
timeC <- timeC[,-(2:5)]

timeC2 <- select(timeC, c(site1, Idotea.resecata:Alia.carinata))


## get site x sample totals for each species:
### got stuck here trying to access each site x sample. maybe just number each sample and try that, losing site identity for now?

site.samples <- length(unique(timeC2$site1))
site.sam <- timeC2$site1
totals.plots <- colSums(timeC[timeC2$site1 == site.sam[i],2:47])
totals.plots <- 
  for (i in length(site.sam)) {
  rbind(totals.plots, colSums(timeC2[timeC2$site1 == site.sam[i],2:47]))
}

totals.plots <- rbind(totals.plots, colSums(timeC2[timeC2$site1 == 'RP',2:47]))
totals.plots <- rbind(totals.plots, colSums(timeC2[timeC2$site1 == 'WI',2:47]))
totals.plots <- rbind(totals.plots, colSums(timeC2[timeC2$site1 == 'NB',2:47]))
totals.plots <- rbind(totals.plots, colSums(timeC2[timeC2$site1 == 'CB',2:47]))
totals.plots <- rbind(totals.plots, colSums(timeC2[timeC2$site1 == 'EI',2:47]))
totals.plots <- rbind(totals.plots, colSums(timeC2[timeC2$site1 == 'CC',2:47]))
totals.plots <- rbind(totals.plots, colSums(timeC2[timeC2$site1 == 'BI',2:47]))
totals.plots <- rbind(totals.plots, colSums(timeC2[timeC2$site1 == 'BE',2:47]))
totals.plots <- as.data.frame(totals.plots)

totals.plots <- rbind(totals.plots, colSums(timeC2[,2:47]))
rownames(totals.plots)<- c('DC', 'RP', 'WI', 'NB', 'CB', 'EI', 'CC', 'BI', 'BE', 'Tot')

totals.pl1 <- totals.plots[,c(order(-totals.plots[10,]))]
totals.pl1 <- totals.pl1[1:9,1:35]

## convert to presence absence: 

totals.pl1[totals.pl1 > 0] <- 1
totals.pl1

test1 <- Coherence(totals1)

Metacommunity(totals1)