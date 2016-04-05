## Metacom for Ross's data

library(metacom)
library(plyr)
library(broom)
library(reshape2)
library(Matrix)

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

melted <- melt(timeC2, id.vars = c("site1"))


## get site x sample totals for each species:

totals.plots <- ddply(melted, .(site1, variable), summarise, sum(value))
totals.plots$site1 <- as.factor(totals.plots$site1)
totals.plots2 <- reshape(totals.plots,  direction = "wide", timevar = "variable", idvar = "site1")

rownames(totals.plots2) <- totals.plots2$site1

totals.plots3 <- rbind(totals.plots2, rowSums(totals.plots2[,2:47])) # no zeros, so don't need this line. 
totals.plots3$sum <- cbind(rowSums(totals.plots3[,2:47]))
totals.plots3 <- totals.plots3[(totals.plots3$sum != 0),]

## convert to presence absence: 
totals.pl2 <- totals.plots3[-146, -1] 
totals.pl2 <- totals.pl2[, -47] 
totals.pl2[totals.pl2 > 0] <- 1
totals.pl2
dim(totals.pl2)


test1 <- Coherence(totals.pl2)

Metacommunity(totals.pl2)

MetaImportance(totals.pl2)
