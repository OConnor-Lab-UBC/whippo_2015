## Metacom for Ross's data

library(metacom)
library(dplyr)
library(broom)
library(reshape2)
library(Matrix)

data <- read.csv("rawcomm.csv")
traits <- read.csv("grazertraits.csv")

data5 <- data[c((data$site == 'RB'),(data$site == 'WI'),(data$site == 'CB'),(data$site == 'NB'),(data$site == 'DC')),]

# exclude the smallest individuals
data.s <- data[(data$Sieve > '2'),]
data <- data.s

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
# for big stuff only: totals1 <- totals1[1:9, 1:25]

## convert to presence absence: 

totals1[totals1 > 0] <- 1
totals1

test1 <- Coherence(totals1)

Metacommunity(totals1)

## Analysis #2: don't collapse samples

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
totals.plots2 <- totals.plots2[,-1]
totals.plots2[146,] <- colSums(totals.plots2[,c(1:46)])

# x <- colSums(totals.plots2[,c(2:47)]) 

#if zeros, then use: 
totals.plots3 <- totals.plots2[,c(order(-as.numeric(totals.plots2[146,])))]
totals.plots3$sum <- cbind(rowSums(totals.plots3[,2:46]))
totals.plots3 <- totals.plots3[(totals.plots3$sum != 0),]

## convert to presence absence: 
totals.pl2 <- totals.plots3[-145, 1:35] # for all sizes
totals.pl2[totals.pl2 > 0] <- 1
#totals.pl2
dim(totals.pl2)


test1 <- Coherence(totals.pl2)

Metacommunity(totals.pl2)

MetaImportance(totals.pl2)



####### TIMEC 5 SITES
## step 1: create a presence absence matrix for the 5 meadows (or for all plots?)
totals5 <- totals[c(1:5,10),]
totals1 <- totals5[,c(order(-totals5[6,]))]
totals1 <- totals1[1:5,1:35]

## convert to presence absence: 

totals1[totals1 > 0] <- 1
totals1

Metacommunity(totals1)




## Analysis #2: don't collapse samples

timeC5 <- data5[(data5$Time.Code2 == 'C'),]
dim(data)
dim(timeC5)

timeC5$site1 <- paste(timeC5$site, timeC5$Sample)
timeC5 <- timeC5[,-52]
timeC5 <- timeC5[,-(2:5)]

timeC25 <- select(timeC5, c(site1, Idotea.resecata:Alia.carinata))

melted <- melt(timeC25, id.vars = c("site1"))


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




##### do it for time A

## step 1: create a presence absence matrix for the 9 meadows (or for all plots?)

timeA <- data5[(data5$Time.Code2 == 'A'),]
dim(data5)
dim(timeA)

timeA <- timeA[,-52]
timeA <- timeA[,-(2:5)]

## Analysis #1: collapse samples within sites
## get site totals for each species and size class: 
totals <- colSums(timeA[timeA$site == 'DC',2:47])
totals <- rbind(totals, colSums(timeA[timeA$site == 'RP',2:47]))
totals <- rbind(totals, colSums(timeA[timeA$site == 'WI',2:47]))
totals <- rbind(totals, colSums(timeA[timeA$site == 'NB',2:47]))
totals <- rbind(totals, colSums(timeA[timeA$site == 'CB',2:47]))
totals <- as.data.frame(totals)

totals <- rbind(totals, colSums(timeA[,2:47]))
rownames(totals)<- c('DC', 'RP', 'WI', 'NB', 'CB', 'Tot')

totals1 <- totals[,c(order(-totals[6,]))]
totals1 <- totals1[1:5,1:36]

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





##### do it for time E

## step 1: create a presence absence matrix for the 9 meadows (or for all plots?)

timeE <- data5[(data5$Time.Code2 == 'E'),]
dim(data5)
dim(timeE)

timeE <- timeE[,-52]
timeE <- timeE[,-(2:5)]

## Analysis #1: collapse samples within sites
## get site totals for each species and size class: 
totals <- colSums(timeE[timeE$site == 'DC',2:47])
totals <- rbind(totals, colSums(timeE[timeE$site == 'RP',2:47]))
totals <- rbind(totals, colSums(timeE[timeE$site == 'WI',2:47]))
totals <- rbind(totals, colSums(timeE[timeE$site == 'NB',2:47]))
totals <- rbind(totals, colSums(timeE[timeE$site == 'CB',2:47]))
totals <- as.data.frame(totals)

totals <- rbind(totals, colSums(timeE[,2:47]))
rownames(totals)<- c('DC', 'RP', 'WI', 'NB', 'CB', 'Tot')

totals1 <- totals[,c(order(-totals[6,]))]
totals1 <- totals1[1:5,1:36]