## Metacom for Ross's data

library(metacom)
library(dplyr)
library(broom)
library(reshape2)
library(Matrix)

data <- read.csv("rawcomm.csv")
#traits <- read.csv("grazertraits.csv")

data5 <- data[(data$site == 'WI'),]
data5 <- rbind(data5, data[(data$site == 'RP'),])
data5 <- rbind(data5, data[(data$site == 'CB'),])
data5 <- rbind(data5, data[(data$site == 'DC'),])
data5 <- rbind(data5, data[(data$site == 'NB'),])


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
## we just need to have this code select the non-zero values of row 10. 

#totals1 <- totals1[1:9,1:35]

#this will automatically trim the species with 0s. hooray!!
totals1[length(totals1[,1]),] <- ifelse(totals1[length(totals1[,1]),]<1, 'zero', 'one')
totals2 <- subset(totals1, select = -grep("zero", totals1[length(totals1[,1]),]))  
totals2 <- totals2[-10,]

## convert to presence absence: 

totals2[totals2 > 0] <- 1
totals2

Metacommunity(totals2)



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
totals.plots2 <- totals.plots2[,-1] # removes column of sites now that rows are named
# calculate sums for each species
totals.plots2[length(totals.plots2[,1])+1,] <- colSums(totals.plots2[,c(1:46)])

#if zeros, then use: 
totals.plots3 <- totals.plots2[,c(order(-as.numeric(totals.plots2[146,])))]

### count individuals in each site
totals.plots3$sum <- cbind(rowSums(totals.plots3[,1:46]))
totals.plots3 <- totals.plots3[(totals.plots3$sum != 0),] # remove any plots with no species

## convert to presence absence: 
tail(totals.plots3)
totals.pl2 <- totals.plots3[-145, -(36:47)] # for all sizes
tail(totals.pl2)
totals.pl2[totals.pl2 > 0] <- 1
#totals.pl2
dim(totals.pl2)

Metacommunity(totals.pl2)

MetaImportance(totals.pl2)



####### TIMEC 5 SITES
## step 1: create a presence absence matrix for the 5 meadows (or for all plots?)
totals5 <- totals[c(1:5),]
totals5[6,] <- colSums(totals5[,c(1:46)])
totals1 <- totals5[,c(order(-totals5[6,]))]
tail(totals1)
totals1 <- totals1[1:5,1:32]


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

totals.plots2 <- totals.plots2[,-1] # remove column of site names now that rows are named

### count individuals in each species
totals.plots2[(length(totals.plots2[,1])+1),] <- colSums(totals.plots2[,c(1:46)])

#if zeros, then use: 
totals.plots3 <- totals.plots2[,c(order(-as.numeric(totals.plots2[length(totals.plots2[,1]),])))]
### count individuals in each site
totals.plots3$sum <- cbind(rowSums(totals.plots3[,1:46]))
totals.plots3 <- totals.plots3[(totals.plots3$sum != 0),]


## convert to presence absence: 
tail(totals.plots3) #identify rows to be removed b/c 
totals.pl2 <- totals.plots3[, -(33:47)] 
tail(totals.pl2) # check to make sure the right cols were removed
totals.pl2 <- totals.pl2[-length(totals.pl2[,1]),] 
totals.pl2[totals.pl2 > 0] <- 1
totals.pl2
dim(totals.pl2)


Metacommunity(totals.pl2)




##### do it for time A
##################################

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

totals$sum <- cbind(rowSums(totals[,1:46]))

totals1 <- totals[,c(order(-totals[6,]))]
tail(totals1)
totals1a <- totals1[1:5,1:37]
tail(totals1a)


## convert to presence absence: 

totals1a[totals1a > 0] <- 1
totals1a

Metacommunity(totals1a)



## Analysis #1: don't collapse samples

timeA5 <- data5[(data5$Time.Code2 == 'A'),]
dim(data)
dim(timeA5)

timeA5$site1 <- paste(timeA5$site, timeA5$Sample)
timeA5 <- timeA5[,-52]
timeA5 <- timeA5[,-(2:5)]

timeA25 <- select(timeA5, c(site1, Idotea.resecata:Alia.carinata))

melted <- melt(timeA25, id.vars = c("site1"))


## get site x sample totals for each species:

totals.plots <- ddply(melted, .(site1, variable), summarise, sum(value))
totals.plots$site1 <- as.factor(totals.plots$site1)
totals.plots2 <- reshape(totals.plots,  direction = "wide", timevar = "variable", idvar = "site1")
rownames(totals.plots2) <- totals.plots2$site1

totals.plots2 <- totals.plots2[,-1] # remove column of site names now that rows are named

### count individuals in each species
totals.plots2[(length(totals.plots2[,1])+1),] <- colSums(totals.plots2[,c(1:46)])

#if zeros, then use: 
totals.plots3 <- totals.plots2[,c(order(-as.numeric(totals.plots2[length(totals.plots2[,1]),])))]
### count individuals in each site
totals.plots3$sum <- cbind(rowSums(totals.plots3[,1:46]))
totals.plots3 <- totals.plots3[(totals.plots3$sum != 0),]


## convert to presence absence: 
tail(totals.plots3) #identify rows to be removed b/c 
totals.pl2 <- totals.plots3[, -(37:47)] 
tail(totals.pl2) # check to make sure the right cols were removed
totals.pl2 <- totals.pl2[-length(totals.pl2[,1]),] 
totals.pl2[totals.pl2 > 0] <- 1
totals.pl2
dim(totals.pl2)


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
tail(totals1)
totals1 <- totals1[1:5,1:43]
tail(totals1)

## convert to presence absence: 

totals1[totals1 > 0] <- 1
totals1

Metacommunity(totals1)





## Analysis #1: don't collapse samples

timeE5 <- data5[(data5$Time.Code2 == 'E'),]
dim(data)
dim(timeE5)

timeE5$site1 <- paste(timeE5$site, timeE5$Sample)
timeE5 <- timeE5[,-52]
timeE5 <- timeE5[,-(2:5)]

timeE25 <- select(timeE5, c(site1, Idotea.resecata:Alia.carinata))

melted <- melt(timeE25, id.vars = c("site1"))


## get site x sample totals for each species:

totals.plots <- ddply(melted, .(site1, variable), summarise, sum(value))
totals.plots$site1 <- as.factor(totals.plots$site1)
totals.plots2 <- reshape(totals.plots,  direction = "wide", timevar = "variable", idvar = "site1")
rownames(totals.plots2) <- totals.plots2$site1

totals.plots2 <- totals.plots2[,-1] # remove column of site names now that rows are named

### count individuals in each species
totals.plots2[(length(totals.plots2[,1])+1),] <- colSums(totals.plots2[,c(1:46)])

#if zeros, then use: 
totals.plots3 <- totals.plots2[,c(order(-as.numeric(totals.plots2[length(totals.plots2[,1]),])))]
### count individuals in each site
totals.plots3$sum <- cbind(rowSums(totals.plots3[,1:46]))
totals.plots3 <- totals.plots3[(totals.plots3$sum != 0),]


## convert to presence absence: 
tail(totals.plots3) #identify rows to be removed b/c 
totals.pl2 <- totals.plots3[, -(44:47)] 
tail(totals.pl2) # check to make sure the right cols were removed
totals.pl2 <- totals.pl2[-length(totals.pl2[,1]),] 
totals.pl2[totals.pl2 > 0] <- 1
totals.pl2
dim(totals.pl2)


Metacommunity(totals.pl2)
