## Metacom for Ross's data

library(metacom)
library(plyr)
library(dplyr)
library(broom)
library(reshape2)
library(Matrix)
library(lattice)

data <- read.csv("rawcomm.csv")
traits <- read.csv("grazertraits.csv")
head(traits)
head(data)


#### CREATING DATASETS WE MIGHT USE

## All 5 resurveyed sites.
data5 <- data[(data$site == 'WI'),]
data5 <- rbind(data5, data[(data$site == 'RP'),])
data5 <- rbind(data5, data[(data$site == 'CB'),])
data5 <- rbind(data5, data[(data$site == 'DC'),])
data5 <- rbind(data5, data[(data$site == 'NB'),])

# exclude the smallest individuals: not a good idea, you lose 10 species this way, not just juveniles
#data.s <- data[(data$Sieve > '2'),]
#data <- data.s

## TIME C, 9 MEADOWS, SITES  
timeC <- data[(data$Time.Code2 == 'C'),]
dim(timeC)
timeC <- timeC[,-52] # get rid of columns that are not needed for this analysis
timeC <- timeC[,-(2:5)]

## collapse samples within sites: get site totals for each species and size class: 
datafile <- timeC ## name of datafile you're using #length(datafile[1,])
dim(datafile)
totals <- colSums(datafile[datafile$site == 'DC',2:47])
totals <- rbind(totals, colSums(datafile[datafile$site == 'RP',2:47]))
totals <- rbind(totals, colSums(datafile[datafile$site == 'WI',2:47]))
totals <- rbind(totals, colSums(datafile[datafile$site == 'NB',2:47]))
totals <- rbind(totals, colSums(datafile[datafile$site == 'CB',2:47]))
totals <- rbind(totals, colSums(datafile[datafile$site == 'EI',2:47]))
totals <- rbind(totals, colSums(datafile[datafile$site == 'CC',2:47]))
totals <- rbind(totals, colSums(datafile[datafile$site == 'BI',2:47]))
totals <- rbind(totals, colSums(datafile[datafile$site == 'BE',2:47]))
totals <- as.data.frame(totals)

totals <- rbind(totals, colSums(timeC[,2:47]))
rownames(totals)<- c('DC', 'RP', 'WI', 'NB', 'CB', 'EI', 'CC', 'BI', 'BE', 'Tot')

#Trim the species with 0 occurrances in this dataset
totals[length(totals[,1]),] <- ifelse(totals[length(totals[,1]),]<1, 'zero', 'one')
totals2 <- subset(totals, select = -grep("zero", totals[length(totals[,1]),])) 
## get rid of row of species sums
totals2 <- (totals2[-10,])

## convert to presence absence: 
totals2[totals2 > 0] <- 1
totals2

Metacommunity(totals2) -> meta
meta
#MetaImportance(totals2, margin = 2)

a <- as.data.frame(meta[1])

pdf('9sitesCallsizes.pdf', width = 7, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="Community 9 sites C",  border="black")
dev.off()


## Analysis #2: don't collapse samples

## create dataset of interest, here subsetting time C from the full data.
timeC <- data[(data$Time.Code2 == 'C'),]
dim(data)
dim(timeC)

## collapse site and sample into a single plot-level identifier
timeC$site1 <- paste(timeC$site, timeC$Sample)
timeC <- timeC[,-52]
timeC <- timeC[,-(2:5)]

## isolate identifier column 
timeC2 <- select(timeC, c(site1, Idotea.resecata:Alia.carinata))

melted <- melt(timeC2, id.vars = c("site1"))

## get site x sample totals for each species:
totals.plots <- ddply(melted, .(site1, variable), summarise, sum(value))
totals.plots$site1 <- as.factor(totals.plots$site1)
totals.plots2 <- reshape(totals.plots,  direction = "wide", timevar = "variable", idvar = "site1")

rownames(totals.plots2) <- totals.plots2$site1
totals.plots2 <- totals.plots2[,-1] # removes column of sites, now that rows are named

# calculate sums for each species
totals.plots2[length(totals.plots2[,1])+1,] <- colSums(totals.plots2[,c(1:length(totals.plots2[1,]))])
#calculate sums for each site
totals.plots2$sum <- cbind(rowSums(totals.plots2)) 

#this will automatically trim the species with 0s. hooray!!
# remove any plots with no species
totals.plots2 <- totals.plots2[(totals.plots2$sum != 0),] 
totals.plots2[length(totals.plots2[,1]),] <- ifelse(totals.plots2[length(totals.plots2[,1]),]<1, 'zero', 'one')
totals.plots3 <- subset(totals.plots2, select = -grep("zero", totals.plots2[length(totals.plots2[,1]),]))

## convert to presence absence: 
tail(totals.plots3)
totals.pl2 <- totals.plots3[-length(totals.plots3[,1]), -length(totals.plots3[1,])]
tail(totals.pl2)
totals.pl2[totals.pl2 > 0] <- 1
dim(totals.pl2)

#Metacommunity(totals.pl2)

Metacommunity(totals.pl2) -> meta
meta
#MetaImportance(totals.pl2, margin = 2)

a <- as.data.frame(meta[1])

pdf('9plotsCallsizes.pdf', width = 65, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="Community 9 sites C",  border="black")
dev.off()



####### 5 SITES
## step 1: create a presence absence matrix for the 5 meadows (or for all plots?)
timeC <- data5[(data5$Time.Code2 == 'C'),]
dim(data5)
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
totals <- as.data.frame(totals)

totals <- rbind(totals, colSums(timeC[,2:47]))
rownames(totals)<- c('DC', 'RP', 'WI', 'NB', 'CB', 'Tot')

#calculate sums for each site
totals$sum <- cbind(rowSums(totals))

#this will automatically trim the species with 0s. hooray!!
totals[length(totals[,1]),] <- ifelse(totals[length(totals[,1]),]<1, 'zero', 'one')
totals2 <- subset(totals, select = -grep("zero", totals[length(totals[,1]),]))  
## get rid of row and col of species sums
totals2 <- totals2[-length(totals2[,1]), -length(totals2[1,])]

## convert to presence absence: 
totals2[totals2 > 0] <- 1

Metacommunity(totals2) -> meta
#MetaImportance(totals2, margin = 2)

a <- as.data.frame(meta[1])

pdf('9sitesCallsizes.pdf', width = 7, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="Community 9 sites C",  border="black")
dev.off()


## Analysis #2: ALL PLOTS, 5 SITES

timeC5 <- data5[(data5$Time.Code2 == 'C'),]
dim(data)
dim(timeC5)

## collapse site and sample into a single plot-level identifier
timeC5$site1 <- paste(timeC5$site, timeC5$Sample)
timeC5 <- timeC5[,-52]
timeC5 <- timeC5[,-(2:5)]

## isolate identifier column 
timeC25 <- select(timeC5, c(site1, Idotea.resecata:Alia.carinata))

## get site x sample totals for each species:
melted <- melt(timeC25, id.vars = c("site1"))
totals.plots <- ddply(melted, .(site1, variable), summarise, sum(value))
totals.plots$site1 <- as.factor(totals.plots$site1)
totals.plots2 <- reshape(totals.plots,  direction = "wide", timevar = "variable", idvar = "site1")
rownames(totals.plots2) <- totals.plots2$site1

totals.plots2 <- totals.plots2[,-1] # remove column of site names now that rows are named

# calculate sums for each species
totals.plots2[length(totals.plots2[,1])+1,] <- colSums(totals.plots2[,c(1:length(totals.plots2[1,]))])
#calculate sums for each site
totals.plots2$sum <- cbind(rowSums(totals.plots2)) 

#this will automatically trim the species with 0s. hooray!!
# remove any plots with no species
totals.plots2 <- totals.plots2[(totals.plots2$sum != 0),] 
totals.plots2[length(totals.plots2[,1]),] <- ifelse(totals.plots2[length(totals.plots2[,1]),]<1, 'zero', 'one')
totals.plots3 <- subset(totals.plots2, select = -grep("zero", totals.plots2[length(totals.plots2[,1]),]))
totals.plots3 <- totals.plots3[-length(totals.plots3[,1]), -length(totals.plots3[1,])]

## convert to presence absence: 
totals.plots3[totals.plots3 > 0] <- 1
Metacommunity(totals.plots3) -> meta

a <- as.data.frame(meta[1])

pdf('9plotsCallsizes.pdf', width = 65, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="Community 9 sites C",  border="black")
dev.off()




##################################
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

#calculate sums for each site
totals$sum <- cbind(rowSums(totals))

#this will automatically trim the species with 0s. hooray!!
totals[length(totals[,1]),] <- ifelse(totals[length(totals[,1]),]<1, 'zero', 'one')
totals2 <- subset(totals, select = -grep("zero", totals[length(totals[,1]),]))  
## get rid of row and col of species sums
totals2 <- totals2[-length(totals2[,1]), -length(totals2[1,])]

## convert to presence absence: 
totals2[totals2 > 0] <- 1

Metacommunity(totals2) -> meta

a <- as.data.frame(meta[1])

pdf('9plotsCallsizes.pdf', width = 65, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="Community 9 sites C",  border="black")
dev.off()


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
totals.plots2 <- totals.plots2[,-1] # removes column of sites, now that rows are named

# calculate sums for each species
totals.plots2[length(totals.plots2[,1])+1,] <- colSums(totals.plots2[,c(1:length(totals.plots2[1,]))])
#calculate sums for each site
totals.plots2$sum <- cbind(rowSums(totals.plots2)) 

#this will automatically trim the species with 0s. hooray!!
# remove any plots with no species
totals.plots2 <- totals.plots2[(totals.plots2$sum != 0),] 
totals.plots2[length(totals.plots2[,1]),] <- ifelse(totals.plots2[length(totals.plots2[,1]),]<1, 'zero', 'one')
totals.plots3 <- subset(totals.plots2, select = -grep("zero", totals.plots2[length(totals.plots2[,1]),]))

## convert to presence absence: 
tail(totals.plots3)
totals.pl2 <- totals.plots3[-length(totals.plots3[,1]), -length(totals.plots3[1,])]
tail(totals.pl2)
totals.pl2[totals.pl2 > 0] <- 1
dim(totals.pl2)

Metacommunity(totals.pl2) -> meta
meta
#MetaImportance(totals.pl2, margin = 2)

a <- as.data.frame(meta[1])

pdf('5plotsAallsizes.pdf', width = 65, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="Community 9 sites C",  border="black")
dev.off()





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

#calculate sums for each site
totals$sum <- cbind(rowSums(totals))

#this will automatically trim the species with 0s. hooray!!
totals[length(totals[,1]),] <- ifelse(totals[length(totals[,1]),]<1, 'zero', 'one')
totals2 <- subset(totals, select = -grep("zero", totals[length(totals[,1]),]))  
## get rid of row and col of species sums
totals2 <- totals2[-length(totals2[,1]), -length(totals2[1,])]

## convert to presence absence: 
totals2[totals2 > 0] <- 1

Metacommunity(totals2) -> meta

a <- as.data.frame(meta[1])

pdf('9plotsCallsizes.pdf', width = 65, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="Community 9 sites C",  border="black")
dev.off()




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
totals.plots2 <- totals.plots2[,-1] # removes column of sites, now that rows are named

# calculate sums for each species
totals.plots2[length(totals.plots2[,1])+1,] <- colSums(totals.plots2[,c(1:length(totals.plots2[1,]))])
#calculate sums for each site
totals.plots2$sum <- cbind(rowSums(totals.plots2)) 

#this will automatically trim the species with 0s. hooray!!
# remove any plots with no species
totals.plots2 <- totals.plots2[(totals.plots2$sum != 0),] 
totals.plots2[length(totals.plots2[,1]),] <- ifelse(totals.plots2[length(totals.plots2[,1]),]<1, 'zero', 'one')
totals.plots3 <- subset(totals.plots2, select = -grep("zero", totals.plots2[length(totals.plots2[,1]),]))

## convert to presence absence: 
tail(totals.plots3)
totals.pl2 <- totals.plots3[-length(totals.plots3[,1]), -length(totals.plots3[1,])]
tail(totals.pl2)
totals.pl2[totals.pl2 > 0] <- 1
dim(totals.pl2)

Metacommunity(totals.pl2) -> meta
meta

a <- as.data.frame(meta[1])
pdf('5plotsEallsizes.pdf', width = 65, height = 20)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="Community 5 all plots E",  border="black")
dev.off()