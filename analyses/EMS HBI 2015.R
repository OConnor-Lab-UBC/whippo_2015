########################
### Whippo et al eelgrass epifaunal data for Barkley Sound
### elements of metacommunity structure (EMS)
### Code written by Mary O'Connor (with help from Domonik Bahlburg)
### adapting Whippo code to Hakai Code, June 2016
#########################

library(metacom)
library(plyr)
library(dplyr)
#library(broom)
library(reshape2)
library(Matrix)
library(lattice)

data <- read.csv("./data/Hakai2015Data.csv", sep = ",")
traits <- read.csv("./data/grazertraits3.csv")
sites <- read.csv("./data/site.info.csv")

dim(data)
head(data)

## melt and recast so the datafile goes from long to wide (species as columns)
data <- data[,-1]
data.m <- melt(data, id = c(1,2,3))
head(data.m)

## some cleaning of species names
#levels(data.m$variable)[levels(data.m$variable)== "Bittium.spp."] <- "Lirobittium.spp."
#levels(data.m$variable)[levels(data.m$variable)== "Olivella.sp."] <- "Callianax.sp."
#levels(data.m$variable)[levels(data.m$variable)== "Cypricercus."] <- "Cyprideis.beaconensis"
#levels(data.m$variable)[levels(data.m$variable)== "Odontosyllis"] <- "Polychaete1"

## merge with site info
#data.s <- merge(data.m, sites, by = "site")

## sum across size classes within plots (samples)
data.p <- ddply(data.m, .(site, Sample, variable), summarise, sum(value))
names(data.p) <- c("site", "Sample", "species","abundance")
#data.p1 <- data.p[-which(is.na(data.p$totals)),] 


## merge with traits and sort by taxa or functional groups
#data.tr <- merge(data.p, traits[,-1], by.x = "species", by.y = "species.names", all.x = TRUE, all.y = FALSE)

## remove all taxa that are not epifauna, because they were not evenly sampled across samples and meadows: 
#levels(data.tr$eelgrss.epifauna)
#data.e <- data.tr %>% filter(eelgrss.epifauna == c("yes", "sometimes"))
#data.y <- data.tr %>% filter(eelgrss.epifauna == "yes")
#data.s <- data.tr %>% filter(patch == "yes")
#data.g <- data.tr %>% filter(function. == "grazer")
#data.tr <- data.y

## 1. create site-level data by collapsing across plots
start.data <- data.p # dataMAY, data.mp, dataAUG, dataJULY, data3times
data.ms <- ddply(start.data, .(site, species), summarise, sum(abundance, na.rm = TRUE)) #order
#data.ms <- data.ms[-nrow(data.ms),]
data2 <- dcast(data.ms, site ~ species, mean) #order

head(data2)
dim(data2)

## order file to impose environmental gradient, if desired
#data2 <- data2[order(data2$area),]  #order
#data2 <- data2[,-2]

## clean out any empty species or empty sites
rowSums(data2[-1]) -> data2$totals
data2 <- data2[which(data2$totals != "0"),] 
#data3 <- data2[-which(is.na(data2$totals)),] 
rownames(data2) <- data2[,1]
data3 <- data2[, -c(1, ncol(data2))]
#data3[data3 > 0] <- 1
cols.to.delete <- which(colSums(data3) == '0')
data4 <- data3[, -(cols.to.delete)]

#remove NaN and Nas, not needed now.
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
data2[is.nan.data.frame(data2)] <- 0
data2[is.na(data2)] <- 0
data2 <- data2[,-ncol(data2)]


### run metacommunity analysis 
Metacommunity(data4, verbose = TRUE) -> meta
meta[2:4]

a <- as.data.frame(meta[1])

pdf('HBI 2015.pdf', width = 7, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="Hakai 2015",  border="black", scales = list(cex = c(0.4, 0.4), x = list(rot = c(90))))
dev.off()

