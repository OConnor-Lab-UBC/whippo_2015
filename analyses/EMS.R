########################
### Whippo et al eelgrass epifaunal data for Barkley Sound
### elements of metacommunity structure (EMS)
### Code written by Mary O'Connor (with help from Domonik Bahlburg)
### started March 2016
#########################

library(metacom)
library(plyr)
library(dplyr)
library(broom)
library(reshape2)
library(Matrix)
library(lattice)

data <- read.csv("./data/rawcomm.csv")
traits <- read.csv("./data/grazertraits3.csv")
sites <- read.csv("./data/site.info.csv")

## brief visualization:
plot(sites$area ~ sites$dfw)
plot(sites$salinity ~ sites$dfw)
# ok, approximated area not correlated with dfw
#head(traits)
head(data)
dim(data)

### what needs to happen here is: 
# 1a. create a dataframe of species (columns) pooled across sizes for each site. 
# 1b. create a dataframe of species (columns) pooled across sizes for each plot. 
# 2. in this dataframe, include potential gradients (watershed position)

## melt and recast so the datafile goes from long to wide (species as columns)
data.m <- melt(data, id = c(1,2,3,4,5,52))

## some cleaning of species names
levels(data.m$variable)[levels(data.m$variable)== "Bittium.spp."] <- "Lirobittium.spp."
levels(data.m$variable)[levels(data.m$variable)== "Olivella.sp."] <- "Callianax.sp."
levels(data.m$variable)[levels(data.m$variable)== "Cypricercus."] <- "Cyprideis.beaconensis"
levels(data.m$variable)[levels(data.m$variable)== "Odontosyllis"] <- "Polychaete1"

# clean up time code issue
levels(unique(data$Time.Code2))
levels(data.m$Time.Code)
data.m$Time.Code <- as.character(data.m$Time.Code)
data.m$Time.Code[data.m$Time.Code == "C "] <- "C"
data.m$Time.Code <- as.factor(data.m$Time.Code)
#data.m$Sieve <- as.factor(data.m$Sieve)
data.m$value <- as.numeric(data.m$value)
levels(data.m$Time.Code)

## merge with site info
data.s <- merge(data.m, sites, by = "site")

## sum across size classes within plots (samples)

data.p <- ddply(data.s, .(site, Time.Code, Sample, Time.Code2, variable, dfw,order.dfw,area,salinity, shoot.density), summarise, sum(value))

data.p$time.ID <- paste(data.p$site, data.p$Time.Code2, sep = '.') #could look at finer time resolution by using Time.Code here
names(data.p) <- c("site", "Time.Code", "Sample", "Time.Code2", "species", "dfw","order","area","salinity","shoot.density","abundance", "TimeID")

## merge with traits and sort by taxa or functional groups
data.tr <- merge(data.p, traits[,-1], by.x = "species", by.y = "species.names", all.x = TRUE, all.y = FALSE)

## remove all taxa that are not epifauna, because they were not evenly sampled across samples and meadows: 
levels(data.tr$eelgrss.epifauna)
data.e <- data.tr %>% filter(eelgrss.epifauna == c("yes", "sometimes"))
data.y <- data.tr %>% filter(eelgrss.epifauna == "yes")
data.s <- data.tr %>% filter(patch == "yes")
data.g <- data.tr %>% filter(function. == "grazer")
data.tr <- data.y


## from here, create different subsets for different analyses; MOVE TO EMS SUBGROUP FILE IF YOU WANT TO LOOK AT SUBGROUPS. 

## group by sampling times
dataMAY <- data.tr[(data.tr$Time.Code2=="A"),]
dataJULY <- data.tr[(data.tr$Time.Code2=="C" & data.tr$site!="BE" & data.tr$site!="EI" & data.tr$site!="CC" & data.tr$site!="BI"),]
dataAUG <- data.tr[(data.tr$Time.Code2=="E"),]
data3times <- data.tr[(data.tr$site!="BE" & data.tr$site!="EI" & data.tr$site!="CC" & data.tr$site!="BI"),]
dataJULY9 <- data.tr[(data.tr$Time.Code2=="C"),]


## 1. create site-level data by collapsing across plots
start.data <- dataJULY9 # dataMAY, data.mp, dataAUG, dataJULY, data3times
data.ms <- ddply(start.data, .(TimeID, area, species), summarise, sum(abundance)) #order
#data.ms <- data.ms[-nrow(data.ms),]
data2 <- dcast(data.ms, TimeID ~ species, mean) #order


# or collapse across times and plots at sites
#data.pooled <- ddply(start.data, .(site, species), summarise, sum(abundance))
#data2p <- dcast(data.pooled, site ~ species, mean)
#data2 <- data2[-10,-ncol(data2)]
data2 <- data2p

head(data2)
dim(data2)

## order file to impose environmental gradient, if desired
data2 <- data2[order(data2$area),]  #order
data2 <- data2[,-2]

## clean out any empty species or empty sites
rowSums(data2[-1]) -> data2$totals
data2 <- data2[which(data2$totals != "0"),] 
#data3 <- data2[-which(is.na(data2$totals)),] 
rownames(data2) <- data2[,1]
data3 <- data2[, -c(1, ncol(data2))]
data3[data3 > 0] <- 1
cols.to.delete <- which(colSums(data3) == '0')
data4 <- data3[, -(cols.to.delete)]

#remove NaN and Nas. 
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
data4[is.nan.data.frame(data4)] <- 0
data4[is.na(data4)] <- 0
data5 <- data4[,-ncol(data4)]


### run metacommunity analysis 
Metacommunity(data4, verbose = TRUE, order = FALSE) -> meta
meta[2:4]

a <- as.data.frame(meta[1])

pdf('July 9 sites specialists.pdf', width = 7, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="July 9 epifauna YES",  border="black", scales = list(cex = c(0.4, 0.4), x = list(rot = c(90))))
dev.off()

