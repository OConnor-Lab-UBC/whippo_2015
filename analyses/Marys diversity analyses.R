### Mary's attempt to look at richness metrics
### April 12 2016
### for Whippo et al, seagrass biodiversity paper

library(vegan)
library(BiodiversityR)
library(plyr)
library(reshape2)
library(dplyr)
library(MuMIn)

data <- read.csv("./data/rawcomm.csv")
#data_old <- read.csv("./data/plot_data_copy.csv")
traits <- read.csv("./data/grazertraits3.csv")
sites <- read.csv("./data/site.info.csv")

traits <- traits[,-c(3,8:10)]

## brief visualization:
plot(log(sites$area) ~ sites$dfw)
plot(sites$salinity ~ sites$dfw)

plot(sites$epiphytes~ sites$fetch.est)
plot(sites$shoot.density~ sites$fetch.est)
plot(sites$fetch.est~ sites$dfw)
plot(log(sites$area) ~ sites$fetch.est)

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

data.p <- ddply(data.s, .(site, Time.Code, Sample, Time.Code2, variable, dfw,order.dfw,area,salinity, shoot.density, fetch.est), summarise, sum(value))

data.p$time.ID <- paste(data.p$site, data.p$Time.Code2, sep = '.') #could look at finer time resolution by using Time.Code here
names(data.p) <- c("site", "Time.Code", "Sample", "Time.Code2", "species", "dfw","order","area","salinity","shoot.density","fetch","abundance", "TimeID")


## merge with traits and sort by taxa or functional groups
data.tr <- merge(data.p, traits[,-1], by.x = "species", by.y = "species.names", all.x = TRUE, all.y = FALSE)

## remove all taxa that are not epifauna, because they were not evenly sampled across samples and meadows: 
levels(data.tr$eelgrss.epifauna)
data.e <- data.tr %>% filter(eelgrss.epifauna == c("yes", "sometimes"))
data.y <- data.tr %>% filter(eelgrss.epifauna == "yes")
data.tr <- data.e

## group by sampling times
dataMAY <- data.tr[(data.tr$Time.Code2=="A"),]
dataJULY <- data.tr[(data.tr$Time.Code2=="C" & data.tr$site!="BE" & data.tr$site!="EI" & data.tr$site!="CC" & data.tr$site!="BI"),]
dataAUG <- data.tr[(data.tr$Time.Code2=="E"),]
data3times <- data.tr[(data.tr$site!="BE" & data.tr$site!="EI" & data.tr$site!="CC" & data.tr$site!="BI"),]
dataJULY9 <- data.tr[(data.tr$Time.Code2=="C"),]

## for each plot, estimate relative abundance of grazers
## first remove filter feeders and predators:
data.t <- data3times # dataJULY9
data7 <- subset(data.t, function. != "unknown", select = c(1:17))
data7 <- subset(data.t, group != "echinoderm", select = c(1:17))
#data7 <- subset(data7, function. != "predator", select = c(1:11))
#data7 <- subset(data7, function. != "unknown", select = c(1:11))
#[data7$species!='Caprella.spp.',]
data8 <- ddply(data.t, .(site, Sample, order, dfw, area, salinity, fetch, function., Time.Code2), summarise, sum(abundance)) #
data9 <- dcast(data8, site + Sample + fetch + area + order + dfw + Time.Code2~ function., sum) #order
data9$total <- (data9$detritovore + data9$suspension.feeder + data9$grazer + data9$predator + data9$filter.feeder) #+ data9$unknown 
data9$pgrazer <- data9$grazer/data9$total
data9$pdet <- data9$detritovore/data9$total

par(mfrow=(c(2,2)))
plot((data9$pgrazer) ~ data9$fetch, pch = 19, xlab = 'Fetch', ylab = 'grazers / total', main = 'abundance / 0.28m2') #, ylim = c(0,1)
plot(log(data9$total) ~ data9$fetch, pch = 19, xlab = 'Fetch', ylab = 'ln(total inverts)', main = 'abundance / 0.28m2')
plot(log(data9$grazer+1) ~ data9$fetch, pch = 19, xlab = 'Fetch', ylab = 'ln(grazers+1)', main = 'abundance / 0.28m2')
plot(log(data9$filter.feeder + 1) ~ data9$fetch, pch = 19, xlab = 'Distance from Freshwater', ylab = 'ln(filter feeders)', main = 'abundance / 0.28m2')

plot(I(log(data9$grazer/(data9$filter.feeder+data9$suspension.feeder))) ~ data9$fetch, pch = 19, xlab = 'Fetch', ylab = 'ln(grazers/filter feeders)', main = 'abundance / 0.28m2')


## some stats for july9sites
hist(log(data9$total))
mod1a <- lm(log(data9$total+1) ~ data9$fetch)
mod2a <- lm(log(data9$grazer+1) ~ data9$fetch)

mod1 <- lm(log(data9$total+1) ~ log(data9$area) * data9$fetch)
mod2 <- lm(log(data9$grazer+1) ~ log(data9$area) * data9$fetch)

mod1b <- lm(log(data9$total+1) ~ data9$dfw * data9$fetch)
mod2b <- lm(log(data9$grazer+1) ~ data9$dfw * data9$fetch)

mod1c <- lm(log(data9$total+1) ~ log(data9$area) + data9$dfw * data9$fetch)
mod2c <- lm(log(data9$grazer+1) ~ log(data9$area) + data9$dfw * data9$fetch)

mod1d <- lm(log(data9$total+1) ~ log(data9$area))
mod2d <- lm(log(data9$grazer+1) ~ log(data9$area))

mod1e <- lm(log(data9$total+1) ~ log(data9$area) + data9$fetch)
mod2e <- lm(log(data9$grazer+1) ~ log(data9$area) + data9$fetch)

mod1f <- lm(log(data9$total+1) ~ data9$dfw + data9$fetch)
mod2f <- lm(log(data9$grazer+1) ~ data9$dfw + data9$fetch)

model.sel(mod1a, mod1b, mod1, mod1c, mod1e, mod1f)
model.sel(mod2a, mod2b, mod2, mod2c, mod2e, mod2f)

## some stats for resample sites
hist(log(data9$total))
plot(log(data9$total + 1) ~ data9$Time.Code2)
plot(log(data9$total + 1) ~ log(data9$area))
plot(log(data9$total + 1) ~ (data9$dfw))
plot(log(data9$total + 1) ~ (data9$fetch))

mod1a <- lm(log(data9$total+1) ~ data9$fetch)
mod2a <- lm(log(data9$grazer+1) ~ data9$fetch)

mod1 <- lm(log(data9$total+1) ~ log(data9$area) * data9$fetch)
mod2 <- lm(log(data9$grazer+1) ~ log(data9$area) * data9$fetch)

mod1b <- lm(log(data9$total+1) ~ data9$dfw * data9$fetch)
mod2b <- lm(log(data9$grazer+1) ~ data9$dfw * data9$fetch)

mod1c <- lm(log(data9$total+1) ~ log(data9$area) + data9$dfw * data9$fetch)
mod2c <- lm(log(data9$grazer+1) ~ log(data9$area) + data9$dfw * data9$fetch)

mod1d <- lm(log(data9$total+1) ~ log(data9$area))
mod2d <- lm(log(data9$grazer+1) ~ log(data9$area))

mod1e <- lm(log(data9$total+1) ~ log(data9$area) + data9$fetch)
mod2e <- lm(log(data9$grazer+1) ~ log(data9$area) + data9$fetch)

mod1f <- lm(log(data9$total+1) ~ data9$dfw + data9$fetch)
mod2f <- lm(log(data9$grazer+1) ~ data9$dfw + data9$fetch)

mod1g <- lm(log(data9$total+1) ~ log(data9$area) + data9$dfw * data9$fetch + data9$Time.Code2)
mod2g <- lm(log(data9$grazer+1) ~ log(data9$area) + data9$dfw * data9$fetch + data9$Time.Code2)

mod1h <- lm(log(data9$total+1) ~ log(data9$area) + data9$dfw * data9$fetch*data9$Time.Code2)
mod2h <- lm(log(data9$grazer+1) ~ log(data9$area) + data9$dfw * data9$fetch*data9$Time.Code2)

mod1i <- lm(log(data9$total+1) ~ log(data9$area) + data9$dfw + data9$fetch*data9$Time.Code2)
mod2i <- lm(log(data9$grazer+1) ~ log(data9$area) + data9$dfw + data9$fetch*data9$Time.Code2)

mod1j <- lm(log(data9$total+1) ~ log(data9$area) + data9$fetch + data9$dfw * data9$Time.Code2)
mod2j <- lm(log(data9$grazer+1) ~ log(data9$area) + data9$fetch + data9$dfw * data9$Time.Code2)

mod1k <- lm(log(data9$total+1) ~ log(data9$area)* data9$Time.Code2*data9$fetch + data9$dfw)

model.sel(mod1a, mod1b, mod1, mod1c, mod1e, mod1f, mod1g, mod1h, mod1i, mod1j, mod1k)
model.sel(mod2a, mod2b, mod2, mod2c, mod2e, mod2f, mod2g, mod2h, mod2i, mod2j)

## by now these are grossly overfit models. we have 15 populations (5 meadows x 3 times) so if time is a factor we get one other thing. fetch? area? dfw? not sure how to do this. 

# order is not contributing, but area is. try dfw, and then all sites and add time, and then I think we're done.

boxplot(log(dataJULY9tots[,7]) ~ dataJULY9tots$dfw, pch = 19, xlab = 'Distance from Freshwater', ylab = 'total inverts', main = 'raw comm data')
boxplot(log(data_old[(data_old$Time == '2.5'),3]) ~ data_old[(data_old$Time == '2.5'),10], pch = 19, xlab = 'Distance from Freshwater', ylab = 'total inverts', main = 'plot data')

mod1 <- lm(log(dataJULY9tots[,7]+0.01) ~ dataJULY9tots$dfw, na.action = na.omit)



#### some statistics
mod1 < lm()
