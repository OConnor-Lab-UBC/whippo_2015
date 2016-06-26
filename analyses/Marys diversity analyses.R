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

data$date1 <- mdy(data$date)

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
data.m <- melt(data, id = c(1,2,3,4,5,52, 53))

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

data.p <- ddply(data.s, .(site, date1, Sample, Time.Code2, variable, dfw,order.dfw,area,salinity, shoot.density, fetch.est), summarise, sum(value))

data.p$time.ID <- paste(data.p$site, data.p$Time.Code2, sep = '.') #could look at finer time resolution by using Time.Code here
names(data.p) <- c("site", "Date", "Sample", "Time.Code2", "species", "dfw","order","area","salinity","shoot.density","fetch","abundance", "TimeID")


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
data8 <- ddply(data.t, .(site, Sample, order, dfw, area, salinity, fetch, function., Date, Time.Code2), summarise, sum(abundance)) #
data9 <- dcast(data8, site + Sample + fetch + area + order + dfw + Date + Time.Code2 ~ function., sum) #order
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
library(RColorBrewer)
cols <- brewer.pal(n = 9, name="Set1")
cols_dfw <- cols[round(data9$dfw)*0.5]
plot(log(data9$total + 1) ~ data9$Date, pch = 19, col = cols_dfw)
plot(log(data9$total + 1) ~ log(data9$area), pch = 19, col = cols_dfw)
plot(log(data9$total + 1) ~ (data9$dfw), pch = 19, col = data9$Time.Code2)
plot(log(data9$grazer + 1) ~ (data9$dfw), pch = 19, col = data9$Time.Code2)
plot(log(data9$total + 1) ~ (data9$fetch))

mod1a <- lm(log(data9$total+1) ~ data9$fetch*data9$Date)
mod2a <- lm(log(data9$grazer+1) ~ data9$fetch*data9$Date)

mod1 <- lm(log(data9$total+1) ~ log(data9$area) *data9$Date)
mod2 <- lm(log(data9$grazer+1) ~ log(data9$area) *data9$Date)

mod1b <- lm(log(data9$total+1) ~ data9$dfw * data9$Date)
mod2b <- lm(log(data9$grazer+1) ~ data9$dfw * data9$Date)

mod1c <- lm(log(data9$total+1) ~ data9$dfw * data9$fetch)
mod2c <- lm(log(data9$grazer+1) ~ data9$dfw * data9$fetch)

model.sel(mod1a, mod1b, mod1, mod1c)
model.sel(mod2a, mod2b, mod2, mod2c)

## by now these are grossly overfit models. we have 15 populations (5 meadows x 3 times) so if time is a factor we get one other thing. fetch? area? dfw? not sure how to do this. 
### diversity analyses
div.data <- dcast(data.t[,c(1:4,12)], site + Date + Sample ~ species, sum)

H <- diversity(div.data[,-(c(1:3))], index ="shannon")
S <- diversity(div.data[,-(c(1:3))], index ="simpson")

div.summary <- cbind(div.data[c(1:3)], H, S)
div.summary2 <- merge(div.summary, sites, by.x = c("site"), by.y = c("site"))

plot(div.summary2$H ~ div.summary2$Date, pch = 19, col = div.summary2$dfw)
plot(div.summary2$H ~ div.summary2$area, pch = 19, col = div.summary2$dfw)
plot(div.summary2$H ~ div.summary2$fetch.est, pch = 19, col = div.summary2$dfw)
plot(div.summary2$H ~ div.summary2$dfw, pch = 19, col = div.summary2$Date)
plot(div.summary2$H ~ div.summary2$area)

mod1a <- lm(div.summary2$H ~ div.summary2$fetch*div.summary2$Date)
mod2a <- lm(div.summary2$S ~ div.summary2$fetch*div.summary2$Date)

mod1 <- lm(div.summary2$H ~ log(div.summary2$area) * div.summary2$Date)
mod2 <- lm(div.summary2$S  ~ log(div.summary2$area) * div.summary2$Date)

mod1b <- lm(div.summary2$H ~ div.summary2$dfw * div.summary2$Date)
mod2b <- lm(div.summary2$S ~ div.summary2$dfw * div.summary2$Date)

mod1c <- lm(div.summary2$H ~ div.summary2$dfw * div.summary2$fetch)
mod2c <- lm(div.summary2$S ~ div.summary2$dfw * div.summary2$fetch)

model.sel(mod1a, mod1b, mod1, mod1c)
model.sel(mod2a, mod2b, mod2, mod2c)