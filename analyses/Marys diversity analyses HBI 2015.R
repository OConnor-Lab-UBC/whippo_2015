### Mary's attempt to look at richness metrics
### April 12 2016
### adapting Whippo et al code for HBI 2015 data

library(vegan)
library(BiodiversityR)
library(plyr)
library(reshape2)
library(dplyr)
library(MuMIn)

data <- read.csv("./data/Hakai2015Data.csv", sep = ",")
#data_old <- read.csv("./data/plot_data_copy.csv")
#traits <- read.csv("./data/grazertraits3.csv")
#sites <- read.csv("./data/site.info.csv")

#traits <- traits[,-c(3,8:10)]

#data$date1 <- mdy(data$date)

## brief visualization:
#plot(log(sites$area) ~ sites$dfw)
#plot(sites$salinity ~ sites$dfw)

#plot(sites$epiphytes~ sites$fetch.est)
#plot(sites$shoot.density~ sites$fetch.est)
#plot(sites$fetch.est~ sites$dfw)
#plot(log(sites$area) ~ sites$fetch.est)

## melt and recast so the datafile goes from long to wide (species as columns)
data <- data[,-1]
data.m <- melt(data, id = c(1,2,3))
head(data.m)

## some cleaning of site names
levels(data.m$site)[levels(data.m$site)== "McMullins N"] <- "MM.N"
levels(data.m$site)[levels(data.m$site)== "McMullins S"] <- "MM.S"
levels(data.m$site)[levels(data.m$site)== "Goose N"] <- "Gse.N"
levels(data.m$site)[levels(data.m$site)== "Goose W"] <- "Gse.W"
levels(data.m$site)[levels(data.m$site)== "Goose E"] <- "Gse.E"
levels(data.m$site)[levels(data.m$site)== "Triquet N"] <- "Trq.N"
levels(data.m$site)[levels(data.m$site)== "Triquet S"] <- "Trq.S"
levels(data.m$site)[levels(data.m$site)== "Choked Pass, S. Pigu"] <- "Chk.SP"
levels(data.m$site)[levels(data.m$site)== "Sandspit"] <- "Chk.S"

## merge with site info
#data.s <- merge(data.m, sites, by = "site")

## sum across size classes within plots (samples)
data.p <- ddply(data.m, .(site, Sample, variable), summarise, sum(value))
names(data.p) <- c("site", "Sample", "species","abundance")

## merge with traits and sort by taxa or functional groups
#data.tr <- merge(data.p, traits[,-1], by.x = "species", by.y = "species.names", all.x = TRUE, all.y = FALSE)

## remove all taxa that are not epifauna, because they were not evenly sampled across samples and meadows: 
#levels(data.tr$eelgrss.epifauna)
#data.e <- data.tr %>% filter(eelgrss.epifauna == c("yes", "sometimes"))
#data.y <- data.tr %>% filter(eelgrss.epifauna == "yes")
#data.s <- data.tr %>% filter(patch == "yes")
#data.g <- data.tr %>% filter(function. == "grazer")
#data.tr <- data.y

### diversity analyses
## 1. create site-level data by collapsing across plots
start.data <- data.p 
data.ms <- ddply(start.data, .(site, species, Sample), summarise, sum(abundance, na.rm = TRUE)) #order
#data.ms <- data.ms[-nrow(data.ms),]
#data2 <- dcast(data.ms, site + Sample ~ species, mean) #order

head(data2)
dim(data2)

div.data <- dcast(data.ms, site + Sample ~ species, mean)
#div.data2 <- div.data[div.data!= 'NA'] 


H <- diversity(div.data[,-c(1:2)], index ="shannon")
S <- diversity(div.data[,-(c(1:2))], index ="simpson")

div.summary <- cbind(div.data[c(1:2)], H, S)
#div.summary2 <- merge(div.summary, sites, by.x = c("site"), by.y = c("site"))

par(mfrow = c(2,1))
plot(div.summary$H ~ div.summary$site, pch = 19, ylim = c(0,3), xlab = 'Site', ylab = 'Diversity (Shannon)')
plot(div.summary$S ~ div.summary$site, pch = 19, ylim = c(0,1), xlab = 'Site', ylab = 'Evenness (Simpson)')

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