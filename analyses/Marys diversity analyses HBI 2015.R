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

R <- specnumber(div.data[,-c(1:2)])
H <- diversity(div.data[,-c(1:2)], index ="shannon")
S <- diversity(div.data[,-(c(1:2))], index ="simpson")

div.summary <- cbind(div.data[c(1:2)], H, S, R)
#div.summary2 <- merge(div.summary, sites, by.x = c("site"), by.y = c("site"))

jpg('HBI 2015 diversity.jpg', width = 7, height = 9)
par(mfrow = c(3,1))
plot(div.summary$R ~ div.summary$site, pch = 19, ylim = c(0,20), axes = FALSE, xlab = 'Site', ylab = 'Observed Species Richness', main = 'Hakai epifauna 2015')
axis(1, at = c(1,2,3,4,5,6,7,8,9), lab = c(unique(levels(div.summary$site))), pos = c(0,-1), cex.lab = 0.8, las = 2)
axis(2, c(0,5,10,15,20), pos = c(0,0), cex.lab = 0.8, las = 2)

plot(div.summary$H ~ div.summary$site, pch = 19, ylim = c(0,3), axes = FALSE, xlab = 'Site', ylab = 'Diversity (Shannon)', main = 'Hakai epifauna 2015')
axis(1, at = c(1,2,3,4,5,6,7,8,9), lab = c(unique(levels(div.summary$site))), pos = c(0,-1), cex.lab = 0.8, las = 2)
axis(2, c(0,1,2,3), pos = c(0,0), cex.lab = 0.8, las = 2)

plot(div.summary$S ~ div.summary$site, pch = 19, ylim = c(0,1), axes = FALSE, xlab = 'Site', ylab = 'Evenness (Simpson)')
axis(1, at = c(1,2,3,4,5,6,7,8,9), lab = c(unique(levels(div.summary$site))), pos = c(0,-1), cex.lab = 0.8, las = 2)
axis(2, c(0,0.5,1), pos = c(0,0), cex.lab = 0.8, las = 2)

dev.off()

mod1 <- lm(div.summary$R ~ div.summary$site)
mod2 <- lm(div.summary$H ~ div.summary$site)

#### do it for 2014

data14 <- read.csv("./data/Hakai2014Data.csv", sep = ",")

## melt and recast so the datafile goes from long to wide (species as columns)
data14 <- data14[,-1]
data14.m <- melt(data14, id = c(1,2,3))
head(data14.m)

## some cleaning of site names
levels(data14.m$Site)[levels(data14.m$Site)== "McMullin North"] <- "MM.N"
levels(data14.m$Site)[levels(data14.m$Site)== "McMullin"] <- "MM.S"
levels(data14.m$Site)[levels(data14.m$Site)== "Lower"] <- "Lowr"
levels(data14.m$Site)[levels(data14.m$Site)== "Goose"] <- "Gse.W"
levels(data14.m$Site)[levels(data14.m$Site)== "Goose East"] <- "Gse.E"
levels(data14.m$Site)[levels(data14.m$Site)== "Triquet"] <- "Trq.N"
levels(data14.m$Site)[levels(data14.m$Site)== "Triquet/No Name Cove"] <- "Trq.S"
levels(data14.m$Site)[levels(data14.m$Site)== "Choked"] <- "Chk"
levels(data14.m$Site)[levels(data14.m$Site)== "Lower Choked"] <- "Chk.L"
levels(data14.m$Site)[levels(data14.m$Site)== ""] <- "1"

## merge with site info
#data.s <- merge(data.m, sites, by = "site")

## sum across size classes within plots (samples)
data14.p <- ddply(data14.m, .(Site, Sample.number, variable), summarise, sum(value))
names(data14.p) <- c("site", "Sample", "species","abundance")

### diversity analyses
## 1. create site-level data by collapsing across plots
start.data <- data14.p 
data.ms <- ddply(start.data, .(site, species, Sample), summarise, sum(abundance, na.rm = TRUE)) #order

div.data <- dcast(data.ms, site + Sample ~ species, mean)
#div.data2 <- div.data[div.data!= 'NA'] 


H <- diversity(div.data[,-c(1:2)], index ="shannon")
S <- diversity(div.data[,-(c(1:2))], index ="simpson")

div.summary <- cbind(div.data[c(1:2)], H, S)
#div.summary2 <- merge(div.summary, sites, by.x = c("site"), by.y = c("site"))

pdf('HBI 2014.pdf', width = 7, height = 9)
par(mfrow = c(2,1))
plot(div.summary$H ~ div.summary$site, pch = 19, ylim = c(0,3), axes = FALSE, xlab = 'Site', ylab = 'Diversity (Shannon)', main = 'Hakai epifauna 2014')
axis(1, at = c(1,2,3,4,5,6,7,8,9,10), lab = c(unique(levels(div.summary$site))), pos = c(0,-1), cex.lab = 0.8, las = 2)
axis(2, c(0,1,2,3), pos = c(0,0), cex.lab = 0.8, las = 2)

plot(div.summary$S ~ div.summary$site, pch = 19, ylim = c(0,1), axes = FALSE, xlab = 'Site', ylab = 'Evenness (Simpson)')
axis(1, at = c(1,2,3,4,5,6,7,8,9,10), lab = c(unique(levels(div.summary$site))), pos = c(0,-1), cex.lab = 0.8, las = 2)
axis(2, c(0,0.5,1), pos = c(0,0), cex.lab = 0.8, las = 2)

dev.off()
