### Mary's attempt to look at dominant species
### September 20, 2015
### for Whippo et al, seagrass biodiversity paper

library(vegan)
library(BiodiversityR)
library(plyr)


## first merge the datafile with the species traits file

data <- read.csv("rawcomm.csv")
traits <- read.csv("grazertraits.csv")

## get site totals for each species and size class: 
totals <- colSums(data[data$site == 'DC',6:51])
totals <- rbind(totals, colSums(data[data$site == 'RP',6:51]))
totals <- rbind(totals, colSums(data[data$site == 'WI',6:51]))
totals <- rbind(totals, colSums(data[data$site == 'NB',6:51]))
totals <- rbind(totals, colSums(data[data$site == 'CB',6:51]))
rownames(totals)<- c('DC', 'RP', 'WI', 'NB', 'CB')
totals <- as.data.frame(totals)
env.data2 <- as.data.frame(c('DC', 'RP', 'WI', 'NB', 'CB'))
colnames(env.data2) <- 'site'

## this works!
DCrad <- rankabundance(totals, y = env.data2, factor = "site", level = 'DC')
RPrad <- rankabundance(totals, y = env.data2, factor = "site", level = 'RP')
WIrad <- rankabundance(totals, y = env.data2, factor = "site", level = 'WI')
NBrad <- rankabundance(totals, y = env.data2, factor = "site", level = 'NB')
CBrad <- rankabundance(totals, y = env.data2, factor = "site", level = 'CB')

par(mfrow=c(2, 2), omi = c(.5, .5, .5, .5))
rankabunplot(DCrad, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "DC")
rankabunplot(RPrad, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "RP")
rankabunplot(NBrad, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "NB")
rankabunplot(CBrad, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "CB")

rankabunplot(WIrad, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "WI")





## get site totals for each species and size class, for all sites time C: 
totalsC <- colSums(data[((data$site == 'DC')&(data$Time.Code2 == "C")),6:51])
totalsC <- rbind(totalsC, colSums(data[((data$site == 'RP')&(data$Time.Code2 == "C")),6:51]))
totalsC <- rbind(totalsC, colSums(data[((data$site == 'WI')&(data$Time.Code2 == "C")),6:51]))
totalsC <- rbind(totalsC, colSums(data[((data$site == 'NB')&(data$Time.Code2 == "C")),6:51]))
totalsC <- rbind(totalsC, colSums(data[((data$site == 'CC')&(data$Time.Code2 == "C")),6:51]))
totalsC <- rbind(totalsC, colSums(data[((data$site == 'BE')&(data$Time.Code2 == "C")),6:51]))
totalsC <- rbind(totalsC, colSums(data[((data$site == 'BI')&(data$Time.Code2 == "C")),6:51]))
totalsC <- rbind(totalsC, colSums(data[((data$site == 'CB')&(data$Time.Code2 == "C")),6:51]))
totalsC <- rbind(totalsC, colSums(data[((data$site == 'EI')&(data$Time.Code2 == "C")),6:51]))
rownames(totalsC)<- c('DC', 'RP', 'WI', 'NB', 'CC', 'BE', 'BI', 'CB', 'EI')
totalsC <- as.data.frame(totalsC, rownames = c('DC', 'RP', 'WI', 'NB', 'CC', 'BE', 'BI', 'CB', 'EI'))

env.data <- as.data.frame(c('DC', 'RP', 'WI', 'NB', 'CC', 'BE', 'BI', 'CB', 'EI'))
colnames(env.data) <- 'site'

## this works!
DCradc <- rankabundance(totalsC, y = env.data, factor = "site", level = 'DC')
RPradc <- rankabundance(totalsC, y = env.data, factor = "site", level = 'RP')
WIradc <- rankabundance(totalsC, y = env.data, factor = "site", level = 'WI')
NBradc <- rankabundance(totalsC, y = env.data, factor = "site", level = 'NB')
CCradc <- rankabundance(totalsC, y = env.data, factor = "site", level = 'CC')
BEradc <- rankabundance(totalsC, y = env.data, factor = "site", level = 'BE')
BIradc <- rankabundance(totalsC, y = env.data, factor = "site", level = 'BI')
CBradc <- rankabundance(totalsC, y = env.data, factor = "site", level = 'CB')
EIradc <- rankabundance(totalsC, y = env.data, factor = "site", level = 'EI')

par(mfrow=c(2, 2), omi = c(.4, .4, .4, .4))
rankabunplot(DCradc, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "DC - C")
rankabunplot(RPradc, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "RP - C")
rankabunplot(EIradc, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "EI - C")
rankabunplot(BEradc, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "BE - C")

par(mfrow=c(2, 2), omi = c(.5, .5, .5, .5))
rankabunplot(BIradc, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "BI - C")
rankabunplot(CCradc, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "CC - C")
rankabunplot(NBradc, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "NB - C")
rankabunplot(CBradc, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "CB - C")

rankabunplot(WIradc, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000),main = "WI - C")





## get site totals for each species and size class, for all sites time A: 
totalsA <- colSums(data[((data$site == 'DC')&(data$Time.Code2 == "A")),6:51])
totalsA <- rbind(totalsA, colSums(data[((data$site == 'RP')&(data$Time.Code2 == "A")),6:51]))
totalsA <- rbind(totalsA, colSums(data[((data$site == 'WI')&(data$Time.Code2 == "A")),6:51]))
totalsA <- rbind(totalsA, colSums(data[((data$site == 'NB')&(data$Time.Code2 == "A")),6:51]))
totalsA <- rbind(totalsA, colSums(data[((data$site == 'CB')&(data$Time.Code2 == "A")),6:51]))
rownames(totalsA)<- c('DC', 'RP', 'WI', 'NB', 'CB')
totalsA <- as.data.frame(totalsA, rownames = c('DC', 'RP', 'WI', 'NB', 'CB'))

DCradA <- rankabundance(totalsA, y = env.data2, factor = "site", level = 'DC')
RPradA <- rankabundance(totalsA, y = env.data2, factor = "site", level = 'RP')
WIradA <- rankabundance(totalsA, y = env.data2, factor = "site", level = 'WI')
NBradA <- rankabundance(totalsA, y = env.data2, factor = "site", level = 'NB')
CBradA <- rankabundance(totalsA, y = env.data2, factor = "site", level = 'CB')

par(mfrow=c(2, 2), omi = c(.5, .5, .5, .5))
rankabunplot(DCradA, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "DCA")
rankabunplot(RPradA, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "RPA")
rankabunplot(NBradA, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "NBA")
rankabunplot(CBradA, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "CBA")

rankabunplot(WIradA, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "WIA")


## get site totals for each species and size class, for all sites time E: 
totalsE <- colSums(data[((data$site == 'DC')&(data$Time.Code2 == "E")),6:51])
totalsE <- rbind(totalsE, colSums(data[((data$site == 'RP')&(data$Time.Code2 == "E")),6:51]))
totalsE <- rbind(totalsE, colSums(data[((data$site == 'WI')&(data$Time.Code2 == "E")),6:51]))
totalsE <- rbind(totalsE, colSums(data[((data$site == 'NB')&(data$Time.Code2 == "E")),6:51]))
totalsE <- rbind(totalsE, colSums(data[((data$site == 'CB')&(data$Time.Code2 == "E")),6:51]))
rownames(totalsE)<- c('DC', 'RP', 'WI', 'NB', 'CB')
totalsE <- as.data.frame(totalsE, rownames = c('DC', 'RP', 'WI', 'NB', 'CB'))

DCradE <- rankabundance(totalsE, y = env.data2, factor = "site", level = 'DC')
RPradE <- rankabundance(totalsE, y = env.data2, factor = "site", level = 'RP')
WIradE <- rankabundance(totalsE, y = env.data2, factor = "site", level = 'WI')
NBradE <- rankabundance(totalsE, y = env.data2, factor = "site", level = 'NB')
CBradE <- rankabundance(totalsE, y = env.data2, factor = "site", level = 'CB')

par(mfrow=c(2, 2), omi = c(.5, .5, .5, .5))
rankabunplot(DCradE, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "DCE")
rankabunplot(RPradE, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "RPE")
rankabunplot(NBradE, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "NBE")
rankabunplot(CBradE, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "CBE")

rankabunplot(WIradE, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "WIE")

### Wizard all times:
par(mfrow=c(2, 2), omi = c(.5, .5, .5, .5))

rankabunplot(WIradA, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "WIA")

rankabunplot(WIradc, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000),main = "WI - C")

rankabunplot(WIradE, scale = "logabun", specnames=c(1:46), xlim = c(0, 46), ylim = c(0, 12000), main = "WIE")
