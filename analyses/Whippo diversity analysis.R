########################
### Epifaunal diversity patterns among eelgrass meadows suggest metacommunity structure
### R. Whippo, N. Knight, C. Prentice, J. Cristiani, M. Siegle, M. I. O'Connor
### Analysis: elements of metacommunity structure (EMS)
### Code written by Mary O'Connor and Ross Whippo
### 2016-2017
### Provides figures 2 and univariate diversity statistics for Tables 1-3, S5
#########################


# load libraries and data -------------------------------------------------

library(vegan)
library(BiodiversityR)
library(plyr)
library(reshape2)
library(dplyr)
library(MuMIn)
library(lubridate)
library(ggplot2)
library(purrr)

data <- read.csv("./data/rawcomm.csv")
traits <- read.csv("./data/grazertraits3.csv")
sites <- read.csv("./data/site.info.csv")

# delete, add and redefine columns ----------------------------------------
traits <- traits[,-c(3,8:10)]

data$date1 <- mdy(data$date)

## brief visualization:
plot(log(sites$area) ~ sites$dfw)
plot(sites$salinity ~ sites$dfw)

plot(sites$epiphytes~ sites$fetch.jc)
plot(sites$fetch.est ~ sites$fetch.jc)
plot(sites$shoot.density~ sites$fetch.jc)
plot(sites$fetch.jc~ sites$dfw)
plot(log(sites$area) ~ sites$fetch.jc)

### what needs to happen here is: 
# 1a. create a dataframe of species (columns) pooled across sizes for each site. 
# 1b. create a dataframe of species (columns) pooled across sizes for each plot. 
# 2. in this dataframe, include potential gradients (watershed position)


# data preparation ---------------------------------------------------
## Melt and recast so the datafile goes from long to wide (species as columns)
data.m <- melt(data, id = c(1,2,3,4,5,52, 53))

## Clean and correct species names
levels(data.m$variable)[levels(data.m$variable)== "Bittium.spp."] <- "Lirobittium.spp."
levels(data.m$variable)[levels(data.m$variable)== "Olivella.sp."] <- "Callianax.sp."
levels(data.m$variable)[levels(data.m$variable)== "Cypricercus."] <- "Cyprideis.beaconensis"
levels(data.m$variable)[levels(data.m$variable)== "Odontosyllis"] <- "Polychaete1"

# Clean up time code issue
levels(unique(data$Time.Code2))
levels(data.m$Time.Code)
data.m$Time.Code <- as.character(data.m$Time.Code)
data.m$Time.Code[data.m$Time.Code == "C "] <- "C"
data.m$Time.Code <- as.factor(data.m$Time.Code)
data.m$value <- as.numeric(data.m$value)
levels(data.m$Time.Code)

## Merge with site info
data.s <- merge(data.m, sites, by = "site")

## Sum across size classes within plots (samples)
data.p <- ddply(data.s, .(site, date1, Sample, Time.Code2, variable, dfw,order.dfw,area,salinity, shoot.density, fetch.jc), summarise, sum(value))

data.p$time.ID <- paste(data.p$site, data.p$Time.Code2, sep = '.') #could look at finer time resolution by using Time.Code here
names(data.p) <- c("site", "Date", "Sample", "Time.Code2", "species", "dfw","order","area","salinity","shoot.density","fetch","abundance", "TimeID")

## Merge with traits and sort by taxa or functional groups
data.tr <- merge(data.p, traits[,-1], by.x = "species", by.y = "species.names", all.x = TRUE, all.y = FALSE)


# Create datafiles for taxa and times -------------------------------------

## Remove all taxa that are not epifauna: 
levels(data.tr$eelgrss.epifauna)
data.e <- data.tr %>% filter(eelgrss.epifauna == c("yes", "sometimes"))
data.y <- data.tr %>% filter(eelgrss.epifauna == "yes")
data.s <- data.tr %>% filter(patch == "yes")
data.g <- data.tr %>% filter(function. == "grazer")
data.c <- data.tr %>% filter(group == "crustacean")
data.ga <- data.tr %>% filter(group == "gastropod")
data.tr <- data.e #reset to data.e

## group by sampling times
dataMAY <- data.tr[(data.tr$Time.Code2=="A"),]
dataJULY <- data.tr[(data.tr$Time.Code2=="C" & data.tr$site!="BE" & data.tr$site!="EI" & data.tr$site!="CC" & data.tr$site!="BI"),]
dataAUG <- data.tr[(data.tr$Time.Code2=="E"),]
data3times <- data.tr[(data.tr$site!="BE" & data.tr$site!="EI" & data.tr$site!="CC" & data.tr$site!="BI"),]
dataJULY9 <- data.tr[(data.tr$Time.Code2=="C"),]

## create site-level data by collapsing across plots
start.data <- dataJULY9 # dataMAY, data.mp, dataAUG, dataJULY, data3times
data.ms <- ddply(start.data, .(TimeID, area, species), summarise, sum(abundance)) #order
data2 <- dcast(data.ms, TimeID ~ species, mean) #order

# estimate number of species per site: 
dim(data2)
data.alpha <- data2
data.alpha <- specnumber(data.alpha[,1:31])
site.alpha <- as.data.frame(cbind(data2$TimeID, data.alpha))
names(site.alpha) <- c("site.time", "alpha")
site.alpha$site <- c("BE", "BI", "CB", "CC", "DC", "EI", "NB", "RP", "WI")

## for each plot, estimate relative abundance of grazers
## first remove filter feeders and predators:
data.t <- dataJULY9 #data3times # dataJULY9
data7 <- subset(data.t, function. != "unknown", select = c(1:17))
data7 <- subset(data.t, group != "echinoderm", select = c(1:17))
data8 <- ddply(data.t, .(site, Sample, order, dfw, area, salinity, fetch, function., Date, Time.Code2), summarise, sum(abundance)) #
data9 <- dcast(data8, site + Sample + fetch + area + order + dfw + Date + Time.Code2 ~ function., sum) 
data9$total <- (data9$detritovore + data9$suspension.feeder + data9$grazer + data9$predator + data9$filter.feeder) 
data9$pgrazer <- data9$grazer/data9$total
data9$pdet <- data9$detritovore/data9$total


# plots to visualize patterns; these plots are not presented in the manuscript figures 
par(mfrow=(c(2,2)))
plot((data9$pgrazer) ~ data9$fetch, pch = 19, xlab = 'Fetch', ylab = 'grazers / total', main = 'abundance / 0.28m2') #, ylim = c(0,1)
plot(log(data9$total) ~ data9$fetch, pch = 19, xlab = 'Fetch', ylab = 'ln(total inverts)', main = 'abundance / 0.28m2')
plot(log(data9$grazer+1) ~ data9$fetch, pch = 19, xlab = 'Fetch', ylab = 'ln(grazers+1)', main = 'abundance / 0.28m2')
plot(log(data9$filter.feeder + 1) ~ data9$fetch, pch = 19, xlab = 'Distance from Freshwater', ylab = 'ln(filter feeders)', main = 'abundance / 0.28m2')

plot(I(log(data9$grazer/(data9$filter.feeder+data9$suspension.feeder))) ~ data9$fetch, pch = 19, xlab = 'Fetch', ylab = 'ln(grazers/filter feeders)', main = 'abundance / 0.28m2')


# STATISTICAL ANALYSIS ON JULY ABUNDANCE PATTERNS 
## TABLE S5
hist(log(data9$total))
mod10 <- lm(log(data9$total+1) ~ data9$site)
mod20 <- lm(log(data9$grazer+1) ~ data9$site) 

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

model.sel(mod1a, mod1b, mod1, mod1c, mod1e, mod1f, mod10)
model.sel(mod2a, mod2b, mod2, mod2c, mod2e, mod2f, mod20)


## Abundance analyses for trends over time in sites sampled 3 times.
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



# Diversity analyses ------------------------------------------------------
## assemble diversity indices for Figure 2, Tables 1 and 3

div.data <- dcast(data.t[,c(1:4,12)], site + Date + Sample ~ species, sum)

H <- diversity(div.data[,-(c(1:3))], index ="shannon")
S <- diversity(div.data[,-(c(1:3))], index ="simpson")
I <- dispindmorisita(div.data[,-(c(1:3))], unique.rm = TRUE)
div.data$alpha.p <- specnumber(div.data[,4:33])
div.data$N <- rowSums(div.data[,(4:33)])

div.data %>%
group_by(site) %>% 
  select(5:10) %>% 
  as.matrix(.) %>% 
  map_df(., dispindmorisita, unique.rm = TRUE, na.rm = TRUE)

## calculate rarified richness for samples with more than 5 individuals
RR.data <- div.data %>%
  filter(N > 4) 
  
RR.data$RR = rarefy(RR.data[,-(c(1:3))], 5)


### Morisita's I within meadows for Table 1
I.BE <- dispindmorisita(div.data[(div.data$site=='BE'),-(c(1:3,34:35))], unique.rm = TRUE, na.rm = TRUE)
str(I.BE)
I.RP <- dispindmorisita(div.data[(div.data$site=='RP'),-(c(1:3,34:35))], unique.rm = TRUE, na.rm = TRUE)
I.DC <- dispindmorisita(div.data[(div.data$site=='DC'),-(c(1:3,34:35))], unique.rm = TRUE, na.rm = TRUE)
I.WI <- dispindmorisita(div.data[(div.data$site=='WI'),-(c(1:3,34:35))], unique.rm = TRUE, na.rm = TRUE)
I.CB <- dispindmorisita(div.data[(div.data$site=='CB'),-(c(1:3,34:35))], unique.rm = TRUE, na.rm = TRUE)
I.CC <- dispindmorisita(div.data[(div.data$site=='CC'),-(c(1:3,34:35))], unique.rm = TRUE, na.rm = TRUE)
I.EI <- dispindmorisita(div.data[(div.data$site=='EI'),-(c(1:3,34:35))], unique.rm = TRUE, na.rm = TRUE)
I.NB <- dispindmorisita(div.data[(div.data$site=='NB'),-(c(1:3,34:35))], unique.rm = TRUE, na.rm = TRUE)
I.BI <- dispindmorisita(div.data[(div.data$site=='BI'),-(c(1:3,34:35))], unique.rm = TRUE, na.rm = TRUE)

## assemble mean I values, and confidence intervals
means <- c(mean(I.DC[,4], na.rm = TRUE), mean(I.WI[,4], na.rm = TRUE), mean(I.BE[,4], na.rm = TRUE), mean(I.EI[,4], na.rm = TRUE), mean(I.RP[,4], na.rm = TRUE), mean(I.NB[,4], na.rm = TRUE), mean(I.CB[,4], na.rm = TRUE), mean(I.BI[,4], na.rm = TRUE), mean(I.CC[,4], na.rm = TRUE))

ci.upper <- function(x) mean(x) + 1.96*sd(x)/sqrt(length(x))
ci.lower <- function(x) mean(x) - 1.96*sd(x)/sqrt(length(x))

ci.lowers <- c(ci.lower(I.DC[,4]), ci.lower(I.WI[,4]), ci.lower(I.BE[,4]), ci.lower(I.EI[,4]), ci.lower(I.RP[,4]), ci.lower(I.NB[,4]), ci.lower(I.CB[,4]), ci.lower(I.BI[,4]), ci.lower(I.CC[,4]))

ci.uppers <- c(ci.upper(I.DC[,4]), ci.upper(I.WI[,4]), ci.upper(I.BE[,4]), ci.upper(I.EI[,4]), ci.upper(I.RP[,4]), ci.upper(I.NB[,4]), ci.upper(I.CB[,4]), ci.upper(I.BI[,4]), ci.upper(I.CC[,4]))

I.means <- cbind(means, ci.lowers, ci.uppers, c('DC', 'WI', 'BE', 'EI', 'RP', 'NB', 'CB', 'BI', 'CC'))
I.means <- as.data.frame(I.means)
names(I.means) <- c("I", "lower", "upper","site")

# Compile indices into one dataframe
div.summary <- cbind(div.data[c(1:3, 34:35)], H, S)
div.summaryE <- merge(div.summary, RR.data, by.x = c("site", "Date","alpha.p", "N", "Sample"), by.y = c("site", "Date","alpha.p","N", "Sample"))
div.summary2 <- merge(div.summaryE, sites, by.x = c("site"), by.y = c("site"))
div.summaryI <- merge(div.summary2, I.means, by.x = c("site"), by.y = c("site"))
div.summaryT <- merge(div.summaryI, site.alpha, by.x = c("site"), by.y = c("site")) # get site alpha from EMS code (this can be changed to reference code above in this file)
head(div.summaryT)
div.summaryT$apB <- as.numeric(as.character(div.summaryT$alpha)) - div.summaryT$alpha.p

## isolate each index
site.H <- ddply(div.summaryI, .(site), summarise, mean(H))
site.S <- ddply(div.summaryI, .(site), summarise, mean(S))
site.N <- ddply(div.data, .(site), summarise, mean(N))
site.Beta <- ddply(div.summaryT, .(site), summarise, mean(apB))
site.RR <- ddply(div.summaryT, .(site), summarise, mean(RR))

# ANALYZE UNIVARIATE DIVERSITY INDICES ----------------------------------------------

# Get and plot observed site means

# does plot level alpha diversity differ among meadows?
mods4 <- lm(div.summary2$alpha.p ~ div.summary2$site)
mods0 <-  lm(div.summary2$alpha.p ~ 1)
anova(mods4, mods0)
model.sel(mods4, mods0)
summary(mods4)

mods1 <- lm(div.summary2$H ~ div.summary2$site)
mods2 <- lm(div.summary2$S ~ div.summary2$site)
mods3 <- lm(div.summary2$N ~ div.summary2$site)
mods4 <- lm(div.summary2$alpha.p ~ div.summary2$site)

# i had date as a predictor, but it seems to be confounded with meadow ID, so I'm not going to include it as a predictor, within the June sampling time.
mod1a <- lm(div.summary2$H ~ div.summary2$fetch.jc)
mod2a <- lm(div.summary2$S ~ div.summary2$fetch.jc)
mod3a <- lm(div.summary2$alpha.p ~ div.summary2$fetch.jc)
mod4a <- lm(div.summary2$RR ~ div.summary2$fetch.jc)

mod1b <- lm(div.summary2$H ~ div.summary2$dfw)
mod2b <- lm(div.summary2$S ~ div.summary2$dfw)
mod3b <- lm(div.summary2$alpha.p ~ div.summary2$dfw)
mod4b <- lm(div.summary2$RR ~ div.summary2$dfw)

mod1c <- lm(div.summary2$H ~ div.summary2$dfw * div.summary2$fetch.jc)
mod2c <- lm(div.summary2$S ~ div.summary2$dfw * div.summary2$fetch.jc)
mod3c <- lm(div.summary2$alpha.p ~ div.summary2$dfw * div.summary2$fetch.jc)
mod4c <- lm(div.summary2$RR ~ div.summary2$dfw * div.summary2$fetch.jc)

mod1f <- lm(div.summary2$H ~ div.summary2$area * div.summary2$fetch.jc)
mod2f <- lm(div.summary2$S ~ div.summary2$area * div.summary2$fetch.jc)
mod3f <- lm(div.summary2$alpha.p ~ div.summary2$area * div.summary2$fetch.jc)
mod4f <- lm(div.summary2$RR ~ div.summary2$area * div.summary2$fetch.jc)

mod1d <- lm(div.summary2$H ~ 1)
mod2d <- lm(div.summary2$S ~ 1)
mod3d <- lm(div.summary2$alpha.p ~ 1)
mod4d <- lm(div.summary2$RR ~ 1)

mod1g <- lm(div.summary2$H ~ div.summary2$site)
mod2g <- lm(div.summary2$S ~ div.summary2$site)
mod3g <- lm(div.summary2$alpha.p ~ div.summary2$site)
mod4g <- lm(div.summary2$RR ~ div.summary2$site)

## TABLE S5
model.sel(mod1a, mod1b, mod1, mod1c, mod1d, mod1f, mod1g)
model.sel(mod2a, mod2b, mod2, mod2c, mod2d, mod2f, mod2g)
model.sel(mod3a, mod3b, mod3, mod3c, mod3d, mod3f, mod3g)
model.sel(mod4a, mod4b, mod4, mod4c, mod4d, mod4f, mod4g)


# FIGURE 2 ----------------------------------------------------------------

div.plot <- ggplot(data = div.summary2, aes(reorder(site, -order.dfw), y = alpha.p, ymin = 0, ymax = 12)) + 
  theme_bw() +
  geom_boxplot() +
  geom_point(x = 5, y = 10.5, pch = '*', size = 8, colour = "gray50") +
  geom_point(x = 7, y = 10.5, pch = '*', size = 8, colour = "gray50") +
  #scale_x_discrete(limits = order.dfw) +
  xlab("Site") +
  ylab("Species Richness")

div.plot
ggsave("Jan2017alphaplot.png", device = "png", width = 4, height = 2.5)

H.plot <- ggplot(data = div.summary2, aes(reorder(site, -order.dfw), y = H, ymin = 0, ymax = 2)) + 
  theme_bw() +
  geom_boxplot() +
  geom_point(x = 5, y = 1.9, pch = '*', size = 8, colour = "gray50") +
  xlab("Site") +
  ylab("Shannon Diversity")

H.plot
ggsave("Jan2017Hplot.png", device = "png", width = 4, height = 2.5)

S.plot <- ggplot(data = div.summary2, aes(reorder(site, -order.dfw), y = S, ymin = 0, ymax = 1)) + 
  theme_bw() +
  geom_boxplot() +
  geom_point(x = 5, y = 1, pch = '*', size = 8, colour = "gray50") +
  geom_point(x = 1, y = 1, pch = '*', size = 8, colour = "gray50") +
  xlab("Site") +
  ylab("Simpson Evenness")

S.plot
ggsave("Jan2017Splot.png", device = "png", width = 4, height = 2.5)


RR.plot <- ggplot(data = div.summary2, aes(reorder(site, -order.dfw), y = RR, ymin = 0, ymax = 4)) + 
  theme_bw() +
  geom_boxplot() +
  geom_point(x = 6, y = 3.5, pch = '*', size = 8, colour = "gray50") +
  geom_point(x = 4, y = 3.5, pch = '*', size = 8, colour = "gray50") +
  geom_point(x = 1, y = 3.5, pch = '*', size = 8, colour = "gray50") +
  xlab("Site") +
  ylab("RR") 

RR.plot
ggsave("Jan2017RRplot.png", device = "png", width = 4, height = 2.5)


## SOME OTHER THINGS NOT IN THE PAPER
## rank abundance curves
rankBE <- as.data.frame(rankabundance(div.data[(div.data$site=='BE'),-(c(1:3,34:35))]))
rankRP <- as.data.frame(rankabundance(div.data[(div.data$site=='RP'),-(c(1:3,34:35))]))
rankDC <- as.data.frame(rankabundance(div.data[(div.data$site=='DC'),-(c(1:3,34:35))]))
rankWI <- as.data.frame(rankabundance(div.data[(div.data$site=='WI'),-(c(1:3,34:35))]))
rankCB <- as.data.frame(rankabundance(div.data[(div.data$site=='CB'),-(c(1:3,34:35))]))
rankCC <- as.data.frame(rankabundance(div.data[(div.data$site=='CC'),-(c(1:3,34:35))]))
rankNB <- as.data.frame(rankabundance(div.data[(div.data$site=='NB'),-(c(1:3,34:35))]))
rankEI <- as.data.frame(rankabundance(div.data[(div.data$site=='EI'),-(c(1:3,34:35))]))
rankBI <- as.data.frame(rankabundance(div.data[(div.data$site=='BI'),-(c(1:3,34:35))]))

I.BE[order(I.BE$imor, decreasing = TRUE),]
I.RP[order(I.RP$imor, decreasing = TRUE),]
I.DC[order(I.DC$imor, decreasing = TRUE),]
I.WI[order(I.WI$imor, decreasing = TRUE),]
I.CB[order(I.CB$imor, decreasing = TRUE),]
I.CC[order(I.CC$imor, decreasing = TRUE),]
I.NB[order(I.NB$imor, decreasing = TRUE),]
I.EI[order(I.EI$imor, decreasing = TRUE),]
I.BI[order(I.BI$imor, decreasing = TRUE),]

## now rank by Pchisq
I.BE[order(I.BE$pchisq),]
I.RP[order(I.RP$pchisq),]
I.DC[order(I.DC$pchisq),]
I.WI[order(I.WI$pchisq),]
I.CB[order(I.CB$pchisq),]
I.CC[order(I.CC$pchisq),]
I.NB[order(I.NB$pchisq),]
I.EI[order(I.EI$pchisq),]
I.BI[order(I.BI$pchisq),]

## for the whole dataset: 
div.data <- dcast(data.tr[,c(1:4,12)], site + Date + Sample ~ species, sum)
rank.all <- as.data.frame(rankabundance(div.data[,-(c(1:3,34:35))]))


