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
library(Rcmdr)
library(BiodiversityR)
library(plyr)
library(reshape2)
library(dplyr)
library(MuMIn)
library(lubridate)
library(ggplot2)
library(purrr)

data.tr <- read.csv("./data/Whippodata.csv")
data.tr <- data.tr[,-1]
sites <- read.csv("./data/site.info.csv")

# Create datafiles for taxa and times -------------------------------------

## Remove all taxa that are not epifauna: 
levels(data.tr$eelgrss.epifauna)
data.e <- data.tr %>% filter(eelgrss.epifauna == "yes" | eelgrss.epifauna == "sometimes")
data.n <- data.tr %>% filter(eelgrss.epifauna == "no")
data.tr <- data.e #reset to data.e

## create data subsets for different sampling times
dataMAY <- data.tr[(data.tr$Time.Code2=="A"),]
dataJULY <- data.tr[(data.tr$Time.Code2=="C" & data.tr$site!="BE" & data.tr$site!="EI" & data.tr$site!="CC" & data.tr$site!="BI"),]
dataAUG <- data.tr[(data.tr$Time.Code2=="E"),]
data3times <- data.tr[(data.tr$site!="BE" & data.tr$site!="EI" & data.tr$site!="CC" & data.tr$site!="BI"),]
dataJULY9 <- data.tr[(data.tr$Time.Code2=="C"),]

## create site-level data by collapsing across plots
start.data <- dataJULY9 # dataMAY, dataAUG, dataJULY, data3times
data.ms <- ddply(start.data, .(TimeID, area, species), summarise, sum(abundance)) #order
data2 <- dcast(data.ms, TimeID ~ species, mean) #order

## Table 1: gamma results
# observed number of species per site for month identified in start.data: 
dim(data2)
data.alpha <- data2
data.alpha <- specnumber(data.alpha[,2:35])
site.alpha <- as.data.frame(cbind(as.character(data2$TimeID), data.alpha))
names(site.alpha) <- c("site.time", "alpha")
site.alpha$site <- c("BE", "BI", "CB", "CC", "DC", "EI", "NB", "RP", "WI")

#info for Table 1: 
site.alpha

## for each plot, estimate relative abundance of grazers for period in start.data
data.t <- start.data 

# STATISTICAL ANALYSIS ON JULY ABUNDANCE PATTERNS 
# Diversity analyses ------------------------------------------------------
## assemble diversity indices for Tables 1 and 3, Figure 2

div.data <- dcast(data.t[,c(1:4,12)], site + Date + Sample ~ species, sum)

H <- diversity(div.data[,-(c(1:3))], index ="shannon") #sample-level H
S <- diversity(div.data[,-(c(1:3))], index ="simpson")
I <- dispindmorisita(div.data[,-(c(1:3))], unique.rm = TRUE)
div.data$alpha.p <- specnumber(div.data[,4:37])
div.data$N <- rowSums(div.data[,(4:37)])


### Morisita's I within meadows for Table 1
I.BE <- dispindmorisita(div.data[(div.data$site=='BE'),-(c(1:3,38:39))], unique.rm = TRUE, na.rm = TRUE)
I.RP <- dispindmorisita(div.data[(div.data$site=='RP'),-(c(1:3,38:39))], unique.rm = TRUE, na.rm = TRUE)
I.DC <- dispindmorisita(div.data[(div.data$site=='DC'),-(c(1:3,38:39))], unique.rm = TRUE, na.rm = TRUE)
I.WI <- dispindmorisita(div.data[(div.data$site=='WI'),-(c(1:3,38:39))], unique.rm = TRUE, na.rm = TRUE)
I.CB <- dispindmorisita(div.data[(div.data$site=='CB'),-(c(1:3,38:39))], unique.rm = TRUE, na.rm = TRUE)
I.CC <- dispindmorisita(div.data[(div.data$site=='CC'),-(c(1:3,34:35))], unique.rm = TRUE, na.rm = TRUE)
I.EI <- dispindmorisita(div.data[(div.data$site=='EI'),-(c(1:3,38:39))], unique.rm = TRUE, na.rm = TRUE)
I.NB <- dispindmorisita(div.data[(div.data$site=='NB'),-(c(1:3,38:39))], unique.rm = TRUE, na.rm = TRUE)
I.BI <- dispindmorisita(div.data[(div.data$site=='BI'),-(c(1:3,38:39))], unique.rm = TRUE, na.rm = TRUE)

## proportion of species with significant I
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

## assemble mean I values, and confidence intervals for TABLE 1
means <- c(mean(I.DC[,4], na.rm = TRUE), mean(I.WI[,4], na.rm = TRUE), mean(I.BE[,4], na.rm = TRUE), mean(I.EI[,4], na.rm = TRUE), mean(I.RP[,4], na.rm = TRUE), mean(I.NB[,4], na.rm = TRUE), mean(I.CB[,4], na.rm = TRUE), mean(I.BI[,4], na.rm = TRUE), mean(I.CC[,4], na.rm = TRUE))

ci.upper <- function(x) mean(x) + 1.96*sd(x)/sqrt(length(x))
ci.lower <- function(x) mean(x) - 1.96*sd(x)/sqrt(length(x))

ci.lowers <- c(ci.lower(I.DC[(I.DC$imst != "NaN"),4]), ci.lower(I.WI[(I.WI$imst != "NaN"),4]), ci.lower(I.BE[(I.BE$imst != "NaN"),4]), ci.lower(I.EI[(I.EI$imst != "NaN"),4]), ci.lower(I.RP[(I.RP$imst != "NaN"),4]), ci.lower(I.NB[(I.NB$imst != "NaN"),4]), ci.lower(I.CB[(I.CB$imst != "NaN"),4]), ci.lower(I.BI[(I.BI$imst != "NaN"),4]), ci.lower(I.CC[(I.CC$imst != "NaN"),4]))

ci.uppers <- c(ci.upper(I.DC[(I.DC$imst != "NaN"),4]), ci.upper(I.WI[(I.WI$imst != "NaN"),4]), ci.upper(I.BE[(I.BE$imst != "NaN"),4]), ci.upper(I.EI[(I.EI$imst != "NaN"),4]), ci.upper(I.RP[(I.RP$imst != "NaN"),4]), ci.upper(I.NB[(I.NB$imst != "NaN"),4]), ci.upper(I.CB[(I.CB$imst != "NaN"),4]), ci.upper(I.BI[(I.BI$imst != "NaN"),4]), ci.upper(I.CC[(I.CC$imst != "NaN"),4]))

Nsig <- c((length(I.DC[(I.DC$pchisq != "NaN" & I.DC$pchisq < 0.01),1])/length(I.DC[(I.DC$imor != "NaN"),1])), (length(I.WI[(I.WI$pchisq != "NaN" & I.WI$pchisq < 0.01),1])/length(I.WI[(I.WI$imor != "NaN"),1])), (length(I.BE[(I.BE$pchisq != "NaN" & I.BE$pchisq < 0.01),1])/length(I.BE[(I.BE$imor != "NaN"),1])), (length(I.EI[(I.EI$pchisq != "NaN" & I.EI$pchisq < 0.01),1])/length(I.EI[(I.EI$imor != "NaN"),1])), (length(I.RP[(I.RP$pchisq != "NaN" & I.RP$pchisq < 0.01),1])/length(I.RP[(I.RP$imor != "NaN"),1])), (length(I.NB[(I.NB$pchisq != "NaN" & I.NB$pchisq < 0.01),1])/length(I.NB[(I.NB$imor != "NaN"),1])), (length(I.CB[(I.CB$pchisq != "NaN" & I.CB$pchisq < 0.01),1])/length(I.CB[(I.CB$imor != "NaN"),1])), (length(I.BI[(I.BI$pchisq != "NaN" & I.BI$pchisq < 0.01),1])/length(I.BI[(I.BI$imor != "NaN"),1])), (length(I.CC[(I.CC$pchisq != "NaN" & I.CC$pchisq < 0.01),1])/length(I.CC[(I.CC$imor != "NaN"),1])))

I.means <- cbind(means, ci.lowers, ci.uppers, Nsig, c('DC', 'WI', 'BE', 'EI', 'RP', 'NB', 'CB', 'BI', 'CC'))
I.means <- as.data.frame(I.means)
names(I.means) <- c("I", "lower", "upper","Nsig","site")

#information for Table 1:
I.means

### CODE FOR TABLE 2 HERE
## Order species by abundance. There is a command for this in BiodiversityR, but it is no longer working R. So doing it by hand.

BEsp <- div.data %>%
  filter(site == "BE") %>%
  select(., -(c(1:3,38:39))) %>%
  summarise_all(., funs(sum)) %>%
  tidyr::gather(., "species", "N") %>%
  arrange(., desc(N)) %>%
  filter(N > 0)

rankBE <- BEsp %>%
  mutate(., rank = min_rank(desc(N))) %>%
  View(.)

DCsp <- div.data %>%
  filter(site == "DC") %>%
  select(., -(c(1:3,38:39))) %>%
  summarise_all(., funs(sum)) %>%
  tidyr::gather(., "species", "N") %>%
  arrange(., desc(N)) %>%
  filter(N > 0)

rankDC <- DCsp %>%
  mutate(., rank = min_rank(desc(N))) %>%
  View(.)

WIsp <- div.data %>%
  filter(site == "WI") %>%
  select(., -(c(1:3,38:39))) %>%
  summarise_all(., funs(sum)) %>%
  tidyr::gather(., "species", "N") %>%
  arrange(., desc(N)) %>%
  filter(N > 0)

rankWI <- WIsp %>%
  mutate(., rank = min_rank(desc(N))) %>%
  View(.)

EIsp <- div.data %>%
  filter(site == "EI") %>%
  select(., -(c(1:3,38:39))) %>%
  summarise_all(., funs(sum)) %>%
  tidyr::gather(., "species", "N") %>%
  arrange(., desc(N)) %>%
  filter(N > 0)

rankEI <- EIsp %>%
  mutate(., rank = min_rank(desc(N))) %>%
  View(.)

RPsp <- div.data %>%
  filter(site == "RP") %>%
  select(., -(c(1:3,38:39))) %>%
  summarise_all(., funs(sum)) %>%
  tidyr::gather(., "species", "N") %>%
  arrange(., desc(N)) %>%
  filter(N > 0)

rankRP <- RPsp %>%
  mutate(., rank = min_rank(desc(N))) %>%
  View(.)

NBsp <- div.data %>%
  filter(site == "NB") %>%
  select(., -(c(1:3,38:39))) %>%
  summarise_all(., funs(sum)) %>%
  tidyr::gather(., "species", "N") %>%
  arrange(., desc(N)) %>%
  filter(N > 0)

rankNB <- NBsp %>%
  mutate(., rank = min_rank(desc(N))) %>%
  View(.)

CBsp <- div.data %>%
  filter(site == "CB") %>%
  select(., -(c(1:3,38:39))) %>%
  summarise_all(., funs(sum)) %>%
  tidyr::gather(., "species", "N") %>%
  arrange(., desc(N)) %>%
  filter(N > 0)

rankCB <- CBsp %>%
  mutate(., rank = min_rank(desc(N))) %>%
  View(.)

BIsp <- div.data %>%
  filter(site == "BI") %>%
  select(., -(c(1:3,38:39))) %>%
  summarise_all(., funs(sum)) %>%
  tidyr::gather(., "species", "N") %>%
  arrange(., desc(N)) %>%
  filter(N > 0)

rankBI <- BIsp %>%
  mutate(., rank = min_rank(desc(N))) %>%
  View(.)

CCsp <- div.data %>%
  filter(site == "CC") %>%
  select(., -(c(1:3,38:39))) %>%
  summarise_all(., funs(sum)) %>%
  tidyr::gather(., "species", "N") %>%
  arrange(., desc(N)) %>%
  filter(N > 0)

rankCC <- CCsp %>%
  mutate(., rank = min_rank(desc(N))) %>%
  View(.)

ALLsp <- div.data %>%
  #filter(site == "CC") %>%
  select(., -(c(1:3,38:39))) %>%
  summarise_all(., funs(sum)) %>%
  tidyr::gather(., "species", "N") %>%
  arrange(., desc(N)) %>%
  filter(N > 0)

rankALL <- ALLsp %>%
  mutate(., rank = min_rank(desc(N))) %>%
  View(.)

### code for rarified richness here below, after ranks are generated:
##
library(iNEXT)
library(ggplot2)

#create data input: 
RBE <- rankBE[,2]
RRP <- rankRP[,2]
RDC <- rankDC[,2]
RWI <- rankWI[,2]
RCC <- rankCC[,2]
RCB <- rankCB[,2]
RNB <- rankNB[,2]
REI <- rankEI[,2]
RBI <- rankBI[,2]






# Compile indices into one dataframe to make FIGURE 2
div.summary <- cbind(div.data[c(1:3, 38:39)], H, S)
#div.summaryE <- merge(div.summary, RR.data, by.x = c("site", "Date","Sample", "alpha.p", "N"), by.y = c("site", "Date","alpha.p","N", "Sample"))
div.summary2 <- merge(div.summary, sites, by.x = c("site"), by.y = c("site"))
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

## check distributions: 
hist(div.summary2$alpha.p)

# does plot level alpha diversity differ among meadows?
mods4 <- lm(div.summary2$alpha.p ~ div.summary2$site)
mods0 <-  lm(div.summary2$alpha.p ~ 1)
anova(mods4, mods0)
model.sel(mods4, mods0)
summary(mods4)

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

## TABLE 3
model.sel(mod1a, mod1b, mod1c, mod1d, mod1f, mod1g)
model.sel(mod2a, mod2b, mod2c, mod2d, mod2f, mod2g)
model.sel(mod3a, mod3b, mod3c, mod3d, mod3f, mod3g)
model.sel(mod4a, mod4b, mod4c, mod4d, mod4f, mod4g)



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
ggsave("Julyalphaplot.png", device = "png", width = 4, height = 2.5)

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


ggsave("histogram", device = "png", width = 4, height = 4)


# Table S5 ----------------------------------------------------------------
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




