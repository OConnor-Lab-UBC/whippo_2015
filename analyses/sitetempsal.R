### temp sal situation

sites <- read.csv("./data/site.info.csv")
abiotic <- read.csv("./data/TempSalData.csv")

library(dplyr)
library(tidyr)
library(nlme)
library(MuMIn)

head(abiotic)

temps <- abiotic %>%
  select(Site, Time.Code, Temp.a, Sal.a) %>%
  group_by(Site, Time.Code) #%>%
  summarise(avgT = mean(Temp.a)) 
 
View(temps)

plot(temps$Temp.a ~ temps$Sal.a)

sites2 <- merge.data.frame(temps %>%
                             filter(Time.Code == "B"), 
                           sites, by.x = "Site", by.y = "site") 

names(sites2)

plot(sites2$Temp.a ~ sites2$Sal.a)


## ok i think we need to do a model comparison with site-level predictors and see predictors are needed in the biotic data (lai, etc?)

## abiotic properties first: is salinity predicted by temp, dfw or fetch?
plot(sites2$Sal.a ~ sites2$Temp.a)
plot(sites2$Sal.a ~ sites2$dfw)
plot(sites2$Sal.a ~ sites2$fetch.jc)
plot(sites2$Temp.a ~ sites2$dfw)

### Salinity
mod1 <- lme(Sal.a ~ I(log(Temp.a)-mean(log(Temp.a))), random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod2 <- lme(Sal.a ~ I(log(Temp.a)-mean(log(Temp.a))) + dfw, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod3 <- lme(Sal.a ~ I(log(Temp.a)-mean(log(Temp.a))) + dfw + fetch.jc, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod4 <- lme(Sal.a ~ 1, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod5 <- lme(Sal.a ~ dfw, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod6 <- lme(Sal.a ~ I(log(Temp.a)-mean(log(Temp.a)))*dfw, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod7 <- lme(Sal.a ~ fetch.jc, random = ~ 1 | Site, data = sites2, na.action = na.omit)

model.sel(mod1, mod2, mod3, mod4, mod5, mod6, mod7)

model.avg(mod6, mod2) -> mod.sal
confint(mod.sal)
summary(mod.sal)


### Temperature
mod1 <- lme(log(Temp.a) ~ I(Sal.a - mean(Sal.a)), random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod2 <- lme(log(Temp.a) ~ I(Sal.a - mean(Sal.a)) + dfw, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod3 <- lme(log(Temp.a) ~ I(Sal.a - mean(Sal.a)) + dfw + fetch.jc, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod4 <- lme(log(Temp.a) ~ 1, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod5 <- lme(log(Temp.a) ~ dfw, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod6 <- lme(log(Temp.a) ~ I(Sal.a - mean(Sal.a))*dfw, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod7 <- lme(log(Temp.a) ~ fetch.jc, random = ~ 1 | Site, data = sites2, na.action = na.omit)

model.sel(mod1, mod2, mod3, mod4, mod5, mod6, mod7)

summary(mod4)
confint(mod4)

## how correlated are temp and salinity?
mod0 <- lme(log(Temp.a) ~ I(Sal.a - mean(Sal.a)), random = ~ Sal.a | Site, data = sites2, method = "REML")
mod1 <- lme(log(Temp.a) ~ I(Sal.a - mean(Sal.a)), random = ~ 1 | Site, data = sites2, method = "REML")
model.sel(mod1, mod0)

mod1 <- lme(log(Temp.a) ~ I(Sal.a - mean(Sal.a)), random = ~ 1 | Site, data = sites2, method = "ML")
mod2 <- lme(log(Temp.a) ~ 1, random = ~ 1 | Site, data = sites2, method = "ML")
model.sel(mod1, mod2)

summary(mod1) #summary tells us these are correlated, but the relationship is quite weak: -0.007 C / 1 PPT change in salinity.

## to use temp and salinity as predictors of 

temps <- abiotic %>%
  select(Site, Time.Code, Temp.a, Sal.a) %>%
  group_by(Site, Time.Code) %>%
summarise(avgT = mean(Temp.a))

sals <- abiotic %>%
  select(Site, Time.Code, Temp.a, Sal.a) %>%
  group_by(Site, Time.Code) %>%
  summarise(avgS = mean(Sal.a))

sites3 <- merge.data.frame(sites,
                           temps %>% 
                           filter(Time.Code == "B"), by.x = "site", by.y = "Site") 
