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
mod1 <- lme(Sal.a ~ Temp.a, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod2 <- lme(Sal.a ~ Temp.a + dfw, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod3 <- lme(Sal.a ~ Temp.a + dfw + fetch.jc, random = ~ 1 | sites2, data = sites2, na.action = na.omit)
mod4 <- lme(Sal.a ~ 1, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod5 <- lme(Sal.a ~ dfw, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod6 <- lme(Sal.a ~ Temp.a*dfw, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod7 <- lme(Sal.a ~ fetch.jc, random = ~ 1 | Site, data = sites2, na.action = na.omit)

model.sel(mod1, mod2, mod3, mod4, mod5, mod6, mod7)

model.avg(mod5, mod2) -> mod.sal
confint(mod.sal)

summary(mod1)


### Temperature
mod1 <- lme(Temp.a ~ Sal.a, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod2 <- lme(Temp.a ~ Sal.a + dfw, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod3 <- lme(Temp.a ~ Sal.a + dfw + fetch.jc, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod4 <- lme(Temp.a ~ 1, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod5 <- lme(Temp.a ~ dfw, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod6 <- lme(Temp.a ~ Sal.a*dfw, random = ~ 1 | Site, data = sites2, na.action = na.omit)
mod7 <- lme(Temp.a ~ fetch.jc, random = ~ 1 | Site, data = sites2, na.action = na.omit)

model.sel(mod1, mod2, mod3, mod4, mod5, mod6, mod7)

model.avg(mod5, mod4, mod1) -> mod.t
summary(mod.t)
confint(mod.t)



## other code I'm no longer using
Sal <- abiotic %>%
  select(Site, Time.Code, Temp.a, Sal.a) %>%
  group_by(Site, Time.Code) %>%
  summarise(avgS = mean(Sal.a))

plot(temps$avgT ~ temps$Site)

View(sites)

plot(abiotic$Sal.a ~ abiotic$Temp.a)



