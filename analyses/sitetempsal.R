### temp sal situation

sites <- read.csv("./data/site.info.csv")
abiotic <- read.csv("./data/TempSalData.csv")

library(dplyr)
library(tidyr)
library(nlme)

head(abiotic)

temps <- abiotic %>%
  select(Site, Time.Code, Temp.a, Sal.a) %>%
  group_by(Site, Time.Code) %>%
  summarise(avgT = mean(Temp.a)) 

  
View(temps)

Sal <- abiotic %>%
  select(Site, Time.Code, Temp.a, Sal.a) %>%
  group_by(Site, Time.Code) %>%
  summarise(avgS = mean(Sal.a))

sites2 <- merge.data.frame(temps %>%
                             filter(Time.Code == "B"), 
                           sites, by = "Site")

plot(temps$avgT ~ temps$Site)

View(sites)

plot(abiotic$Sal.a ~ abiotic$Temp.a)


## ok i think we need to do a model comparison with site-level predictors (so merge site dataframe), and see predictors are needed in the biotic data (lai, etc?)

mod1 <- lme(Sal.a ~ Time.Code + Temp.a, random = ~ 1 | Site, data = abiotic, na.action = na.omit)
summary(mod1)
