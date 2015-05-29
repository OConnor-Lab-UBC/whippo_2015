library(nlme)


#Read in data

data=read.csv("./plot_data.csv")

data.cs <- data[(data$Site!='BI')&(data$Site!='CC')&(data$Site!='EI')&(data$Site!='BE'),]

head(data)

#basic anova:

mod1 <- lm(Rarefied_Richness ~ Site*Time, data = data.cs, na.action = na.omit)

summary(mod1)


# or treat site as a random effect

mod2 <- lme(Rarefied_Richness ~ Time, random = ~ 1 + Time | Site, data = data.cs, na.action = na.omit)

mod2a <- lme(Rarefied_Richness ~ Time, random = ~ 1 | Site, data = data.cs, na.action = na.omit)

mod2b <- lm(Rarefied_Richness ~ Time, data = data.cs, na.action = na.omit)

anova(mod2, mod2a)