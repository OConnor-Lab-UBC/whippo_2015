
library(vegan)
#library(plyr)
library(dplyr)
library(ggplot2)

#read in data
data <- read.csv("rawcomm.csv", header = TRUE)
traits <- read.csv("grazertraits.csv", header = TRUE)
sites <- read.csv("site.info.csv", header = TRUE)

### I want to look at patterns: 
### within grazers
### within > 2 mm

### to do this, I think the surest way is to 
### a) collapse the rawcomm to have no size classes
### b) merge the rawcomm with the grazer traits.
### c) then create a subsetted dataset for only species that I want. 

## a) 
## create an id for each sample at each site, time. 
idmaker = function(x) { return(paste(x, collapse=""))}
data$ID <- apply(data[,c("site", "Time.Code2", "Sample")], 1, idmaker)

## now treat each sample as a group collapse across sizes to sum taxa
data.s <- group_by(data, ID)
func <- function(x) { sum(x) }
data.sum <- as.data.frame(summarise_each(data.s, funs(sum), (Idotea.resecata:Alia.carinata)))


### b) merge the rawcomm with the grazer traits. [subset the data by functional group in this step!]
library(reshape2)
data3 <- melt(data.sum, id.vars = "ID", variable.name = 'species')
data3.1 <- merge(data3, traits, by.x = 'variable', by.y = 'species.names', all = TRUE)

## subset out grazers and estimate diversity
# data.gr <- data3.1[which(data3.1$function. == "grazer"),]
data.gr <- filter(data3.1, function. == "grazer")
data.grz <- data.frame(cast(data.gr, ID ~ variable, list))

data.grz$H <- diversity(data.grz[,2:17], index = "shannon") #not working... not sure why.
data.grz$S <- diversity(data.grz[,2:17], index = "simpson")
data.grz$N <- apply(data.grz[,2:17], 1, function(x) sum(x))
head(data.grz)

## add in sample information 
data.desc <- ddply(data.s, .(site, Time.Code2, Sample), summarise, unique(ID))
data.grz2 <- merge(data.grz, data.desc, by.x = "ID", by.y = "unique(ID)", all.x = TRUE, all.y = FALSE)

## add in site information
data.grz3 <- merge(data.grz2, sites, by.x = 'site', by.y = 'site')

## look at it
plot(log(data.grz3$N) ~ data.grz3$order*data.grz3$Time.Code2)

grzN <- ggplot(data.grz3, aes(factor(order), log(N)))
grzN_abund <- grzN + geom_boxplot(aes(fill = factor(Time.Code2))) 

grzS <- ggplot(data.grz3, aes(factor(order), S))
grzS_plot <- grzS + geom_boxplot(aes(fill = factor(Time.Code2))) 





######################
## extra crap
######################

## estimating diversity indices
#H <- diversity(data[,6:51], index = "shannon")
#S <- diversity(data[,6:51], index = "simpson")
#N <- apply(data[,6:51], 1, function(x) sum(x))

#data$H <- H
#data$S <- S
#data$N <- N
#head(data)

#hist(data[(data$Sieve == 1),]$N)
#hist(data[(data$Sieve == 2),]$N)
#hist(data[(data$Sieve == 4),]$N)
#hist(data[(data$Sieve == 8),]$N)

#plot(data$N, data$Sieve)


data.sum3 <- merge(data.sum2, sites, by.x = "site", by.y = "site", all.x = TRUE, all.y = FALSE)
head(data.sum3)

data5 <- cast(data3, ID ~ Species, list)



## simpler way to transpose: 
tdata.sum2 <- t(data.sum2)
head(tdata.sum2)
tdata.sum2df <- data.frame(tdata.sum2)
names(tdata.sum2df) <- as.character(tdata.sum2df[1,])
data3.1 <- merge(tdata.sum2df, traits, by.x = 'ID', by.y = 'species.names', all = TRUE)
## but struggling to get column names to stick; can't figure this out.
tdata.sum2dfa <- tdata.sum2df[-1,]
colnames(tdata.sum2dfa) <- c(tdata.sum2df$ID)
tdata.sum2dfa[1:10, 1:10]






## exploring dplyr



data2 <- ddply(data, .(site, date, sample, Time.Code, Time.Code2), summarise, cbind(sum(data$Idotea.resecata), sum(data$Phllaplysia.taylori)))

data[(data$site == "DC" & data$Time.Code == "A"),]

filter(data, site == 'DC', Time.Code2 == 'A')

select(data, 6:51)

mutate(data, 
       N = sum(select(data, Idotea.resecata:Alia.carinata)))

summarise_each_(data, ID = data$site + data$Time.Code2 + data$sample)


data.g <- group_by()
