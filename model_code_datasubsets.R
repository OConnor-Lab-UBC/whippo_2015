
library(vegan)
library(plyr)

#read in data
data <- read.csv("rawcomm.csv")
traits <- read.csv("grazertraits.csv")

H <- diversity(data[,6:51], index = "shannon")
S <- diversity(data[,6:51], index = "simpson")
N <- ddply(data, .(site, date, Sample), summarise, rowSums(data[,6:51]))
N <- apply(data[,6:51], 1, function(x) sum(x))

data$H <- H
data$S <- S
data$N <- N
head(data)

hist(data[(data$Sieve == 1),]$N)
hist(data[(data$Sieve == 2),]$N)
hist(data[(data$Sieve == 4),]$N)
hist(data[(data$Sieve == 8),]$N)

plot(data$N, data$Sieve)

### I want to look at patterns: 
### within grazers
### within > 2 mm

### to do this, I think the surest way is to 
### a) collapse the rawcomm to have no size classes
### b) merge the rawcomm with the grazer traits.
### c) then create a subsetted dataset for only species that I want. 
