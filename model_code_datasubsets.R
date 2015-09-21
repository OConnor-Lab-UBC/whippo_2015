
library(vegan)
library(plyr)
library(dplyr)

#read in data
data <- read.csv("rawcomm.csv")
traits <- read.csv("grazertraits.csv")

H <- diversity(data[,6:51], index = "shannon")
S <- diversity(data[,6:51], index = "simpson")
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

## a) 
## create an id for each sample at each site, time. 
idmaker = function(x) { return(paste(x, collapse=""))}
data$ID <- apply(data[,c("site", "Time.Code2", "Sample")], 1, idmaker)

## now treat each sample as a group collapse across sizes to sum taxa
data.s <- group_by(data, ID)
func <- function(x) { sum(x) }
data.sum <- summarise_each(data.s, funs(sum), (Idotea.resecata:Alia.carinata))
data.desc <- ddply(data.s, .(site, Time.Code2, Sample), summarise, unique(ID))
data.sum2 <- merge(data.sum, data.desc, by.x = "ID", by.y = "unique(ID)", all.x = TRUE, all.y = FALSE)
head(data.sum2)

### b) merge the rawcomm with the grazer traits.
library(reshape)
data3 <- melt(data.sum2, id = "ID", variable = 'Species')
data3.1 <- merge(data3, traits, by.x = 'Species', by.y = 'species.names', all = TRUE)
## start here: remove extra columns, then cast data back to original shape, and calculate diversity indices.
data4[,-(data4$species & data4$native & data4$larvae & data4$tube.building & data4$other.notes)]

data5 <- cast(data4, ID ~ Species, list)
## exploring dplyr



data2 <- ddply(data, .(site, date, sample, Time.Code, Time.Code2), summarise, cbind(sum(data$Idotea.resecata), sum(data$Phllaplysia.taylori)))

data[(data$site == "DC" & data$Time.Code == "A"),]

filter(data, site == 'DC', Time.Code2 == 'A')

select(data, 6:51)

mutate(data, 
       N = sum(select(data, Idotea.resecata:Alia.carinata)))

summarise_each_(data, ID = data$site + data$Time.Code2 + data$sample)


data.g <- group_by()
