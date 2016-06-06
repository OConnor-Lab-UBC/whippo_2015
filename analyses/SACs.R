## whippo et al
### species accumulation curves

data <- read.csv("./analyses/rawcomm.csv")

## melt and recast so the datafile goes from long to wide (species as columns)
data.m <- melt(data, id = c(1,2,3,4,5,52))

levels(data.m$Time.Code)
data.m$Time.Code <- as.character(data.m$Time.Code)
data.m$Time.Code[data.m$Time.Code == "C "] <- "C"
data.m$Time.Code <- as.factor(data.m$Time.Code)
#data.m$Sieve <- as.factor(data.m$Sieve)
data.m$value <- as.numeric(data.m$value)
levels(data.m$Time.Code)

head(data.m)

## sum across size classes within samples
data.mp <- ddply(data.m, .(site, Time.Code, Sample, Time.Code2, variable), summarise, sum(value))
data.mp$time.ID <- paste(data.mp$site, data.mp$Time.Code2, sep = '.') #could look at finer time resolution by using Time.Code here
names(data.mp) <- c("site", "Time.Code", "Sample", "Time.Code2", "species", "abundance", "TimeID")

## merge with traits and sort by taxa or functional groups
data.tr <- merge(data.mp, traits[,-1], by.x = "species", by.y = "species.names", all.x = TRUE, all.y = FALSE)

## remove all taxa that are not epifauna, because they were not evenly sampled across samples and meadows: 
levels(data.tr$eelgrss.epifauna)
data.e <- data.tr %>% filter(eelgrss.epifauna == c("yes", "sometimes"))
data.y <- data.tr %>% filter(eelgrss.epifauna == "yes")
data.tr <- data.e

## diverge from EMS code, and dcast to form species by plot / site / time matrix.
comm <- data.tr
data.ms <- ddply(comm, .(TimeID, Sample, species), summarise, sum(abundance))
data2 <- dcast(data.ms[(data.ms$TimeID == 'BE.C'),], Sample ~ species, mean)

#remove NaN and Nas. 
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
data2[is.nan.data.frame(data2)] <- 0
data2[is.na(data2)] <- 0

sp1 <- specaccum(data2[,-1], "random")




par
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", xlab = 'plot', ylab = 'species', ylim = c(0, 30), main = 'DC.C')

plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", xlab = 'plot', ylab = 'species', ylim = c(0, 30), main = 'RP.C')

plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", xlab = 'plot', ylab = 'species', ylim = c(0, 30), main = 'WI.C')

plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", xlab = 'plot', ylab = 'species', ylim = c(0, 30), main = 'CC.C')

plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", xlab = 'plot', ylab = 'species', ylim = c(0, 30), main = 'NB.C')

plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", xlab = 'plot', ylab = 'species', ylim = c(0, 30), main = 'EI.C')

plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", xlab = 'plot', ylab = 'species', ylim = c(0, 30), main = 'BI.C')

plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", xlab = 'plot', ylab = 'species', ylim = c(0, 30), main = 'CB.C')

plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", xlab = 'plot', ylab = 'species', ylim = c(0, 30), main = 'BE.C')
