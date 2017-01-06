########################
### Whippo et al eelgrass epifaunal data for Barkley Sound
### elements of metacommunity structure (EMS)
### Code written by Mary O'Connor (with help from Domonik Bahlburg)
### started March 2016
### modified file in Jan 2017 to streamline wiht 'Marys diversity analyses.R' file
#########################



# Load libraries ----------------------------------------------------------

library(metacom)
library(broom)

#not these if working with Marys diversity analyses.R
library(plyr)
library(dplyr)
library(Matrix)
library(lattice)
library(reshape2)


# Data preparation --------------------------------------------------------

## start with data2 file:

## clean out any empty species or empty sites
rowSums(data2[-1]) -> data2$totals
data2 <- data2[which(data2$totals != "0"),] 
#data3 <- data2[-which(is.na(data2$totals)),] 
rownames(data2) <- data2[,1]
data3 <- data2[, -c(1, ncol(data2))]
data3[data3 > 0] <- 1
cols.to.delete <- which(colSums(data3) == '0')
data4 <- data3[, -(cols.to.delete)]

#remove NaN and Nas. 
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
data4[is.nan.data.frame(data4)] <- 0
data4[is.na(data4)] <- 0
data5 <- data4[,-ncol(data4)]


### run metacommunity analysis 
Metacommunity(data4, verbose = TRUE, order = FALSE) -> meta
meta[2:4]

a <- as.data.frame(meta[1])

pdf('July 9 sites epifauna.pdf', width = 7, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="July 9 epifauna YES",  border="black", scales = list(cex = c(0.4, 0.4), x = list(rot = c(90))))
dev.off()

