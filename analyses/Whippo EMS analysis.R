########################
### Epifaunal diversity patterns among eelgrass meadows suggest metacommunity structure
### R. Whippo, N. Knight, C. Prentice, J. Cristiani, M. Siegle, M. I. O'Connor
### Analysis: elements of metacommunity structure (EMS)
### Code written by Mary O'Connor (with help from Domonik Bahlburg)
### started March 2016
### modified file in Jan 2017 to streamline with 'Whippo diversity analysis.R' file
#########################



# Load libraries ----------------------------------------------------------

library(metacom)
library(broom)

# Data preparation --------------------------------------------------------

## start with data2 file:
## create site-level data by collapsing across plots; these names refer to datafiles created in the 'Whippo diversity analysis.R' file.
start.data <- dataAUG # dataMAY, data.mp, dataAUG, dataJULY, data3times, dataJULY9
data.ms <- ddply(start.data, .(TimeID, area, species), summarise, sum(abundance))
data2 <- dcast(data.ms, TimeID ~ species, mean) #order

## clean out any empty species or empty sites
rowSums(data2[-1]) -> data2$totals
data2 <- data2[which(data2$totals != "0"),] 
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
Metacommunity(data4, verbose = TRUE, order = FALSE, method = "r1") -> meta
meta[2:4]

a <- as.data.frame(meta[1])

pdf('July 9 sites epifauna.pdf', width = 7, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="July 9 epifauna YES",  border="black", scales = list(cex = c(0.4, 0.4), x = list(rot = c(90))))
dev.off()

IdentifyStructure(meta)
Imagine(data4, sitenames = rownames(data4), speciesnames = colnames(data4), xline = 2)
axis(1, las = 2, cex = 1)
