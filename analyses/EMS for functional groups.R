########################
### Whippo et al eelgrass epifaunal data for Barkley Sound
### elements of metacommunity structure (EMS)
### Code written by Mary O'Connor (with help from Domonik Bahlburg)
### started March 2016
#########################

### code for functional groups; start w EMS code

head(data.tr)

# choose a functional group (use subset, or it causes trouble later)
#grazers <- data.tr[(data.tr$function.=="grazer"), c(1:2,4:12)]
grazers <- subset(data.tr, function. == "grazer", select = c(1:2,4:12))
crust <- subset(data.tr, group == "crustacean", select = c(1:2,4:12))
gast <- subset(data.tr, group == "gastropod", select = c(1:2,4:12))


## choose a sampling time
data.ft <- poly #gast #crust #grazers
#data3timesg <- data.ft[(data.ft$site!="BE" & data.ft$site!="EI" & data.ft$site!="CC" & data.ft$site!="BI"),]
dataJULYg <- data.ft[(data.ft$Time.Code2=="C"),]

## 1. create site-level data by collapsing across plots
start.data <- dataJULYg # dataMAY, data.mp, dataAUG, dataJULY, data3times
data.ms <- ddply(start.data, .(TimeID, order, species), summarise, sum(abundance))
#data.ms <- data.ms[-nrow(data.ms),]
data2 <- dcast(data.ms, TimeID + order ~ species, mean)

## order file to impose environmental gradient, if desired
data2 <- data2[order(data2$order),]
data2 <- data2[,-2]

## clean out any empty species or empty sites
rowSums(data2[-1]) -> data2$totals
data2 <- data2[which(data2$totals != "0"),] 
#data3 <- data2[-which(is.na(data2$totals)),] 
rownames(data2) <- data2[,1]
data3 <- data2[, -c(1, ncol(data2))]
data3[data3 > 0] <- 1
cols.to.delete <- which(colSums(data3) == '0')
data4 <- data3[, -(cols.to.delete)]

### run metacommunity analysis 
Metacommunity(data4, verbose = TRUE, allowEmpty = FALSE, order = TRUE) -> meta.p
meta.p[2:4]

a <- as.data.frame(meta.c[1])

pdf('July crust dfw.pdf', width = 7, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="July crust dfw",  border="black", scales = list(cex = c(0.4, 0.4), x = list(rot = c(90))))
dev.off()
