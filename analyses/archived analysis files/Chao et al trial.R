### testing out Chao et al 2014's methods
library(iNEXT)
library(ggplot2)

data(spider)
str(spider)
iNEXT(spider, q=0, datatype="abundance")


iNEXT(x, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50) -> test


## bring in Ross's data from two sites: 
RBE <- BEsp[,2]
RRP <- rankRP[,2]
RDC <- rankDC[,2]
RWI <- rankWI[,2]
RCC <- rankCC[,2]
RCB <- rankCB[,2]
RNB <- rankNB[,2]
REI <- rankEI[,2]
RBI <- rankBI[,2]

Rd <- cbind(RBE, RRP, RDC, RWI, RCC, RCB, RNB, REI, RBI)
Rd <- Rd[-(45:46),]

## based on a first run of iNEXT, I see that our maxn is 17484. So can I now set m? 
m <- seq(1,17484,17484/100)

iNEXT(Rd, q= c(0), datatype="abundance", size=m, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50) -> test
# goes faster if you just do q = 0

ggiNEXT(test, type=1, se=TRUE, facet.var="site", color.var="site", grey=FALSE)  

ggiNEXT(test, type=1, facet.var="none", color.var="site")

ggiNEXT(test, type=3, facet.var="site")

test

## you choose a base level that is 2x the smallest reference sample size, OR the largest reference sample size, whichever is biggest. so in this case teh base level for comparison is n = 17484. 
## then from iNEXT, we can (somehow) extract these richnesses and their confidence intervals. 
$RBE
test$iNextEst$RBE[103,4:6] -> BE
test$iNextEst$RRP[103,4:6] -> RP
test$iNextEst$RDC[103,4:6] -> DC
test$iNextEst$RWI[103,4:6] -> WI
test$iNextEst$RCC[103,4:6] -> CC
test$iNextEst$RCB[103,4:6] -> CB
test$iNextEst$RNB[103,4:6] -> NB
test$iNextEst$REI[103,4:6] -> EI
test$iNextEst$RBI[103,4:6] -> BI

q1 <- rbind(BE, RP, DC, WI, CC, CB, NB, EI, BI)
rownames(q1) <- c("BE", "RP", "DC", "WI", "CC", "CB", "NB", "EI", "BI")

plot(q1[,1], pch = 19, ylim = c(0, 40), xlim = c(0,10), axes = FALSE, ylab = 'Species Richness', xlab = "Site")
arrows(c(1:9), q1[,1], c(1:9), q1[,3], angle = 90, length = 0.05)
arrows(c(1:9), q1[,1], c(1:9), q1[,2], angle = 90, length = 0.05)
axis(1, at = c(1:9), labels = rownames(q1), pos = 0)
axis(2, las = 2, pos = 0)

q1
