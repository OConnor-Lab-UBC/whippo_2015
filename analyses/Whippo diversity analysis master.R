##### CODE FOR WORK FLOW FROM RAW DATA DO DIVERISTY METRICS #####
##### MARY O'CONNOR SEPTEMBER 2015
##### For Whippo et al, seagrass biodiversity paper


library(vegan)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

#### FUNCTIONS #####

idmaker = function(x) { return(paste(x, collapse=""))}

count.func <- function(x) {
  y = ifelse( x > 0, 1, 0)
  sum(y)
}

#### multiplot function for later: (from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#####


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

#### ALL TAXA ####

## a) 
## create an id for each sample at each site, time. 
data$ID <- apply(data[,c("site", "Time.Code2", "Sample")], 1, idmaker)

## now treat each sample as a group & collapse across sizes to sum taxa
data.s <- group_by(data, ID)
func <- function(x) { sum(x) }
data.sum <- as.data.frame(summarise_each(data.s, funs(sum), (Idotea.resecata:Alia.carinata)))

## add in sample information 
data.desc <- ddply(data.s, .(site, Time.Code2, Sample), summarise, unique(ID)) ## this step could be consolidated later.
data.sum2 <- merge(data.sum, data.desc, by.x = "ID", by.y = "unique(ID)", all.x = TRUE, all.y = FALSE)

## add in site information
data.sum3 <- merge(data.sum2, sites, by.x = 'site', by.y = 'site')

data.sum3$H <- diversity(data.sum2[,2:47], index = "shannon") 
data.sum3$S <- diversity(data.sum2[,2:47], index = "simpson")
data.sum3$N <- apply(data.sum2[,2:47], 1, function(x) sum(x))
data.sum3$R <- apply(data.sum2[,2:47], 1, function(x) count.func(x))
head(data.sum3)

totN <- ggplot(data.sum3, aes(factor(order), log(N)))
totN_abund <- totN + geom_boxplot(aes(fill = factor(Time.Code2))) 

totS <- ggplot(data.sum3, aes(factor(order), S))
totS_plot <- totS + geom_boxplot(aes(fill = factor(Time.Code2))) # no clear signal of time, or space

totH <- ggplot(data.sum3, aes(factor(order), H))
totH_plot <- totH + geom_boxplot(aes(fill = factor(Time.Code2))) # no clear signal of time, or space

totR <- ggplot(data.sum3, aes(factor(order), R))
totR_plot <- totR + geom_boxplot(aes(fill = factor(Time.Code2)))

multiplot(totN_abund, totR_plot, totS_plot, totH_plot)


### GRAZERS ONLY

### b) merge the rawcomm with the grazer traits. 
### first, melt the dataframe to make species into a column so I can merge with species traits file.
data3 <- melt(data.sum, id.vars = "ID", variable.name = 'species')
data3.1 <- merge(data3, traits, by.x = 'species', by.y = 'species.names', all = TRUE)

## remove duplicate species name column
data3.1 <- data3.1[,-4]
## remove other columns we don't need now
data3.1 <- data3.1[,-(8:12)]

## subset out grazers and estimate diversity
data.gr <- filter(data3.1, function. == "grazer")
data.grz <- dcast(data.gr, ID ~ species)

data.grz$H <- diversity(data.grz[,2:18], index = "shannon") 
data.grz$S <- diversity(data.grz[,2:18], index = "simpson")
data.grz$N <- apply(data.grz[,2:18], 1, function(x) sum(x))
data.grz$R <- apply(data.grz[,2:18], 1, function(x) count.func(x))
head(data.grz)

## add in sample information 
data.grz2 <- merge(data.grz, data.desc, by.x = "ID", by.y = "..1", all.x = TRUE, all.y = FALSE)

## add in site information
data.grz3 <- merge(data.grz2, sites, by.x = 'site', by.y = 'site')

## look at it
plot(log(data.grz3$N) ~ data.grz3$order*data.grz3$Time.Code2)

grzN <- ggplot(data.grz3, aes(factor(order), log(N)))
grzN_abund <- grzN + geom_boxplot(aes(fill = factor(Time.Code2))) 

grzS <- ggplot(data.grz3, aes(factor(order), S))
grzS_plot <- grzS + geom_boxplot(aes(fill = factor(Time.Code2))) 

grzH <- ggplot(data.grz3, aes(factor(order), H))
grzH_plot <- grzH + geom_boxplot(aes(fill = factor(Time.Code2))) 

grzR <- ggplot(data.grz3, aes(factor(order), R))
grzR_plot <- grzR + geom_boxplot(aes(fill = factor(Time.Code2))) 


multiplot(grzN_abund, grzR_plot, grzS_plot, grzH_plot)


##### exclude 1 mm sizes ####

## remove the 1mm rows
data1 <- data[(data$Sieve!='1'),]

## create an id for each sample at each site, time. 
idmaker = function(x) { return(paste(x, collapse=""))}
data1$ID <- apply(data1[,c("site", "Time.Code2", "Sample")], 1, idmaker)

## now treat each sample as a group collapse across sizes to sum taxa
data1.s <- group_by(data1, ID)
func <- function(x) { sum(x) }
data1.sum <- as.data.frame(summarise_each(data1.s, funs(sum), (Idotea.resecata:Alia.carinata)))

data1.sum$H <- diversity(data1.sum[,2:47], index = "shannon") 
data1.sum$S <- diversity(data1.sum[,2:47], index = "simpson")
data1.sum$N <- apply(data1.sum[,2:47], 1, function(x) sum(x))
data1.sum$R <- apply(data1.sum[,2:47], 1, function(x) count.func(x))
head(data1.sum)

## add in sample information 
data1.desc <- ddply(data1.s, .(site, Time.Code2, Sample), summarise, unique(ID))
data1.sum2 <- merge(data1.sum, data1.desc, by.x = "ID", by.y = "..1", all.x = TRUE, all.y = FALSE)

## add in site information
data1.sum3 <- merge(data1.sum2, sites, by.x = 'site', by.y = 'site')

## look at it
plot(log(data1.sum3$N) ~ data1.sum3$order) # spatial trend is still there
plot(log(data1.sum3$N) ~ data1.sum3$Time.Code2) # temporal trend is weaker

bigN <- ggplot(data1.sum3, aes(factor(order), log(N)))
bigN_abund <- bigN + geom_boxplot(aes(fill = factor(Time.Code2))) ## sptial trend there, temporal trend there but less dramatic 

bigS <- ggplot(data1.sum3, aes(factor(order), S))
bigS_plot <- bigS + geom_boxplot(aes(fill = factor(Time.Code2))) # no clear signal of time, or space

bigH <- ggplot(data1.sum3, aes(factor(order), H))
bigH_plot <- bigH + geom_boxplot(aes(fill = factor(Time.Code2))) # no clear signal of time, or space

bigR <- ggplot(data1.sum3, aes(factor(order), R))
bigR_plot <- bigR + geom_boxplot(aes(fill = factor(Time.Code2)))

multiplot(bigN_abund, bigR_plot, bigS_plot, bigH_plot)

##### only 1-2 mm sizes ####

## remove the 1mm rows
data.sm <- data[(data$Sieve=='1'),]

## create an id for each sample at each site, time. 
idmaker = function(x) { return(paste(x, collapse=""))}
data.sm$ID <- apply(data.sm[,c("site", "Time.Code2", "Sample")], 1, idmaker)

## now treat each sample as a group collapse across sizes to sum taxa
data.sm.s <- group_by(data.sm, ID)
func <- function(x) { sum(x) }
data.sm.sum <- as.data.frame(summarise_each(data.sm.s, funs(sum), (Idotea.resecata:Alia.carinata)))

data.sm.sum$H <- diversity(data.sm.sum[,2:47], index = "shannon") # interesting; works here but didn't with grazers
data.sm.sum$S <- diversity(data.sm.sum[,2:47], index = "simpson")
data.sm.sum$N <- apply(data.sm.sum[,2:47], 1, function(x) sum(x))
data.sm.sum$R <- apply(data.sm.sum[,2:47], 1, function(x) count.func(x))
head(data.sm.sum)

## add in sample information 
data.sm.desc <- ddply(data.sm.s, .(site, Time.Code2, Sample), summarise, unique(ID))
data.sm.sum2 <- merge(data.sm.sum, data.sm.desc, by.x = "ID", by.y = "..1", all.x = TRUE, all.y = FALSE)

## add in site information
data.sm.sum3 <- merge(data.sm.sum2, sites, by.x = 'site', by.y = 'site')

## look at it
plot(log(data.sm.sum3$N) ~ data.sm.sum3$order)
plot(log(data.sm.sum3$N) ~ data.sm.sum3$Time.Code2)

smN <- ggplot(data.sm.sum3, aes(factor(order), log(N)))
smN_abund <- smN + geom_boxplot(aes(fill = factor(Time.Code2))) ## clear signal of more big things after May at marine sites. This does not translate to a change in diversity; S and H metrics show trends decoupled from abundance, suggesting a role of changes in dominance with time as well. 

smS <- ggplot(data.sm.sum3, aes(factor(order), S))
smS_plot <- smS + geom_boxplot(aes(fill = factor(Time.Code2))) # no clear signal of time, or space

smH <- ggplot(data.sm.sum3, aes(factor(order), H))
smH_plot <- smH + geom_boxplot(aes(fill = factor(Time.Code2))) # no clear signal of time, or space

smR <- ggplot(data.sm.sum3, aes(factor(order), R))
smR_plot <- smR + geom_boxplot(aes(fill = factor(Time.Code2)))

multiplot(smN_abund, smR_plot, smS_plot, smH_plot)


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





### modeling 
library(lme4)
library(MuMIn)

## total inverts
class(data.sum3$order)
mem.a1 <- lmer(log(N+0.01) ~ Time.Code2 + (Time.Code2|site), data = data.sum3, REML = FALSE)
mem.a2 <- lmer(log(N+0.01) ~ Time.Code2 + dfw + (Time.Code2|site), data = data.sum3, REML = FALSE)
mem.a3 <- lmer(log(N+0.01) ~ Time.Code2*dfw + (Time.Code2|site), data = data.sum3, REML = FALSE)
mem.a3.norand <- lmer(log(N+0.01) ~ Time.Code2*dfw + (1|site), data = data.sum3, REML = FALSE)

abund.sel <- model.sel(mem.a1,mem.a2,mem.a3)
anova(mem.a3, mem.a3.norand)
anova(mem.a3, mem.a2)

summary(mem.a2)

mod1 <- lm(log(N+0.01) ~ Time.Code2 + dfw, data = data.sum3)
