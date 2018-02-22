###################################################################################
#                                                                                ##
# Basic Diversity Metrics                                                        ##
# Data are current as of 2018-02-21                                              ##
# Data source: O'Connor Lab - UBC                                                ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2018-02-21                                                        ##
#                                                                                ##
###################################################################################

# SUMMARY: 
# Basic diversity analyses used to address reviewer comments of "Epifaunal 
# diversity patterns among eelgrass meadows suggest metacommunity structure" at
# Ecosphere including rarified species richness, ENS, Shannon Diversity, Simpson
# Evenness, and Hellinger Distance Matrix


# Required Files (check that script is loading latest version):
# epicomm_201802.csv

# Associated Scripts:
# NONE

# TO DO

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# RECENT CHANGES TO SCRIPT                                                        #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# MANIPULATE DATA                                                                 #   
#                                                                                 #
###################################################################################

###################################################################################
# RECENT CHANGES TO SCRIPT                                                        #
###################################################################################

# 2018-02-21 Created script. 

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

library(tidyverse) # manipulate data
library(reshape2) # manipulate data
library(vegan) # diversity analyses
library(viridis) # color palette
library(ggpubr) # combining plots

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# full community data
epicomm <- read_csv("./data/epicomm_201802.csv")

# remove redundant columns and summarise species occurrences
epicomm_summ <- epicomm %>%
  select(-date, -Sieve) %>%
  unite(sitetime, site, 'Time Code', sep = "", remove = TRUE)
# fix CBC typo
epicomm_summ$sitetime <- gsub("CBC ", "CBC", epicomm_summ$sitetime)
# finish summarization for every sample by summing
epicomm_summ <- epicomm_summ %>%
  group_by(sitetime, Sample) %>%
  summarise_all(funs(sum))

# group all primary sites 
primary <- c("DCA", "DCC", "DCE", "WIA", "WIC", "WIE", "RPA", "RPC", "RPE", "NBA", "NBC", "NBE", "CBA", "CBC", "CBE")
epicomm_prim_full <- epicomm_summ %>%
  select(-Sample) %>%
  subset(sitetime %in% primary) 

# separate each sitetime for analysis

DCA <- epicomm_prim_full %>%
  subset(sitetime == "DCA")
DCC <- epicomm_prim_full %>%
  subset(sitetime == "DCC")
DCE <- epicomm_prim_full %>%
  subset(sitetime == "DCE")

WIA <- epicomm_prim_full %>%
  subset(sitetime == "WIA")
WIC <- epicomm_prim_full %>%
  subset(sitetime == "WIC")
WIE <- epicomm_prim_full %>%
  subset(sitetime == "WIE")

RPA <- epicomm_prim_full %>%
  subset(sitetime == "RPA")
RPC <- epicomm_prim_full %>%
  subset(sitetime == "RPC")
RPE <- epicomm_prim_full %>%
  subset(sitetime == "RPE")

NBA <- epicomm_prim_full %>%
  subset(sitetime == "NBA")
NBC <- epicomm_prim_full %>%
  subset(sitetime == "NBC")
NBE <- epicomm_prim_full %>%
  subset(sitetime == "NBE")

CBA <- epicomm_prim_full %>%
  subset(sitetime == "CBA")
CBC <- epicomm_prim_full %>%
  subset(sitetime == "CBC")
CBE <- epicomm_prim_full %>%
  subset(sitetime == "CBE")

###################################################################################
# DIVERSITY MEASURES                                                              #
###################################################################################

############### OBSERVED RICHNESS

DCA_rich <- DCA[,2:31] %>%
  specnumber()
DCC_rich <- DCC[,2:31] %>%
  specnumber()
DCE_rich <- DCE[,2:31] %>%
  specnumber()

WIA_rich <- WIA[,2:31] %>%
  specnumber()
WIC_rich <- WIC[,2:31] %>%
  specnumber()
WIE_rich <- WIE[,2:31] %>%
  specnumber()

RPA_rich <- RPA[,2:31] %>%
  specnumber()
RPC_rich <- RPC[,2:31] %>%
  specnumber()
RPE_rich <- RPE[,2:31] %>%
  specnumber()

NBA_rich <- NBA[,2:31] %>%
  specnumber()
NBC_rich <- NBC[,2:31] %>%
  specnumber()
NBE_rich <- NBE[,2:31] %>%
  specnumber()

CBA_rich <- CBA[,2:31] %>%
  specnumber()
CBC_rich <- CBC[,2:31] %>%
  specnumber()
CBE_rich <- CBE[,2:31] %>%
  specnumber()

# combine into single data frame
# make all same length
length(DCA_rich) <- 17                      
length(DCC_rich) <- 17  
length(DCE_rich) <- 17 
length(WIA_rich) <- 17 
length(WIC_rich) <- 17  
length(WIE_rich) <- 17  
length(RPA_rich) <- 17  
length(RPC_rich) <- 17  
length(RPE_rich) <- 17  
length(NBA_rich) <- 17  
length(NBC_rich) <- 17  
length(NBE_rich) <- 17  
length(CBA_rich) <- 17  
length(CBC_rich) <- 17  
length(CBE_rich) <- 17 
# combine
rich_prim <- melt(data.frame(DCA_rich, DCC_rich, DCE_rich, WIA_rich, WIC_rich, WIE_rich, RPA_rich, RPC_rich, RPE_rich, NBA_rich, NBC_rich, NBE_rich, CBA_rich, CBC_rich, CBE_rich))
# rename columns, reduce sitetime values, and split into site and time
colnames(rich_prim) <- c('sitetime', 'value') 
rich_prim$sitetime <- as.character(rich_prim$sitetime)
rich_prim$sitetime <- substr(rich_prim$sitetime,1,3)
rich_prim <- transform(rich_prim, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
rich_prim <- rich_prim %>%
  select(-sitetime)
# reorder factors
rich_prim$site <- factor(rich_prim$site, levels = c("DC", "WI", "RP", "NB", "CB"))


############### ENS

DCA_ens <- exp(diversity(DCA[,2:31], "shannon"))
DCC_ens <- exp(diversity(DCC[,2:31], "shannon"))
DCE_ens <- exp(diversity(DCE[,2:31], "shannon"))

WIA_ens <- exp(diversity(WIA[,2:31], "shannon"))
WIC_ens <- exp(diversity(WIC[,2:31], "shannon"))
WIE_ens <- exp(diversity(WIE[,2:31], "shannon"))

RPA_ens <- exp(diversity(RPA[,2:31], "shannon"))
RPC_ens <- exp(diversity(RPC[,2:31], "shannon"))
RPE_ens <- exp(diversity(RPE[,2:31], "shannon"))

NBA_ens <- exp(diversity(NBA[,2:31], "shannon"))
NBC_ens <- exp(diversity(NBC[,2:31], "shannon"))
NBE_ens <- exp(diversity(NBE[,2:31], "shannon"))

CBA_ens <- exp(diversity(CBA[,2:31], "shannon"))
CBC_ens <- exp(diversity(CBC[,2:31], "shannon"))
CBE_ens <- exp(diversity(CBE[,2:31], "shannon"))

# combine into single data frame
# make all same length
length(DCA_ens) <- 17                      
length(DCC_ens) <- 17  
length(DCE_ens) <- 17 
length(WIA_ens) <- 17 
length(WIC_ens) <- 17  
length(WIE_ens) <- 17  
length(RPA_ens) <- 17  
length(RPC_ens) <- 17  
length(RPE_ens) <- 17  
length(NBA_ens) <- 17  
length(NBC_ens) <- 17  
length(NBE_ens) <- 17  
length(CBA_ens) <- 17  
length(CBC_ens) <- 17  
length(CBE_ens) <- 17 
# combine
ens_prim <- melt(data.frame(DCA_ens, DCC_ens, DCE_ens, WIA_ens, WIC_ens, WIE_ens, RPA_ens, RPC_ens, RPE_ens, NBA_ens, NBC_ens, NBE_ens, CBA_ens, CBC_ens, CBE_ens))
# rename columns, reduce sitetime values, and split into site and time
colnames(ens_prim) <- c('sitetime', 'value') 
ens_prim$sitetime <- as.character(ens_prim$sitetime)
ens_prim$sitetime <- substr(ens_prim$sitetime,1,3)
ens_prim <- transform(ens_prim, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
ens_prim <- ens_prim %>%
  select(-sitetime)
# reorder factors
ens_prim$site <- factor(ens_prim$site, levels = c("DC", "WI", "RP", "NB", "CB"))


############### SHANNON DIVERSITY

DCA_shannon <- DCA[,2:31] %>%
  diversity("shannon")
DCC_shannon <- DCC[,2:31] %>%
  diversity("shannon")
DCE_shannon <- DCE[,2:31] %>%
  diversity("shannon")

WIA_shannon <- WIA[,2:31] %>%
  diversity("shannon")
WIC_shannon <- WIC[,2:31] %>%
  diversity("shannon")
WIE_shannon <- WIE[,2:31] %>%
  diversity("shannon")

RPA_shannon <- RPA[,2:31] %>%
  diversity("shannon")
RPC_shannon <- RPC[,2:31] %>%
  diversity("shannon")
RPE_shannon <- RPE[,2:31] %>%
  diversity("shannon")

NBA_shannon <- NBA[,2:31] %>%
  diversity("shannon")
NBC_shannon <- NBC[,2:31] %>%
  diversity("shannon")
NBE_shannon <- NBE[,2:31] %>%
  diversity("shannon")

CBA_shannon <- CBA[,2:31] %>%
  diversity("shannon")
CBC_shannon <- CBC[,2:31] %>%
  diversity("shannon")
CBE_shannon <- CBE[,2:31] %>%
  diversity("shannon")

# combine into single data frame
# make all same length
length(DCA_shannon) <- 17                      
length(DCC_shannon) <- 17  
length(DCE_shannon) <- 17 
length(WIA_shannon) <- 17 
length(WIC_shannon) <- 17  
length(WIE_shannon) <- 17  
length(RPA_shannon) <- 17  
length(RPC_shannon) <- 17  
length(RPE_shannon) <- 17  
length(NBA_shannon) <- 17  
length(NBC_shannon) <- 17  
length(NBE_shannon) <- 17  
length(CBA_shannon) <- 17  
length(CBC_shannon) <- 17  
length(CBE_shannon) <- 17 
# combine
shannon_prim <- melt(data.frame(DCA_shannon, DCC_shannon, DCE_shannon, WIA_shannon, WIC_shannon, WIE_shannon, RPA_shannon, RPC_shannon, RPE_shannon, NBA_shannon, NBC_shannon, NBE_shannon, CBA_shannon, CBC_shannon, CBE_shannon))
# rename columns, reduce sitetime values, and split into site and time
colnames(shannon_prim) <- c('sitetime', 'value') 
shannon_prim$sitetime <- as.character(shannon_prim$sitetime)
shannon_prim$sitetime <- substr(shannon_prim$sitetime,1,3)
shannon_prim <- transform(shannon_prim, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
shannon_prim <- shannon_prim %>%
  select(-sitetime)
# reorder factors
shannon_prim$site <- factor(shannon_prim$site, levels = c("DC", "WI", "RP", "NB", "CB"))


############### EVENNESS

DCA_even <- diversity(DCA[,2:31])/log(specnumber(DCA[,2:31]))
DCC_even <- diversity(DCC[,2:31])/log(specnumber(DCC[,2:31]))
DCE_even <- diversity(DCE[,2:31])/log(specnumber(DCE[,2:31]))

WIA_even <- diversity(WIA[,2:31])/log(specnumber(WIA[,2:31]))
WIC_even <- diversity(WIC[,2:31])/log(specnumber(WIC[,2:31]))
WIE_even <- diversity(WIE[,2:31])/log(specnumber(WIE[,2:31]))

RPA_even <- diversity(RPA[,2:31])/log(specnumber(RPA[,2:31]))
RPC_even <- diversity(RPC[,2:31])/log(specnumber(RPC[,2:31]))
RPE_even <- diversity(RPE[,2:31])/log(specnumber(RPE[,2:31]))

NBA_even <- diversity(NBA[,2:31])/log(specnumber(NBA[,2:31]))
NBC_even <- diversity(NBC[,2:31])/log(specnumber(NBC[,2:31]))
NBE_even <- diversity(NBE[,2:31])/log(specnumber(NBE[,2:31]))

CBA_even <- diversity(CBA[,2:31])/log(specnumber(CBA[,2:31]))
CBC_even <- diversity(CBC[,2:31])/log(specnumber(CBC[,2:31]))
CBE_even <- diversity(CBE[,2:31])/log(specnumber(CBE[,2:31]))

# combine into single data frame
# make all same length
length(DCA_even) <- 17                      
length(DCC_even) <- 17  
length(DCE_even) <- 17 
length(WIA_even) <- 17 
length(WIC_even) <- 17  
length(WIE_even) <- 17  
length(RPA_even) <- 17  
length(RPC_even) <- 17  
length(RPE_even) <- 17  
length(NBA_even) <- 17  
length(NBC_even) <- 17  
length(NBE_even) <- 17  
length(CBA_even) <- 17  
length(CBC_even) <- 17  
length(CBE_even) <- 17 
# combine
even_prim <- melt(data.frame(DCA_even, DCC_even, DCE_even, WIA_even, WIC_even, WIE_even, RPA_even, RPC_even, RPE_even, NBA_even, NBC_even, NBE_even, CBA_even, CBC_even, CBE_even))
# rename columns, reduce sitetime values, and split into site and time
colnames(even_prim) <- c('sitetime', 'value') 
even_prim$sitetime <- as.character(even_prim$sitetime)
even_prim$sitetime <- substr(even_prim$sitetime,1,3)
even_prim <- transform(even_prim, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
even_prim <- even_prim %>%
  select(-sitetime)
# reorder factors
even_prim$site <- factor(even_prim$site, levels = c("DC", "WI", "RP", "NB", "CB"))




############### RARIFIED RICHNESS


DCA_rare <- DCA[,2:31] %>%
  rarefy(5)

DCA_rare <- DCA[,2:31] %>%
  rarefy(5)
DCC_rare <- DCC[,2:31] %>%
  rarefy(5)
DCE_rare <- DCE[,2:31] %>%
  rarefy(5)

WIA_rare <- WIA[,2:31] %>%
  rarefy(5)
WIC_rare <- WIC[,2:31] %>%
  rarefy(5)
WIE_rare <- WIE[,2:31] %>%
  rarefy(5)

RPA_rare <- RPA[,2:31] %>%
  rarefy(5)
RPC_rare <- RPC[,2:31] %>%
  rarefy(5)
RPE_rare <- RPE[,2:31] %>%
  rarefy(5)

NBA_rare <- NBA[,2:31] %>%
  rarefy(5)
NBC_rare <- NBC[,2:31] %>%
  rarefy(5)
NBE_rare <- NBE[,2:31] %>%
  rarefy(5)

CBA_rare <- CBA[,2:31] %>%
  rarefy(5)
CBC_rare <- CBC[,2:31] %>%
  rarefy(5)
CBE_rare <- CBE[,2:31] %>%
  rarefy(5)

# combine into single data frame
# make all same length
length(DCA_rare) <- 17                      
length(DCC_rare) <- 17  
length(DCE_rare) <- 17 
length(WIA_rare) <- 17 
length(WIC_rare) <- 17  
length(WIE_rare) <- 17  
length(RPA_rare) <- 17  
length(RPC_rare) <- 17  
length(RPE_rare) <- 17  
length(NBA_rare) <- 17  
length(NBC_rare) <- 17  
length(NBE_rare) <- 17  
length(CBA_rare) <- 17  
length(CBC_rare) <- 17  
length(CBE_rare) <- 17 
# combine
rare_prim <- melt(data.frame(DCA_rare, DCC_rare, DCE_rare, WIA_rare, WIC_rare, WIE_rare, RPA_rare, RPC_rare, RPE_rare, NBA_rare, NBC_rare, NBE_rare, CBA_rare, CBC_rare, CBE_rare))
# rename columns, reduce sitetime values, and split into site and time
colnames(rare_prim) <- c('sitetime', 'value') 
rare_prim$sitetime <- as.character(rare_prim$sitetime)
rare_prim$sitetime <- substr(rare_prim$sitetime,1,3)
rare_prim <- transform(rare_prim, site = substr(sitetime, 1, 2), time = substr(sitetime, 3, 3))
rare_prim <- rare_prim %>%
  select(-sitetime)
# reorder factors
rare_prim$site <- factor(rare_prim$site, levels = c("DC", "WI", "RP", "NB", "CB"))













###################################################################################
# FIGURES                                                                         #
###################################################################################

########### OBSERVED RICHNESS

rich_plot <- ggplot(rich_prim, aes(x = time, y = value, fill = site)) + 
  geom_boxplot() + 
  fill_palette(viridis(5, begin = 0.3)) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Richness")

########### ENS

ens_plot <- ggplot(ens_prim, aes(x = time, y = value, fill = site)) + 
  geom_boxplot() + 
  fill_palette(viridis(5, begin = 0.3)) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="ENS")

########### SHANNON DIVERSITY

shannon_plot <- ggplot(shannon_prim, aes(x = time, y = value, fill = site)) + 
  geom_boxplot() + 
  fill_palette(viridis(5, begin = 0.3)) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x=" May                  June/July                 August", y="Shannon Diversity")

########### EVENNESS

even_plot <- ggplot(even_prim, aes(x = time, y = value, fill = site)) + 
  geom_boxplot() + 
  fill_palette(viridis(5, begin = 0.3)) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x=" May                 June/July                 August", y="Evenness")

########### FULL FIGURE

Figure3 <- ggarrange(rich_plot, ens_plot, shannon_plot, even_plot,
                     labels = c("A", "B", "C", "D"),
                     ncol = 2, nrow = 2,
                     common.legend = TRUE, legend = "right")
#annotate_figure(Figure3, bottom = text_grob("Figure 3: Measures of A) observed richness, B) shannon diversity, and C) effective number of species (ENS) across five seagrass habitats types \n sampled in May, June/July, and August", size = 10))

# best size: ~950x620




#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#


##### SCRATCH PAD

# code to separate sitetime column


separate(shannon_prim, sitetime, 1, c("site", "time"), 2)


S <- specnumber(DCA[,2:31]) # observed number of species
(raremax <- min(rowSums(DCA[,2:31])))
Srare <- rarefy(DCA[,2:31], raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(DCA[,2:31], step = 20, sample = raremax, col = "blue", cex = 0.6)




