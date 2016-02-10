# eelgrass_model
This code is associated with the manuscript of Ross Whippo's masters thesis. 

Temporal variation exceeds spatial variation in grazer biodiversity in eelgrass meadows
 
R. Whippo1,2, N.S. Knight1, C. Prentice3, M. Siegle1, M. I. O’Connor1
 
 1. 	Department of Zoology and Biodiversity Research Centre, University of British Columbia, Vancouver, B.C. V6T 1Z4
2. 	Smithsonian Institution, Tennenbaum Marine Observatories Network, PO Box 37012, National Museum of Natural History, MRC 106, Washington, DC 20013-7012
3. 	Carolyn’s current address: SFU… 

### Data files to use: (summary of all data files below)
rawcomm.csv  # master data file with species abundances at each site.
site.info.csv # meadow area, order, distance from fresh water
grazertraits.csv # information about which species are grazers.


### Code files to use:
Whippo diversity analysis master.R starts with rawcomm data and plots basic diversity estimates from each plot for each meadow. 


### data files - summary
data1 <- read.csv('comm_data.csv') # appears to be all species abundance data without site or time lables (useless?)

data2 <- read.csv('comm_data_main.csv') # appears to be all species abundance data for all sites without site or time lables (useless?)

data3 <- read.csv('plot_data.csv') # summaries of plot-scale diversity indices

data4 <- read.csv('plot_data_cnt.csv') # not sure how this differs from data3

data5 <- read.csv('plot_data_main.csv') # not sure how this differs from data3

data6 <- read.csv('rawcomm.csv') # this appears to be the master data file

data7 <- read.csv('satellite_site_totals.csv') # species abundance data totals for satellite sites

data8 <- read.csv('site.info.csv') # meadow area, order, distance from fresh water

data9 <- read.csv('site_avgs.csv') # diversity stats (alpha, beta, etc) for all sites

data10 <- read.csv('site_avgs_main.csv') # diversity stats (alpha, beta, etc) for main sites

data11 <- read.csv('site_ID.csv')

data12 <- read.csv('site_totals.csv') # specis abundances w/ no site codes.
