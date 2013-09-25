library(RCurl, lib.loc=c('R', 'R-libs'))

#These functions import code and data directly from the github servers:
source("../brooks/code/source_https.r")
source("../brooks/code/load_https.r")

#Import the data and the heatmap code.
source_https("https://raw.github.com/wesesque/paleon/master/code/biomass-import.r")
taxa = c('Cherries', 'Willow', 'Walnuts', 'Hickory', 'Beech', 'Fir', 'Spruce', 'Ironwoods', 'Cedar', 'Hemlock', 'Basswood', 'Ashes', 'Elms', 'Poplar', 'Pine', 'Tamarack', 'Birches', 'Maple', 'Oaks')
args <- commandArgs(trailingOnly = TRUE)
indx = as.numeric(args[2])
sp = taxa[indx]

source_https("https://raw.github.com/wesesque/paleon/master/code/biomass/bootstrap.r")
