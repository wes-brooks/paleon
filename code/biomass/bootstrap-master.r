library(RCurl)
library(devtools)

#These functions import code and data directly from the github servers:
source("../brooks/R/load_https.r")

#Import the data and the heatmap code.
source_url("https://raw.github.com/wesesque/paleon/master/code/biomass/biomass-import.r")
taxa = c('tot', 'Cherries', 'Willow', 'Walnuts', 'Hickory', 'Beech', 'Fir', 'Spruce', 'Ironwoods', 'Cedar', 'Hemlock', 'Basswood', 'Ashes', 'Elms', 'Poplar', 'Pine', 'Tamarack', 'Birches', 'Maple', 'Oaks')

args <- commandArgs(trailingOnly = TRUE)
cluster = as.numeric(args[1])
indx = as.numeric(args[2]) + 1
taxon = taxa[indx]

source_url("https://raw.github.com/wesesque/paleon/master/code/biomass/bootstrap.r")
