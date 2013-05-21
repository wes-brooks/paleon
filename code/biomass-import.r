setwd("~/git/paleon")

#Import the stem density, count, and standard deviation data
composition <- read.csv("data/glo.forest.composition_v1_3alb.csv", header=T)
biomass <- read.csv("data/biomass.tablev1_3alb.csv", header=T)
biomass.wi <- read.csv("data/biomassbysp_1_3_albwisc.csv", header=T)
stemdens <- read.csv("data/density.ba.biomass_v1_3alb.csv", header=T)

#Process the Wisconsin biomass data
p0 = ncol(biomass.wi)
biomass.wi$tot = apply(biomass.wi[,4:p0], 1, sum)
biomass.wi = biomass.wi[which(!is.na(biomass.wi[,4])),]

#Process composition data
p = ncol(composition)
composition$tot = apply(composition[,4:p], 1, sum)
composition[,4:p] = composition[,4:p] / composition$tot
composition = composition[which(!is.na(composition$tot)),]

#Process biomass data
p2 = ncol(biomass)
biomass$tot = apply(biomass[,4:p2], 1, sum)
biomass = biomass[which(!is.na(biomass$tot)),]

#Extract the WI composition data
indx = which( outer(biomass.wi$x, composition$x, "==") & 
       outer(biomass.wi$y, composition$y, "=="), 
       arr.ind=TRUE)
composition.wi = composition[indx[,2],]

#Extract the WI stemdensity data
indx = which( outer(biomass.wi$x, stemdens$x, "==") & 
       outer(biomass.wi$y, stemdens$y, "=="), 
       arr.ind=TRUE)
stemdens.wi = stemdens[indx[,2],]
