#Required libraries:
#require(RCurl)
#require(devtools)
require(evaluate, lib.loc='R-libs')
require(bitops, lib.loc='R-libs')
require(RCurl, lib.loc='R-libs')
require(devtools, lib.loc='R-libs')

#prerequisites for the brooks package
require(plotrix, lib.loc='R-libs')
require(gtable, lib.loc='R-libs')
require(proto, lib.loc='R-libs')
require(plyr, lib.loc='R-libs')
require(reshape2, lib.loc='R-libs')
require(scales, lib.loc='R-libs')
require(ggplot2, lib.loc='R-libs')
require(RColorBrewer, lib.loc='R-libs')
require(xtable, lib.loc='R-libs')
#If the 'brooks' package isnt loaded then import it from github:
if (!'package:brooks' %in% search()) {
    install_github('brooks', 'wesesque', quick=TRUE)
    require(brooks)
}

#Import the stem density, count, and standard deviation data
composition <- load_https("https://raw.github.com/wesesque/paleon/master/data/glo.forest.composition_v1_3alb.csv")
biomass <- load_https("https://raw.github.com/wesesque/paleon/master/data/biomass.tablev1_3alb.csv")
biomass.wi <- load_https("https://raw.github.com/wesesque/paleon/master/data/biomassbysp_1_3_albwisc.csv")
stemdens <- load_https("https://raw.github.com/wesesque/paleon/master/data/density.ba.biomass_v1_3alb.csv")

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

#Divide into hardwoods and softwoods:
hardwood = c("Ashes", "Birches", "Elms", "Maple", "Poplar", "Basswood", "Oaks", "Willow", "Alder", "Ironwoods", "Walnuts", "Hickory", "Beech", "Celtis", "Cherries", "Juniper", "Rose.Trees", "Sycamore", "Hardwood.undif", "Water", "Cornus", "Buckeye")
softwood = c("Tamarack", "Pine", "Fir", "Cedar", "Spruce", "Hemlock")

biomass$hardwood = rowSums(biomass[,hardwood])
biomass$softwood = rowSums(biomass[,softwood])

biomass.wi$hardwood = rowSums(biomass.wi[,hardwood[-which(hardwood=='Buckeye')]])
biomass.wi$softwood = rowSums(biomass.wi[,softwood])


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


#Make sure changes are reflected in the global environment:
assign('biomass', biomass, envir=.GlobalEnv)
assign('biomass.wi', biomass.wi, envir=.GlobalEnv)