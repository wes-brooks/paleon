#Required libraries:
require(bitops, lib.loc='R-libs')
require(RCurl, lib.loc='R-libs')

load_https <- function(url, sep=',', header=TRUE, row.names=NULL, ...) {
  # Import the data:
  read.table(text = getURL(url,
    followlocation=TRUE, cainfo=system.file("CurlSSL", "cacert.pem", package="RCurl")),
    sep=sep, header=header, row.names=row.names, ...)
}

#Import the stem density, count, and standard deviation data
biomass <- load_https("https://raw.github.com/wesesque/paleon/master/data/biomass.tablev1_3alb.csv")
biomass.wi <- load_https("https://raw.github.com/wesesque/paleon/master/data/biomassbysp_1_3_albwisc.csv")

#Process the Wisconsin biomass data
p0 = ncol(biomass.wi)
biomass.wi$tot = apply(biomass.wi[,4:p0], 1, sum)
biomass.wi = biomass.wi[which(!is.na(biomass.wi[,4])),]

#Process biomass data
p2 = ncol(biomass)
biomass$tot = apply(biomass[,4:p2], 1, sum)
biomass = biomass[which(!is.na(biomass$tot)),]

#Divide into hardwoods and softwoods:
hard = c("Ashes", "Birches", "Elms", "Maple", "Poplar", "Basswood", "Oaks", "Willow", "Alder", "Ironwoods", "Walnuts", "Hickory", "Beech", "Celtis", "Cherries", "Juniper", "Rose.Trees", "Sycamore", "Hardwood.undif", "Water", "Cornus", "Buckeye")
soft = c("Tamarack", "Pine", "Fir", "Cedar", "Spruce", "Hemlock")

biomass$hardwood = rowSums(biomass[,hard])
biomass$softwood = rowSums(biomass[,soft])

biomass.wi$hardwood = rowSums(biomass.wi[,hard[-which(hard=='Buckeye')]])
biomass.wi$softwood = rowSums(biomass.wi[,soft])

#Make sure changes are reflected in the global environment:
assign('biomass', biomass, envir=.GlobalEnv)
assign('biomass.wi', biomass.wi, envir=.GlobalEnv)