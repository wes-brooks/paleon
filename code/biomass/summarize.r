require(RCurl)

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

#Divide into hardwoods and softwoods:
hard = c("Ashes", "Birches", "Elms", "Maple", "Poplar", "Basswood", "Oaks", "Willow", "Alder", "Ironwoods", "Walnuts", "Hickory", "Beech", "Celtis", "Cherries", "Juniper", "Rose.Trees", "Sycamore", "Hardwood.undif", "Water", "Cornus", "Buckeye")
soft = c("Tamarack", "Pine", "Fir", "Cedar", "Spruce", "Hemlock")

biomass.wi$hardwood = rowSums(biomass.wi[,hard[-which(hard=='Buckeye')]])
biomass.wi$softwood = rowSums(biomass.wi[,soft])


hardwood = read.csv("~/misc/paleon/logbiomass-NA-hardwood.csv", header=FALSE)
softwood = read.csv("~/misc/paleon/logbiomass-NA-softwood.csv", header=FALSE)

rowSums(exp(softwood)) -> s
rowSums(exp(hardwood)) -> h

pdf("~/git/paleon/figures/biomass/hardwood-histogram.pdf", height=7, width=7)
    hist(h, breaks=50, xlab='biomass', main='') #, main="Histogram of total (smoothed)\nhardwood biomass in Wisconsin", xlab="biomass")
    abline(v=sum(biomass.wi$hardwood), col='red', lty=2, lwd=2)
dev.off()


pdf("~/git/paleon/figures/biomass/softwood-histogram.pdf", height=7, width=7)
    hist(s, breaks=50, xlab='biomass', main='') #main="Histogram of total (smoothed)\nsoftwood biomass in Wisconsin", xlab="biomass")
    abline(v=sum(biomass.wi$softwood), col='red', lty=2, lwd=2)
dev.off()


> dev.new()
> hist(b, breaks=50)
> abline(v=sum(biomass.wi$hardwood), col='red', lty=2, lwd=2)