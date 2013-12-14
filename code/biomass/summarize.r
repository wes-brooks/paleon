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




#Plot the mean and standard deviation of draws for hardwood, softwood biomass:
require(fBasics)
#hard = read.csv("/NO_BACKUP/brooks/paleon/paleon/output/logbiomass-NA-hardwood.csv", header=FALSE)
#soft = read.csv("/NO_BACKUP/brooks/paleon/paleon/output/logbiomass-NA-softwood.csv", header=FALSE)
hard = read.csv("~/Dropbox/logbiomass-WI-hardwood.csv", header=FALSE)
soft = read.csv("~/Dropbox/logbiomass-WI-softwood.csv", header=FALSE)
tot = read.csv("~/misc/paleon/logbiomass-NA-tot-wi.csv", header=FALSE)

hh = colMeans(exp(hard))
ss = colMeans(exp(soft))
tt = colMeans(exp(tot))

hh.sd = colStdevs(exp(hard))
ss.sd = colStdevs(exp(soft))
tt.sd = colStdevs(exp(tot))


#rm(hard)
#rm(soft)

#hard2 = cbind(biomass[,c('x','y')], mean=hh, sd=hh.sd)
#soft2 = cbind(biomass[,c('x','y')], mean=ss, sd=ss.sd)
hard2 = cbind(biomass.wi[,c('x','y')], obs=biomass.wi$hardwood, mean=hh, sd=hh.sd)
soft2 = cbind(biomass.wi[,c('x','y')], obs=biomass.wi$softwood, mean=ss, sd=ss.sd)
tot2 = cbind(biomass.wi[,c('x','y')], obs=biomass.wi$tot, mean=tt, sd=tt.sd)


rr = c(min(min(tot2$mean), min(hard2$mean), min(soft2$mean)), max(max(tot2$mean), max(hard2$mean), max(soft2$mean)))
rr.sd = c(min(min(tot2$sd), min(hard2$sd), min(soft2$sd)), max(max(tot2$sd), max(hard2$sd), max(soft2$sd)))

#Maps of raw biomass:
plots = list(
    total=list(data=tot2, caption='total'),
    hardwood=list(data=hard2, caption='hardwood'),
    softwood=list(data=soft2, caption='softwood')
)

for (item in plots) {
    p0 = ggplot(item[['data']]) +
        aes(x,y) +
        aes_string(color='obs') +
        scale_color_gradient("biomass", trans='log', low='white', high='blue',
            limits=rr,
            #breaks=c(2, 150, 8100, 442400),
            #labels=c("0.002", "0.150", "8.100", "442.4")) +
            breaks=c(1, 20, 400, 8000)) +
            
        geom_point(shape=15, solid=TRUE) #+
        #ggtitle(paste("Mean ", item[['caption']], " biomass", sep=''))

    pdf(paste("~/misc/paleon/figures/", item[['caption']], "-WI-biomass-raw.pdf", sep=''), width=5, height=4)
    print(p0)
    dev.off()

    p1 = ggplot(item[['data']]) +
        aes(x,y) +
        aes_string(color='mean') +
        scale_color_gradient("mean", trans='log', low='white', high='blue',
            limits=rr,
            #breaks=c(2, 150, 8100, 442400),
            #labels=c("0.002", "0.150", "8.100", "442.4")) +
            breaks=c(1, 20, 400, 8000)) +
            
        geom_point(shape=15, solid=TRUE) #+
        #ggtitle(paste("Mean ", item[['caption']], " biomass", sep=''))

    pdf(paste("~/misc/paleon/figures/", item[['caption']], "-WI-biomass-mean.pdf", sep=''), width=5, height=4)
    print(p1)
    dev.off()

    p2 = ggplot(item[['data']]) +
        aes(x,y) +
        aes_string(color='sd') +
        scale_color_gradient("st. dev.", trans='log', low='white', high='blue',
            limits=rr.sd,
            breaks=c(2, 50, 1000)) +
            
        geom_point(shape=15, solid=TRUE) #+
        #ggtitle(paste("Standard deviation of ", item[['caption']], " biomass", sep=''))

    #pdf(paste("/NO_BACKUP/brooks/paleon/paleon/figures/biomass/", item[['caption']], "-sd-biomass.pdf", sep=''), width=10, height=6)
    pdf(paste("~/misc/paleon/figures/", item[['caption']], "-WI-biomass-sd.pdf", sep=''), width=5, height=4)
    print(p2)
    dev.off()
}






#Plot the mean and standard deviation of draws for hardwood, softwood biomass:
require(fBasics)
hard = read.csv("/NO_BACKUP/brooks/paleon/paleon/output/logbiomass-NA-hardwood.csv", header=FALSE)
soft = read.csv("/NO_BACKUP/brooks/paleon/paleon/output/logbiomass-NA-softwood.csv", header=FALSE)

hh = sapply(1:length(hard), function(k) {mean(exp(hard[[k]]))})
ss = sapply(1:length(soft), function(k) {mean(exp(soft[[k]]))})

hh.sd = sapply(1:length(hard), function(k) {sd(exp(hard[[k]]))})
ss.sd = sapply(1:length(soft), function(k) {sd(exp(soft[[k]]))})

rm(hard)
rm(soft)

#Trim out the cells where the biomass was off the scale
indx.h = which(!is.na(hh) & !is.na(hh.sd) & hh<Inf & hh>-Inf & hh.sd<Inf)
hh = hh[indx.h]
hh.sd = hh.sd[indx.h]

indx.s = which(!is.na(ss) & !is.na(ss.sd) & ss<Inf & ss>-Inf & ss.sd<Inf)
ss = ss[indx.s]
ss.sd = ss.sd[indx.s]

hard2 = cbind(biomass[indx.h,c('x','y')], mean=hh, sd=hh.sd)
soft2 = cbind(biomass[indx.s,c('x','y')], mean=ss, sd=ss.sd)

rr = c(min(min(hard2$mean), min(soft2$mean)), max(max(hard2$mean), max(soft2$mean)))
rr.sd = c(min(min(hard2$sd), min(soft2$sd)), max(max(hard2$sd), max(soft2$sd)))

#Maps of raw biomass:
plots = list(hardwood=list(data=hard2, caption='hardwood', meanbreaks=c(50, 1000, 20000, 400000), sdbreaks=c(50, 1000, 20000)),
    softwood=list(data=soft2, caption='softwood', meanbreaks=c(1e23, 1e56, 1e89, 1e122), sdbreaks=c(1e24, 1e57, 1e90, 1e123))
)

for (item in plots) {
    p = ggplot(item[['data']]) +
        aes(x,y) +
        aes_string(color='mean') +
        scale_color_gradient("biomass", trans='log', low='white', high='blue',
            breaks=item[['meanbreaks']]) +
            
        geom_point(shape=15, solid=TRUE) #+
        #ggtitle(paste("Mean ", item[['caption']], " biomass", sep=''))

    pdf(paste("/NO_BACKUP/brooks/paleon/paleon/figures/biomass/", item[['caption']], "-mean-biomass.pdf", sep=''), width=10, height=6)
    #pdf(paste("~/misc/paleon/figures/", item[['caption']], "-mean-wi-biomass.pdf", sep=''), width=5, height=4)
    print(p)
    dev.off()


    p = ggplot(item[['data']]) +
        aes(x,y) +
        aes_string(color='sd') +
        scale_color_gradient("st. dev.", trans='log', low='white', high='blue',
            breaks=item[['sdbreaks']]) +
            
        geom_point(shape=15, solid=TRUE) #+
        #ggtitle(paste("Standard deviation of ", item[['caption']], " biomass", sep=''))

    pdf(paste("/NO_BACKUP/brooks/paleon/paleon/figures/biomass/", item[['caption']], "-sd-biomass.pdf", sep=''), width=10, height=6)
    #pdf(paste("~/misc/paleon/figures/", item[['caption']], "-sd-wi-biomass.pdf", sep=''), width=5, height=4)
    print(p)
    dev.off()
}
