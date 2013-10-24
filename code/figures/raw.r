require(ggplot2)
require(RCurl)
require(devtools)
#If the 'brooks' package isnt loaded then import it from github:
if (!'package:brooks' %in% search()) {
    install_github('wesesque/brooks')
    require(brooks)
}

#import the data
source_url("https://raw.github.com/wesesque/paleon/master/code/biomass/biomass-import.r")


#Map of biomass:
p = ggplot(biomass) +
    aes(x,y) +
    aes(color=tot) +
    scale_color_gradient("biomass (1000s)", trans='log', low='white', high='blue',
        breaks=c(2, 150, 8100, 442400),
        labels=c("0.002", "0.150", "8.100", "442.4")) +
    geom_point() +
    ggtitle("Total tree biomass")

pdf("~/git/paleon/figures/prelim-talk/raw-biomass.pdf", width=10, height=6)
p
dev.off()


#Histogram of biomass (marginal):
pp = ggplot(biomass) + aes(x=tot) + geom_histogram() + geom_histogram(aes(y = ..density..)) + geom_density()