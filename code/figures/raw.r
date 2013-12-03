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


#Maps of raw biomass:
plots = list(tot=list(part='tot', caption='Total tree biomass', file='tot-raw-biomass'),
    hardwood=list(part='hardwood', caption='Total hardwood tree biomass', file='hardwood-raw-biomass'),
    softwood=list(part='softwood', caption='Total softwood tree biomass', file='softwood-raw-biomass')
)

data = biomass

for (item in plots) {
    p = ggplot(data) +
        aes(x,y) +
        aes_string(color=item[['part']]) +
        scale_color_gradient("biomass", trans='log', low='white', high='blue',
            limits=range(data$tot[data$tot>0]),
            breaks=c(2, 150, 8000, 400000)) +
        geom_point(shape=15, solid=TRUE) +
        ggtitle(item[['caption']])

    pdf(paste("~/git/paleon/figures/prelim-talk/", item[['file']], ".pdf", sep=''), width=10, height=6)
    print(p)
    dev.off()
}


#Histogram of biomass (marginal):
pp = ggplot(biomass) + aes(x=tot) + geom_histogram() + geom_histogram(aes(y = ..density..)) + geom_density()