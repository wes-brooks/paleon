#load necessary R libraries
library(mgcv)
library(tweedie)
library(Matrix)
library(plotrix)



#Import the data and the heatmap code.
source("~/git/brooks/code/matplot.r")
source("code/biomass-import.r")
taxa = c('Cherries', 'Willow', 'Walnuts', 'Hickory', 'Beech', 'Fir', 'Spruce', 'Ironwoods', 'Cedar', 'Hemlock', 'Basswood', 'Ashes', 'Elms', 'Poplar', 'Pine', 'Tamarack', 'Birches', 'Maple', 'Oaks')
#args <- commandArgs(trailingOnly = TRUE)
#indx = as.numeric(args[1])

for (sp in c(taxa[1])) {

    #################################################################
    #Data - neighborhood structure (first-order)
    #################################################################
    #Get the first-order neighborhood of each observation
    n = nrow(biomass.wi)
    D = as.matrix(dist(biomass.wi[,c('x','y')], diag=TRUE), n, n)
    diag(D) = Inf
    neighbors = list()
    graphfile = 'data/neighborhood.wi.txt'
    cat(paste(n, '\n', sep=''), file=graphfile)
    for (i in 1:n) {
        neighbors[[i]] = which(D[i,]==min(D[i,]))
        cat(paste(i, length(neighbors[[i]]), paste(as.vector(neighbors[[i]]), collapse=' '), '\n', sep=' '), file=graphfile, append=TRUE)
    }
    neighborhood = Matrix(0, n, n)
    for (i in 1:n) {
        neighborhood[i,neighbors[[i]]] = 1
    }


    #################################################################
    #Modeling - Tweedie
    #################################################################
    #Make the one-stage model
    #Set up the data, including mean composition within the first-order neighborhood:
    modeldata = list(loc=seq(nrow(biomass.wi)), indicator=ifelse(biomass.wi[,sp]>0, 1, 0), biomass=biomass.wi[,sp], logbiomass=log(biomass.wi[,sp]), biomass.nz=ifelse(biomass.wi[,sp]>0,biomass.wi[,sp],NA), composition=composition.wi[,sp], x=composition.wi[,'x'], y=composition.wi[,'y'])
    modeldata = as.data.frame(modeldata)
    neighborhood.composition = vector()
    for (i in 1:n) {
        neighborhood.composition = c(neighborhood.composition, sum(modeldata$composition * neighborhood[i,]) / sum(neighborhood[i,]))
    }
    modeldata$neighborhood.composition = neighborhood.composition
    modeldata$log.neighborhood.composition = ifelse(neighborhood.composition==0,NA,log(neighborhood.composition))




    #################################################################
    #Plotting - Heatmap(observed)
    #################################################################
    #Put the observed biomass into matrices that we can plot as heatmaps
    #Create the matrix for the observed stem density (filled by default with NAs):
    loc = with(modeldata, list(lat=unique(y), long=unique(x)))
    biomat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
    rownames(biomat) <- sort(unique(modeldata$y), decreasing=F)
    colnames(biomat) <- sort(unique(modeldata$x), decreasing=F)

    #Put the observed biomass into its lat-long matrices
    for(row in 1:length(modeldata[['x']])) {
        biomat[as.character(modeldata$y[row]), as.character(modeldata$x[row])] = modeldata$logbiomass[row]
    }
    biomat = ifelse(biomat==-Inf,NA,biomat)

    #Plot the model and the data:
    xx = range(biomat, na.rm=TRUE)
    xx[1] = floor(10*xx[1])/10
    xx[2] = ceiling(10*xx[2])/10
    dev.new()
    par(bty='n')
    gwr.matplot(biomat, c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75", xrange=xx)
    title(paste(sp, " (observed): biomass (log scale)", sep=""))
}