#load necessary R libraries
library(mgcv)
library(tweedie)
library(Matrix)
library(plotrix)


#Import the data and the heatmap code.
source("~/git/brooks/code/matplot.r")
source("code/biomass-import.r")
#taxa = c('Cherries', 'Willow', 'Walnuts', 'Hickory', 'Beech', 'Fir', 'Spruce', 'Ironwoods', 'Cedar', 'Hemlock', 'Basswood', 'Ashes', 'Elms', 'Poplar', 'Pine', 'Tamarack', 'Birches', 'Maple', 'Oaks')
#args <- commandArgs(trailingOnly = TRUE)
#indx = as.numeric(args[1])
taxa = c('Oaks', 'Maple')

for (sp in taxa) {
    ################################################################
    #Data - neighborhood structure (first-order)
    ################################################################
    #Get the first-order neighborhood of each observation for INLA
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

    #Get the first-order neighborhood of each observation
    n = nrow(biomass.wi)
    D = as.matrix(dist(biomass.wi[,c('x','y')], diag=TRUE), n, n)
    #diag(D) = Inf
    neighbors = list()
    for (i in 1:n) {
        neighbors[[i]] = which(D[i,]==min(D[i,]))
    }
    neighborhood = Matrix(0, n, n)
    for (i in 1:n) {
        neighborhood[i,neighbors[[i]]] = 1
    }

    #Set up the data, including mean composition within the first-order neighborhood:
    modeldata = list(loc=seq(nrow(biomass.wi)), indicator=ifelse(biomass.wi[,sp]>0, 1, 0), biomass=biomass.wi[,sp], logbiomass=log(biomass.wi[,sp]), biomass.nz=ifelse(biomass.wi[,sp]>0,biomass.wi[,sp],NA), composition=composition.wi[,sp], stemdensity=stemdens.wi$stem.density, x=composition.wi[,'x'], y=composition.wi[,'y'])
    modeldata = as.data.frame(modeldata)
    neighborhood.composition = vector()
    for (i in 1:n) {
        neighborhood.composition = c(neighborhood.composition, sum(modeldata$composition * neighborhood[i,]) / sum(neighborhood[i,]))
    }
    modeldata$neighborhood.composition = neighborhood.composition
    modeldata$log.neighborhood.composition = ifelse(neighborhood.composition==0,NA,log(neighborhood.composition))
    modeldata$count = stemdens.wi$stem.density * neighborhood.composition





    #################################################################
    #Modeling - Spatial/GAM/Tweedie
    #################################################################
    #Find the optimal Tweedie parameter
    theta = 1.43
    bm.opt = function(theta, data, k1=150, k2=50) {
        result = list()

        data = as.data.frame(data)
        model = gam(biomass ~ s(x,y, k=k1) + s(neighborhood.composition, stemdensity, k=k2), control=glm.control(maxit=100), data=data, family=Tweedie(p=theta, link='log'))
        cat(paste("theta: ", theta, '\n', sep=''))

        #Get the scale (a) and the location (b)
        abs(resid(model, type='deviance')) -> scale
        predict(model, type='link') -> loc
        m = lm(scale~loc)

        return(coef(m)[2]**2)
    }
    k1 = 150
    k2 = 50
    #tuning = optimize(bm.opt, interval=c(1,2), data=modeldata, k1=k1, k2=k2, tol=0.01)
    gam1 = gam(biomass ~ s(x,y,k=k1) + s(neighborhood.composition, stemdensity, k=k2), control=glm.control(maxit=100), data=modeldata, family=Tweedie(p=theta, link='log'))
    gam1.resid = resid(gam1)
    gam1.fits = exp(-summary(gam1)$scale^-1 * exp((fitted(gam1))^(2-theta))/(2-theta))

    #Put the fitted biomass into matrices that we can plot as heatmaps
    #Create the matrix for the observed stem density (filled by default with NAs):
    loc = with(modeldata, list(lat=unique(y), long=unique(x)))
    biomat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
    rownames(biomat) <- sort(unique(modeldata$y), decreasing=F)
    colnames(biomat) <- sort(unique(modeldata$x), decreasing=F)

    #Put the observed biomass into its lat-long matrices
    for(row in 1:length(modeldata[['x']])) {
        biomat[as.character(modeldata$y[row]), as.character(modeldata$x[row])] = gam1$resid[row]
    }
    biomat = ifelse(biomat==-Inf,NA,biomat)

    #Plot the model and the data:
    pdf(paste('figures/spline-tweedie/', sp, "-residuals-heatmap.pdf", sep=""), width=12, height=4)
    layout(t(matrix(1:3)))
    plot(gam1, pers=TRUE, select=1, main='spatial smooth')
    plot(gam1, pers=TRUE, select=2,, xlab='composition', ylab='stem density', main='biomass smooth')
    par(bty='n')
    gwr.matplot(biomat, c(0.6,1), c(0,0), c(0.8,0.2), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75", xcrit=0, col.crit=c(1,1,1))
    title(paste(sp, " biomass residuals (log scale)", sep=""))
    dev.off()


    #################################################################
    #Modeling - Spatial/GAM/Delta
    #################################################################
    #Find the optimal Tweedie parameter
    k1 = 150
    k2 = 50
    gam2 = gam(indicator ~ s(x,y,k=k1) + s(neighborhood.composition, stemdensity, k=k2), control=glm.control(maxit=100), data=modeldata, family='binomial')
    gam2.fits = exp(fitted(gam2))/(1+exp(fitted(gam2)))

    #Put the fitted biomass into matrices that we can plot as heatmaps
    #Create the matrix for the observed stem density (filled by default with NAs):
    loc = with(modeldata, list(lat=unique(y), long=unique(x)))
    biomat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
    rownames(biomat) <- sort(unique(modeldata$y), decreasing=F)
    colnames(biomat) <- sort(unique(modeldata$x), decreasing=F)

    #Put the observed biomass into its lat-long matrices
    for(row in 1:length(modeldata[['x']])) {
        biomat[as.character(modeldata$y[row]), as.character(modeldata$x[row])] = gam1.fits[row]
    }
    biomat = ifelse(biomat==-Inf,NA,biomat)

    #Put the fitted biomass into matrices that we can plot as heatmaps
    #Create the matrix for the observed stem density (filled by default with NAs):
    loc = with(modeldata, list(lat=unique(y), long=unique(x)))
    biomat2 = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
    rownames(biomat2) <- sort(unique(modeldata$y), decreasing=F)
    colnames(biomat2) <- sort(unique(modeldata$x), decreasing=F)

    #Put the observed biomass into its lat-long matrices
    for(row in 1:length(modeldata[['x']])) {
        biomat2[as.character(modeldata$y[row]), as.character(modeldata$x[row])] = gam2.fits[row]
    }
    biomat2 = ifelse(biomat2==-Inf,NA,biomat2)



    #Plot the model and the data:
    pdf(paste('figures/spline-tweedie/', sp, "-comparison   .pdf", sep=""), width=12, height=6)
    layout(t(matrix(1:2)))
    par(bty='n')
    gwr.matplot(biomat, c(1,1), c(1,0), c(1,0.2), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75", xrange=c(0,1))
    title(paste(sp, " P(Y=0) (Tweedie)", sep=""))
    par(bty='n')
    gwr.matplot(biomat2, c(1,1), c(1,0), c(1,0.2), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75", xrange=c(0,1))
    title(paste(sp, " P(Y=0) (Delta)", sep=""))
    dev.off()
}