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

    n.zero = length(which(modeldata$biomass==0))

    #Plot the model and the data:
    pdf(paste('figures/raw/', sp, "-observed-heatmap.pdf", sep=""), width=12, height=6)
    layout(t(matrix(1:2)))
    par(bty='n')
    gwr.matplot(biomat, c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75")
    title(paste('Observed ', sp, " biomass (log scale)", sep=""))
    hist(modeldata$biomass, breaks=50, main=paste(sp, ' biomass', sep=''), xlab=paste("cell biomass (zeros: ", n.zero, ")", sep=''))
    dev.off()


    #################################################################
    #Plotting - vs Composition(observed)
    #################################################################
    #Put the observed biomass into matrices that we can plot as heatmaps
    posdata = modeldata[modeldata$biomass>0,]
    g1 = gam(log(biomass)~s(neighborhood.composition, k=50), data=posdata, family=gaussian, gamma=1.4)

    #Plot the model and the data:
    xx = range(posdata$composition) * 1.1
    yy = range(log(posdata$biomass)) * c(0.8, 1.1)
    fits = fitted(g1)
    pdf(paste('figures/raw/', sp, "-biomass-v-composition.pdf", sep=""), width=5, height=6)
    par(bty='n')
    plot(posdata$neighborhood.composition, posdata$logbiomass, pch=16, cex=0.35, xlim=xx, ylim=yy, bty='n', xlab='composition', ylab='log(biomass)')
    title(paste(sp, " biomass vs. composition", sep=""))
    par(new=TRUE)
    plot(sort(posdata$neighborhood.composition), fits[order(posdata$neighborhood)], type='l', lwd=2, xlim=xx, ylim=yy, ann=FALSE, xaxt='n', yaxt='n', col='purple')
    dev.off()


    #################################################################
    #Modeling - Aspatial/GLM/Tweedie
    #################################################################
    bm.opt = function(theta, data, k=150) {
        result = list()

        data = as.data.frame(data)
        model = glm(biomass ~ neighborhood.composition * stemdensity, control=glm.control(maxit=100), data=data, family=tweedie(var.power=theta, link.power=0))
        cat(paste("theta: ", theta, '\n', sep=''))

        #Get the scale (a) and the location (b)
        abs(resid(model, type='deviance')) -> scale
        predict(model, type='link') -> loc
        m = lm(scale~loc)

        return(coef(m)[2]**2)
    }
    tuning = optimize(bm.opt, interval=c(1,2), data=modeldata, k=15, tol=0.01)


    #################################################################
    #Modeling - Aspatial/GAM/Tweedie
    #################################################################


    #Find the optimal Tweedie parameter
    bm.opt = function(theta, data, k=150) {
        result = list()

        data = as.data.frame(data)
        model = gam(biomass ~ s(neighborhood.composition, stemdensity, k=k), control=glm.control(maxit=100), data=data, family=Tweedie(p=theta, link='log'))
        cat(paste("theta: ", theta, '\n', sep=''))

        #Get the scale (a) and the location (b)
        abs(resid(model, type='deviance')) -> scale
        predict(model, type='link') -> loc
        m = lm(scale~loc)

        return(coef(m)[2]**2)
    }
    k=100
    tuning = optimize(bm.opt, interval=c(1,2), data=modeldata, k=k, tol=0.01)
    gam1 = gam(biomass ~ s(neighborhood.composition, stemdensity, k=k), data=modeldata, family=Tweedie(p=tuning$minimum, link='log'))
    gam1.resid = resid(gam1)

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

    #Plot the inputs and the fitted smooth
    pdf(paste('figures/aspatial-tweedie/', sp, "-plots.pdf", sep=""), width=12, height=4)
    layout(t(matrix(1:3)))
    symbols(modeldata$neighborhood.composition, modeldata$stemdensity, circles=sqrt(modeldata$biomass), fg=NA, bg=rgb(255,0,0,60,maxColorValue=255), inches=0.1, bty='n', xlab=paste(sp, " composition fraction (smoothed)", sep=''), ylab="Aggregate stem density")
    plot(gam1, pers=TRUE, xlab='composition', ylab='stem density', main='biomass smooth')
    par(bty='n')
    gwr.matplot(biomat, c(0.6,1), c(0,0), c(0.8,0.2), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75", xcrit=0, col.crit=c(1,1,1))
    title(paste(sp, " biomass residuals (log scale)", sep=""))
    dev.off()


    #################################################################
    #Modeling - Aspatial/GAM/Delta
    #################################################################
    #Plot the inputs
    pdf(paste('figures/aspatial-delta/', sp, "-composition-v-presence.pdf", sep=""), width=12, height=6)
    plot(modeldata$neighborhood.composition, jitter(ifelse(modeldata$biomass>0,1,0), factor=0.2), pch=16, cex=0.3, bty='n', xlab=paste(sp, " composition fraction (smoothed)", sep=''), ylab=paste(sp, " presence", sep=""))
    dev.off()

    #Fit the GAM model for the gamma density
    k=25
    indx = which(modeldata$biomass>0)
    posdata = as.data.frame(modeldata[indx,])    
    gam1 = gam(biomass ~ s(neighborhood.composition, stemdensity, k=k), control=glm.control(maxit=100), data=posdata, family=Gamma(link='log'))
    gam1.resid = resid(gam1)

    #Put the fitted biomass into matrices that we can plot as heatmaps
    #Create the matrix for the observed stem density (filled by default with NAs):
    loc = with(modeldata, list(lat=unique(y), long=unique(x)))
    biomat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
    rownames(biomat) <- sort(unique(modeldata$y), decreasing=F)
    colnames(biomat) <- sort(unique(modeldata$x), decreasing=F)

    #Put the observed biomass into its lat-long matrices
    for(row in 1:length(posdata[['x']])) {
        biomat[as.character(posdata$y[row]), as.character(posdata$x[row])] = gam1$resid[row]
    }
    biomat = ifelse(biomat==-Inf,NA,biomat)

    #Plot the model and the data:
    pdf(paste('figures/aspatial-delta/', sp, "-gamma-residuals-heatmap.pdf", sep=""), width=6, height=6)
    par(bty='n')
    gwr.matplot(biomat, c(0.6,1), c(0,0), c(0.8,0.2), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75", xcrit=0, col.crit=c(1,1,1))
    title(paste(sp, " biomass residuals (log scale)", sep=""))
    dev.off()


    ##################################################################
    ##Modeling - Aspatial/GAM/Delta (COZIGAM)
    ##################################################################
    ##Fit the GAM model for the gamma density
    #k=25
    #gam1 = zigam(biomass ~ s(neighborhood.composition, k=k), data=modeldata, constraint='proportional', family='Gamma(link=log)')
    #gam1.resid = resid(gam1)
    #
    ##Put the fitted biomass into matrices that we can plot as heatmaps
    ##Create the matrix for the observed stem density (filled by default with NAs):
    #loc = with(modeldata, list(lat=unique(y), long=unique(x)))
    #biomat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
    #rownames(biomat) <- sort(unique(modeldata$y), decreasing=F)
    #colnames(biomat) <- sort(unique(modeldata$x), decreasing=F)
    #
    ##Put the observed biomass into its lat-long matrices
    #for(row in 1:length(modeldata[['x']])) {
    #    biomat[as.character(modeldata$y[row]), as.character(modeldata$x[row])] = gam1$resid[row]
    #}
    #biomat = ifelse(biomat==-Inf,NA,biomat)
    #
    ##Plot the model and the data:
    #pdf(paste('figures/aspatial-cozigam/', sp, "-residuals-heatmap.pdf", sep=""), width=6, height=6)
    #par(bty='n')
    #gwr.matplot(biomat, c(0.6,1), c(0,0), c(0.8,0.2), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75", xcrit=0, col.crit=c(1,1,1))
    #title(paste(sp, " biomass residuals (log scale)", sep=""))
    #dev.off()


    ##################################################################
    ##Modeling - Spatial/GAM/Delta (COZIGAM)
    ##################################################################
    ##Fit the GAM model for the gamma density
    #k1 = 75
    #k2 = 25
    #gam1 = cozigam(biomass ~ s(x,y,k=k1) + s(neighborhood.composition, stemdensity, k=k2), data=modeldata, constraint='proportional', family='Gamma(link=log)')
    #gam1.resid = resid(gam1)
    #
    ##Put the fitted biomass into matrices that we can plot as heatmaps
    ##Create the matrix for the observed stem density (filled by default with NAs):
    #loc = with(modeldata, list(lat=unique(y), long=unique(x)))
    #biomat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
    #rownames(biomat) <- sort(unique(modeldata$y), decreasing=F)
    #colnames(biomat) <- sort(unique(modeldata$x), decreasing=F)
    #
    ##Put the observed biomass into its lat-long matrices
    #for(row in 1:length(modeldata[['x']])) {
    #    biomat[as.character(modeldata$y[row]), as.character(modeldata$x[row])] = gam1$resid[row]
    #}
    #biomat = ifelse(biomat==-Inf,NA,biomat)
    #
    ##Plot the model and the data:
    #pdf(paste('figures/aspatial-cozigam/', sp, "-residuals-heatmap.pdf", sep=""), width=6, height=6)
    #par(bty='n')
    #gwr.matplot(biomat, c(0.6,1), c(0,0), c(0.8,0.2), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75", xcrit=0, col.crit=c(1,1,1))
    #title(paste(sp, " biomass residuals (log scale)", sep=""))
    #dev.off()



    #################################################################
    #Modeling - Spatial/GAM/Tweedie
    #################################################################
    #Find the optimal Tweedie parameter
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
    tuning = optimize(bm.opt, interval=c(1,2), data=modeldata, k1=k1, k2=k2, tol=0.01)
    gam1 = gam(biomass ~ s(x,y,k=k1) + s(neighborhood.composition, stemdensity, k=k2), control=glm.control(maxit=100), data=modeldata, family=Tweedie(p=tuning$minimum, link='log'))
    gam1.resid = resid(gam1)

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

    ##################################################################
    ##Modeling - INLA
    ##################################################################
    ##Make the one-stage model
    #
    mi1 = inla(log(biomass.nz) ~ f(loc, model='besag', graph=graphfile), data=modeldata, verbose=TRUE, family='gaussian')
    mi2 = inla(indicator ~ f(loc, model='besag', graph=graphfile), data=modeldata, verbose=TRUE, family='binomial')
    fits1 = mi1$summary.random$loc$mean
    fits2 = mi2$summary.random$loc$mean

    #Put the fitted biomass into matrices that we can plot as heatmaps
    #Create the matrix for the observed stem density (filled by default with NAs):
    loc = with(modeldata, list(lat=unique(y), long=unique(x)))
    biomat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
    rownames(biomat) <- sort(unique(modeldata$y), decreasing=F)
    colnames(biomat) <- sort(unique(modeldata$x), decreasing=F)

    #Put the observed biomass into its lat-long matrices
    for(row in 1:length(modeldata[['x']])) {
        biomat[as.character(modeldata$y[row]), as.character(modeldata$x[row])] = fits1[row]
    }
    biomat = ifelse(biomat==-Inf,NA,biomat)

    #Create the matrix for the fitted probability of being nonzero
    loc = with(modeldata, list(lat=unique(y), long=unique(x)))
    biomat2 = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
    rownames(biomat2) <- sort(unique(modeldata$y), decreasing=F)
    colnames(biomat2) <- sort(unique(modeldata$x), decreasing=F)

    #Put the observed biomass into its lat-long matrices
    for(row in 1:length(modeldata[['x']])) {
        biomat2[as.character(modeldata$y[row]), as.character(modeldata$x[row])] = fits2[row]
    }
    biomat2 = ifelse(biomat2==-Inf,NA,biomat2)

    #Plot residuals
    resids = log(modeldata$biomass.nz) - mi1$summary.random$loc$mean - mi1$summary.fixed[1]
    #Put the fitted biomass into matrices that we can plot as heatmaps
    #Create the matrix for the observed stem density (filled by default with NAs):
    loc = with(modeldata, list(lat=unique(y), long=unique(x)))
    biomat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
    rownames(biomat) <- sort(unique(modeldata$y), decreasing=F)
    colnames(biomat) <- sort(unique(modeldata$x), decreasing=F)

    #Put the observed biomass into its lat-long matrices
    for(row in 1:length(modeldata[['x']])) {
        biomat[as.character(modeldata$y[row]), as.character(modeldata$x[row])] = resids[row]
    }
    biomat = ifelse(biomat==-Inf,NA,biomat)

    #Plot the model and the data:
    pdf(paste('figures/inla-delta/', sp, "-residuals-heatmap.pdf", sep=""), width=6, height=6)
    par(bty='n')
    gwr.matplot(biomat, c(0.6,1), c(0,0), c(0.8,0.2), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75", xcrit=0, col.crit=c(1,1,1))
    title(paste(sp, " biomass residual (log scale)", sep=""))
    dev.off()
}