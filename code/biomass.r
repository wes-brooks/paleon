setwd("~/git/paleon")

#load necessary R libraries
library(mgcv)
#library(reshape)
#library(plotrix)
#library(gwselect)
#registerCores(n=3)

#Set the working directory and import the code we'll use to make heatmaps.
source("~/git/brooks/code/matplot.r")
models = list()

#Import the stem density, count, and standard deviation data
composition <- read.csv("data/glo.forest.composition_v1_3alb.csv", header=T)
biomass <- read.csv("data/biomass.tablev1_3alb.csv", header=T)

#Process composition data
p = ncol(composition)
composition$tot = apply(composition[,4:p], 1, sum)
composition[,4:p] = composition[,4:p] / composition$tot
composition = composition[which(!is.na(composition$tot)),]

#Process biomass data
p2 = ncol(biomass)
biomass$tot = apply(biomass[,4:p2], 1, sum)
biomass = biomass[which(!is.na(biomass$tot)),]
composition = cbind(composition, biomass[,4:p2])


sp = "Oaks"
models = list()
gm = list()
indx = which(biomass[,sp]>0)
modeldata = list(indicator=ifelse(biomass[biomass=biomass[indx,sp], logbiomass=log(biomass[indx,sp]), composition=composition[indx,sp], x=composition[indx,'x'], y=composition[indx,'y'])
modeldata = list(indicator=ifelse(biomabiomass=biomass[indx,sp], logbiomass=log(biomass[indx,sp]), composition=composition[indx,sp], x=composition[indx,'x'], y=composition[indx,'y'])
paleon.dist = dist(cbind(modeldata$x, modeldata$y))


sp = "Oaks"
models = list()
gm = list()
indx = which(biomass[,sp]>0)
modeldata = list(indicator=ifelse(biomabiomass=biomass[indx,sp], logbiomass=log(biomass[indx,sp]), composition=composition[indx,sp], x=composition[indx,'x'], y=composition[indx,'y'])
modeldata = list(indicator=ifelse(biomabiomass=biomass[indx,sp], logbiomass=log(biomass[indx,sp]), composition=composition[indx,sp], x=composition[indx,'x'], y=composition[indx,'y'])
paleon.dist = dist(cbind(modeldata$x, modeldata$y))
pd2 = as.matrix(paleon.dist)


modeldata = data.frame(biomass=biomass[,sp], logbiomass=log(biomass[,sp]), composition=composition[,sp], x=composition$x, y=composition$y)
modeldata = modeldata[which(modeldata$biomass>0),]
models[[sp]] = gam(logbiomass~s(composition, k=100) + s(x, y, k=200), data=modeldata, family="gaussian")


gm[[sp]] = gam(biomass~s(composition, k=100) + s(x,y, k=100), data=modeldata, family="Gamma")
gm[[sp]] = gam(biomass~s(composition, k=50), data=modeldata, family="gaussian")

gam(modeldata$biomass~s(modeldata$composition, k=50), family="Gamma")



#Set the directory to store plots
plot_dir = "~/git/paleon/figures/"
#model = "gamma"
#genus = "tot"
#xx = range(stemdensity[,genus], na.rm=TRUE)


pdf(paste(plot_dir, sp, "-biomass-spline.pdf", sep=""))
xx = c(0,1)
yy = c(-5,12)
plot(modeldata$composition, log(modeldata$biomass), cex=0.2, pch=20, xlim=xx, ylim=yy, col='purple', bty='n', xlab=paste(sp, " proportion", sep=""), ylab=paste(sp, " biomass", sep=""))
par(new=T)
plot(modeldata$composition, fitted(models[[sp]]), cex=0.2, pch=20, xlim=xx, ylim=yy, col='black', bty='n', ann=F, xaxt='n', yaxt='n')
title(paste(sp, " biomass spline", sep=""))
dev.off()


plot(modeldata$composition, fitted(gm[[sp]]), cex=0.2, pch=20, xlim=xx, ylim=yy, col='black', bty='n', ann=F, xaxt='n', yaxt='n')


#Put the observed stem density and the weight into matrices that we can plot as heatmaps
#Create the matrix for the observed stem density (filled by default with NAs):
loc = with(modeldata, list(lat=unique(y), long=unique(x)))
bmmat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
rownames(bmmat) <- sort(unique(modeldata$y), decreasing=F)
colnames(bmmat) <- sort(unique(modeldata$x), decreasing=F)

#Put the observed stem densities and weights into their lat-long matrices
for(row in 1:dim(modeldata)[1]) {
    bmmat[as.character(modeldata$y[row]), as.character(modeldata$x[row])] = ifelse((!is.nan(modeldata$biomass[row]) & modeldata$biomass[row]>0), modeldata$biomass[row], NA)
}


#Put the model's fitted values into a matrix that we can plot as a heatmap (indexed by Lat, Long):
#First, create a data frame with Lat, Long, and fitted stem density
rows = which(!is.na(modeldata$biomass))
fits = cbind(modeldata[rows, c("y", "x")], fitted=fitted(models[[sp]]))

#Now create a matrix of the appropriate dimensions (filled with NAs by default)
fitted = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
rownames(fitted) <- sort(unique(modeldata$y), decreasing=F)
colnames(fitted) <- sort(unique(modeldata$x), decreasing=F)

#Put the fittes stem densities into the lat-long matrix
for(row in 1:dim(modeldata)[1])
    fitted[as.character(fits$y[row]), as.character(fits$x[row])] = fits$fitted[row]                  


resid = log(bmmat) - fitted


#Make several plots of the model and the data:
pdf(paste(plot_dir, sp, "-biomass-residual.pdf", sep=""))
par(bty='n')
gwr.matplot(resid, c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75")
title(paste(sp, " biomass residual (log scale)", sep=""))
dev.off()


pdf(paste(plot_dir, sp, "-biomass-observed.pdf", sep=""))
par(bty='n')
gwr.matplot(log(bmmat), c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75")
title(paste(sp, " observed biomass (log scale)", sep=""))
dev.off()


pdf(paste(plot_dir, sp, "-biomass-fitted.pdf", sep=""))
par(bty='n')
gwr.matplot(fitted, c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75")
title(paste(sp, " fitted biomass (log scale)", sep=""))
dev.off()
