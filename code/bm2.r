setwd("~/git/paleon")

#load necessary R libraries
library(mgcv)
library(tweedie)
library(Matrix)

#Import the data and the heatmap code.
source("~/git/brooks/code/matplot.r")
source("code/biomass-import.r")
models = list()
sp = "Oaks"


#######################################
#Get the first-order neighborhood of each observation
n = nrow(biomass.wi)
D = as.matrix(dist(biomass.wi[,c('x','y')], diag=TRUE), n, n)
diag(D) = Inf
neighbors = list()
for (i in 1:n) {
    neighbors[[i]] = which(D[i,]==min(D[i,]))
}
neighborhood = Matrix(0, n, n)
for (i in 1:n) {
    neighborhood[i,neighbors[[i]]] = 1
}




#######################################
#Make the one-stage model

#Set up the data, including mean composition within the first-order neighborhood:
modeldata = list(indicator=ifelse(biomass.wi[,sp]>0, 1, 0), biomass=biomass.wi[,sp], logbiomass=log(biomass.wi[,sp]), composition=composition.wi[,sp], x=composition.wi[,'x'], y=composition.wi[,'y'])
neighborhood.composition = vector()
for (i in 1:n) {
    neighborhood.composition = c(neighborhood.composition, sum(modeldata$composition * neighborhood[i,]) / sum(neighborhood[i,]))
}
modeldata$neighborhood.composition = neighborhood.composition

theta = 1.78
gam1 = gam(biomass~s(composition, k=160) + s(x,y,k=100), data=modeldata, gamma=1.4, family=Tweedie(p=theta, link='log'))

#Get the scale (a) and the location (b)
sqrt(abs(resid(gam1, type='deviance'))) -> a
predict(gam1, type='link') -> b

#Analyze the residuals to see whether the error variance is stable across the range of the latent predictor
scatter.smooth(exp(b),a)
summary(lm(a~exp(b)))

#Analyze the residuals to see whether the error variance is stable across the range of the fitted values
dev.new()
scatter.smooth(b,a)
summary(lm(a~b))



########################################
##Two-stage: get the data for first-stage and second-stage models
indx = which(biomass.wi[,sp]>0)
stage1.data = list(indicator=ifelse(biomass.wi[,sp]>0, 1, 0), biomass=biomass.wi[,sp], composition=composition.wi[,sp], x=composition.wi[,'x'], y=composition.wi[,'y'])
stage2.data = list(biomass=biomass.wi[indx,sp], logbiomass=log(biomass.wi[indx,sp]), composition=composition.wi[indx,sp], x=composition.wi[indx,'x'], y=composition.wi[indx,'y'])

#Use the mean composition within the first-order neighborhood:
n1 = nrow(biomass.wi)
neighborhood.composition = vector()
for (i in 1:n1) {
    neighborhood.composition = c(neighborhood.composition, sum(stage1.data$composition * neighborhood[i,]) / sum(neighborhood[i,]))
}
stage1.data$neighborhood.composition = neighborhood.composition
stage2.data$neighborhood.composition = neighborhood.composition[indx]

#Make the two-stage models:
stage1.gam = gam(indicator~s(neighborhood.composition, k=50) + s(x,y,k=100), data=stage1.data, family='binomial', gamma=1.4)
stage2.gam = gam(biomass~s(neighborhood.composition, k=50) + s(x,y,k=100), data=stage2.data, family=Gamma(link='log'), gamma=1.4)

sqrt(abs(resid(stage2.gam, type='deviance'))) -> a2
predict(stage2.gam, type='link') -> b2


f = formula("biomass~composition", data=modeldata)


##Try a GAM model:
#gam1 = gam(biomass~s(composition, k=50), data=modeldata, family=Tweedie(p=theta, link='log'))
#plot(modeldata$composition, fitted(gam1), cex=0.2, pch=20, xlim=xx, ylim=yy, col='black', bty='n', ann=F, xaxt='n', yaxt='n')


##################Model the spatial effect and the composition effect as separate smooths:
##Fit the tweedie model
#gam2 = gam(biomass~s(x, y, k=200) + s(composition, k=100), gamma=1.4, data=modeldata, family=Tweedie(p=theta, link='log'))
#fits = predict(gam2, type='response')
#observed = modeldata$biomass

##Find the probability of zero biomass implied by the Tweedie model:
#lambda = 1/summary(gam2)$dispersion * fits^(2-theta) / (2-theta)
#p.zero = exp(-lambda)
#hist(exp(-lambda), breaks=30)



##################Model the spatial effect and the composition effect together in one smooth:
##Fit the tweedie model
#modeldata.small = modeldata[1:1000,]
#gam.small = gam(biomass~s(x, y, composition, k=200), data=modeldata.small, family=Tweedie(p=theta, link='log'))
#fits.small = predict(gam.small, type='response')
#observed.small = modeldata.small$biomass
#
##Find the probability of zero biomass implied by the Tweedie model:
#lambda.small = 1/summary(gam.small)$dispersion * fits.small^(2-theta) / (2-theta)
p.zero.small = exp(-lambda.small)
hist(exp(-lambda.small), breaks=30)




##Put the observed stem density and the weight into matrices that we can plot as heatmaps
##Create the matrix for the observed stem density (filled by default with NAs):
#loc = with(modeldata, list(lat=unique(y), long=unique(x)))
#pzmat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
#rownames(pzmat) <- sort(unique(modeldata$y), decreasing=F)
#colnames(pzmat) <- sort(unique(modeldata$x), decreasing=F)

##Put the observed stem densities and weights into their lat-long matrices
#for(row in 1:dim(modeldata)[1]) {
#    pzmat[as.character(modeldata$y[row]), as.character(modeldata$x[row])] = p.zero[row]
#}

##Make several plots of the model and the data:
#par(bty='n')
#gwr.matplot(pzmat, c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75")
#title(paste(sp, " biomass residual (log scale)", sep=""))
