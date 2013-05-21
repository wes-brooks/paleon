#################################################################
#Imports (data and libraries)
#################################################################
setwd("~/git/paleon")

#load necessary R libraries
library(mgcv)
library(tweedie)
library(Matrix)
library(plotrix)

#Import the data and the heatmap code.
source("~/git/brooks/code/matplot.r")
source("code/biomass-import.r")
sp = "Oaks"


#################################################################
#Data - neighborhood structure (first-order)
#################################################################
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




#################################################################
#Modeling - Tweedie
#################################################################
#Make the one-stage model
#Set up the data, including mean composition within the first-order neighborhood:
modeldata = list(indicator=ifelse(biomass.wi[,sp]>0, 1, 0), biomass=biomass.wi[,sp], logbiomass=log(biomass.wi[,sp]), composition=composition.wi[,sp], x=composition.wi[,'x'], y=composition.wi[,'y'])
md = as.data.frame(modeldata)
md$logbiomass[which(md$logbiomass==-Inf)] = NaN

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



#################################################################
#Modeling - Tweedie (using neighborhood-mean composition)
#################################################################
#Make the one-stage model
#Set up the data, including mean composition within the first-order neighborhood:
theta2 = 1.38
gam2 = gam(biomass~s(neighborhood.composition, k=25) + s(x,y,k=150), data=modeldata, gamma=1.4, family=Tweedie(p=theta2, link='log'))

#Get the scale (a) and the location (b)
sqrt(abs(resid(gam2, type='deviance'))) -> a2
predict(gam2, type='link') -> b2

#Analyze the residuals to see whether the error variance is stable across the range of the latent predictor
dev.new()
scatter.smooth(exp(b2),a2)
summary(lm(a2~exp(b2)))

#Analyze the residuals to see whether the error variance is stable across the range of the fitted values
dev.new()
scatter.smooth(b2,a2, main="Residual location-scale plot", bty='n')
summary(lm(a2~b2))



#################################################################
#Modeling - two-stage (delta)
#################################################################
#Two-stage: get the data for first-stage and second-stage models
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

sqrt(abs(resid(stage2.gam, type='deviance'))) -> a3
predict(stage2.gam, type='link') -> b3
dev.new()
scatter.smooth(b3,a3)




#################################################################
#Plotting - Heatmap(observed)
#################################################################
#Put the observed biomass into matrices that we can plot as heatmaps
#Create the matrix for the observed stem density (filled by default with NAs):
loc = with(stage2.data, list(lat=unique(y), long=unique(x)))
biomat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
rownames(biomat) <- sort(unique(stage2.data$y), decreasing=F)
colnames(biomat) <- sort(unique(stage2.data$x), decreasing=F)

#Put the observed biomass into its lat-long matrices
for(row in 1:length(stage2.data[['x']])) {
    biomat[as.character(stage2.data$y[row]), as.character(stage2.data$x[row])] = stage2.data$logbiomass[row]
}

#Plot the model and the data:
par(bty='n')
gwr.matplot(biomat, c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75")
title(paste(sp, " (observed): biomass (log scale)", sep=""))



#################################################################
#Plotting - Heatmap(fitted - Tweedie, E[biomass])
#################################################################
#Put the observed biomass into matrices that we can plot as heatmaps
#Create the matrix for the observed stem density (filled by default with NAs):
loc = with(modeldata, list(lat=unique(y), long=unique(x)))
biomat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
rownames(biomat) <- sort(unique(modeldata$y), decreasing=F)
colnames(biomat) <- sort(unique(modeldata$x), decreasing=F)
output = predict(gam1, type='link')

#Put the observed biomass into its lat-long matrices
for(row in 1:length(modeldata[['x']])) {
    biomat[as.character(modeldata$y[row]), as.character(modeldata$x[row])] = output[row]
}

#Plot the model and the data:
par(bty='n')
gwr.matplot(biomat, c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75")
title(paste(sp, " (Tweedie model): E[biomass] (log scale)", sep=""))



#################################################################
#Plotting - Heatmap(fitted - Tweedie, P[biomass=0])
#################################################################
#Put the probability of zero into matrices that we can plot as heatmaps
loc = with(modeldata, list(lat=unique(y), long=unique(x)))
pzmat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
rownames(pzmat) <- sort(unique(modeldata$y), decreasing=F)
colnames(pzmat) <- sort(unique(modeldata$x), decreasing=F)

phi = summary(gam1)$scale
mu = predict(gam1, type='response')
p = theta
lambda = mu^(2-p) / (2-p) / phi
output = -lambda - log(1-exp(-lambda))

#Put the probability of zero into its lat-long matrices
for(row in 1:length(modeldata[['x']])) {
    pzmat[as.character(modeldata$y[row]), as.character(modeldata$x[row])] = output[row]
}

#Plot the model and the data:
par(bty='n')
gwr.matplot(pzmat, c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75")
title(paste(sp, " (Tweedie model): logit[P(biomass=0)]", sep=""))



#################################################################
#Plotting - Heatmap(fitted - Tweedie2, E[biomass])
#################################################################
#Put the observed biomass into matrices that we can plot as heatmaps
#Create the matrix for the observed stem density (filled by default with NAs):
loc = with(modeldata, list(lat=unique(y), long=unique(x)))
biomat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
rownames(biomat) <- sort(unique(modeldata$y), decreasing=F)
colnames(biomat) <- sort(unique(modeldata$x), decreasing=F)
output = predict(gam2, type='link')

#Put the observed biomass into its lat-long matrices
for(row in 1:length(modeldata[['x']])) {
    biomat[as.character(modeldata$y[row]), as.character(modeldata$x[row])] = output[row]
}

#Plot the model and the data:
par(bty='n')
gwr.matplot(biomat, c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75")
title(paste(sp, " (Tweedie model 2): E[biomass] (log scale)", sep=""))


#################################################################
#Plotting - Heatmap(fitted - Tweedie2 (neighborhood composition), P[biomass=0])
#################################################################
#Put the probability of zero into matrices that we can plot as heatmaps
loc = with(modeldata, list(lat=unique(y), long=unique(x)))
pzmat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
rownames(pzmat) <- sort(unique(modeldata$y), decreasing=F)
colnames(pzmat) <- sort(unique(modeldata$x), decreasing=F)

phi = summary(gam2)$scale
mu = predict(gam2, type='response')
p = theta2
lambda = mu^(2-p) / (2-p) / phi
output = -lambda - log(1-exp(-lambda))

#Put the probability of zero into its lat-long matrices
for(row in 1:length(modeldata[['x']])) {
    pzmat[as.character(modeldata$y[row]), as.character(modeldata$x[row])] = output[row]
}

#Plot the model and the data:
par(bty='n')
gwr.matplot(pzmat, c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75")
title(paste(sp, " (Tweedie model 2): logit[P(biomass=0)]", sep=""))




#################################################################
#Plotting - Heatmap(fitted - Delta, P[biomass=0])
#################################################################
#Put the probability of zero into matrices that we can plot as heatmaps
loc = with(stage1.data, list(lat=unique(y), long=unique(x)))
pzmat.1 = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
rownames(pzmat.1) <- sort(unique(stage1.data$y), decreasing=F)
colnames(pzmat.1) <- sort(unique(stage1.data$x), decreasing=F)

output = -predict(stage1.gam, type="link")

#Put the probability of zero into its lat-long matrices
for(row in 1:length(stage1.data[['x']])) {
    pzmat.1[as.character(stage1.data$y[row]), as.character(stage1.data$x[row])] = output[row]
}

#Plot the model and the data:
par(bty='n')
gwr.matplot(pzmat.1, c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75")
title(paste(sp, " (Delta model): logit[P(biomass=0)]", sep=""))



#################################################################
#Plotting - Heatmap(fitted - Tweedie2 v Delta, P[biomass=0])
#################################################################
#Put the probability of zero into matrices that we can plot as heatmaps
phi = summary(gam2)$scale
mu = predict(gam2, type='response')
p = theta2
lambda = mu^(2-p) / (2-p) / phi
output1 = -lambda - log(1-exp(-lambda))
output2 = -predict(stage1.gam, type="link")

plot(output1, output2, main="logit[P(biomass=0)], Tweedie vs. delta", xlab="Tweedie", ylab="delta", bty='n')
abline(a=0,b=1)


#################################################################
#Plotting - Heatmap(fitted - Tweedie2 v Delta, biomass)
#################################################################
#Put the probability of zero into matrices that we can plot as heatmaps
output1 = predict(gam2, type='link')
output2 = predict(stage2.gam, type="link")

plot(output1[indx], output2, main="Oak biomass, Tweedie vs. delta", xlab="Tweedie", ylab="delta", bty='n')
abline(a=0,b=1)


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
#p.zero.small = exp(-lambda.small)
#hist(exp(-lambda.small), breaks=30)

