f = formula("biomass~composition", data=modeldata)

eval.tweedie = function(theta) {
    m = glm(f, data=modeldata, family=tweedie(var.power=theta, link.power=0), control=glm.control(maxit=100))
    return(sum(residuals(m, type="deviance")**2))
}

opt = optimize(eval.tweedie, interval=c(1,2))


#Try a GAM model:
library(mgcv)
gam1 = gam(biomass~s(composition, k=50), data=modeldata, family=Tweedie(p=theta, link='log'))
plot(modeldata$composition, fitted(gam1), cex=0.2, pch=20, xlim=xx, ylim=yy, col='black', bty='n', ann=F, xaxt='n', yaxt='n')


##################Model the spatial effect and the composition effect as separate smooths:
#Fit the tweedie model
gam2 = gam(biomass~s(x, y, k=200) + s(composition, k=100), data=modeldata, family=Tweedie(p=theta, link='log'))
fits = predict(gam2, type='response')
observed = modeldata$biomass

#Find the probability of zero biomass implied by the Tweedie model:
lambda = 1/summary(gam2)$dispersion * fits^(2-theta) / (2-theta)
p.zero = exp(-lambda)
hist(exp(-lambda), breaks=30)



##################Model the spatial effect and the composition effect together in one smooth:
#Fit the tweedie model
gam2 = gam(biomass~s(x, y, k=200) + s(composition, k=100), data=modeldata, family=Tweedie(p=theta, link='log'))
fits = predict(gam2, type='response')
observed = modeldata$biomass

#Find the probability of zero biomass implied by the Tweedie model:
lambda = 1/summary(gam2)$dispersion * fits^(2-theta) / (2-theta)
p.zero = exp(-lambda)
hist(exp(-lambda), breaks=30)




#Put the observed stem density and the weight into matrices that we can plot as heatmaps
#Create the matrix for the observed stem density (filled by default with NAs):
loc = with(modeldata, list(lat=unique(y), long=unique(x)))
pzmat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
rownames(pzmat) <- sort(unique(modeldata$y), decreasing=F)
colnames(pzmat) <- sort(unique(modeldata$x), decreasing=F)

#Put the observed stem densities and weights into their lat-long matrices
for(row in 1:dim(modeldata)[1]) {
    pzmat[as.character(modeldata$y[row]), as.character(modeldata$x[row])] = p.zero[row]
}

#Make several plots of the model and the data:
par(bty='n')
gwr.matplot(pzmat, c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, na.color="grey75")
title(paste(sp, " biomass residual (log scale)", sep=""))
