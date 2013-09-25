#These functions import code and data directly from the github servers:
source("../brooks/code/source_https.r")
source("../brooks/code/load_https.r")

#load necessary R libraries
library(MASS)
library(mgcv)
library(nlme)
library(tweedie)
library(Matrix)
library(plotrix)

#Import the data and the heatmap code.
source_https("https://raw.github.com/wesesque/paleon/master/code/biomass-import.r")
taxa = c('Cherries', 'Willow', 'Walnuts', 'Hickory', 'Beech', 'Fir', 'Spruce', 'Ironwoods', 'Cedar', 'Hemlock', 'Basswood', 'Ashes', 'Elms', 'Poplar', 'Pine', 'Tamarack', 'Birches', 'Maple', 'Oaks')
#args <- commandArgs(trailingOnly = TRUE)
#indx = as.numeric(args[1])
#sp = taxa[indx]
sp='tot'

#Establish the file for output
#sink(paste("output/", sp, ".txt", sep=""))
cat(paste("running for: ", sp, '\n' sep=''))

#################################################################
#Modeling - Tweedie
#################################################################
#Make the one-stage model
#Set up the data, including mean composition within the first-order neighborhood:
modeldata = list(biomass=biomass.wi[,sp], x=composition.wi[,'x'], y=composition.wi[,'y'])

#Function to set the optimial Tweedie theta:
bm.opt = function(theta, data, k=150) {
	result = list()

	data = as.data.frame(data)
	model = gam(biomass~s(x,y,k=k), data=data, gamma=1.4, family=Tweedie(p=theta, link='log'))

	#Get the scale (a) and the location (b)
	sqrt(abs(resid(model, type='deviance'))) -> scale
	predict(model, type='link') -> loc
	m = lm(scale~loc)

    cat(paste("Tuning. theta: ", round(theta, 3), ", slope: ", round(abs(coef(m)[2]),4), '\n', sep=''))
	return(abs(coef(m)[2]))
}

k=250
tuning = optimize(bm.opt, interval=c(1,2), data=modeldata, k=k, tol=0.01)
cat(paste("\ntheta: ", round(tuning$minimum, 3), ", slope: ", round(tuning$objective,4), '\n', sep=''))

#Produce a model with the 'optimal' theta:
bm = gam(biomass~s(x,y,k=k), data=modeldata, gamma=1.4, family=Tweedie(p=tuning$minimum, link='log'))

###################
#Draw from the "posterior" distribution of biomass:
###################
Xp = predict(bm, type='lpmatrix')
br = mvrnorm(n=n, coef(bm), bm$Vp)

f = fitted(bm)
modeldata.bs = modeldata
n = nrow(biomass.wi)

mean.biomass = vector()
br = matrix(0,nrow=0, ncol=length(coef(bm)))
ss = list()
theta = vector()
s2 = vector()

S = 19
for (i in 1:S) {
    y = rtweedie(n, mu=f, phi=bm$sig2, power=tuning$minimum)
    modeldata.bs$biomass = y

    tuning.bs = optimize(bm.opt, interval=c(1,2), data=modeldata.bs, k=k, tol=0.01)
    sp = gam(biomass~s(x,y,k=k), data=modeldata.bs, gamma=1.4, family=Tweedie(p=tuning.bs$minimum, link='log'))$sp
    

    bm2 = gam(biomass~s(x,y,k=k), data=modeldata, gamma=1.4, sp=sp, family=Tweedie(p=tuning.bs$minimum, link='log'))

    br = rbind(br, mvrnorm(n=100, coef(bm2), bm2$Vp))

    s2 = c(s2, bm2$sig2)
    ss[[i]] = sp
    theta = c(theta, tuning.bs$minimum)
}

#Add some draws from the original model:
br = rbind(br, mvrnorm(n=100, coef(bm), bm$Vp))

#Get the estimated total mean biomass for each simulation:
lp = Xp %*% t(br)
mean.biomass.2 = colSums(exp(lp))
