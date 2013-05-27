#load necessary R libraries
library(mgcv)
library(tweedie)
library(Matrix)
library(plotrix)



#Import the data and the heatmap code.
source("code/biomass-import.r")
taxa = c('Cherries', 'Willow', 'Walnuts', 'Hickory', 'Beech', 'Fir', 'Spruce', 'Ironwoods', 'Cedar', 'Hemlock', 'Basswood', 'Ashes', 'Elms', 'Poplar', 'Pine', 'Tamarack', 'Birches', 'Maple', 'Oaks')
args <- commandArgs(trailingOnly = TRUE)
indx = as.numeric(args[1])
sp = taxa[indx]

#Establish teh file for output
sink(paste("output/", sp, ".txt", sep=""))
print("running for: ")
print(sp)

#################################################################
#Data - neighborhood structure (first-order)
#################################################################
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


#################################################################
#Modeling - Tweedie
#################################################################
#Make the one-stage model
#Set up the data, including mean composition within the first-order neighborhood:
modeldata = list(indicator=ifelse(biomass.wi[,sp]>0, 1, 0), biomass=biomass.wi[,sp], logbiomass=log(biomass.wi[,sp]), composition=composition.wi[,sp], x=composition.wi[,'x'], y=composition.wi[,'y'])
neighborhood.composition = vector()
for (i in 1:n) {
    neighborhood.composition = c(neighborhood.composition, sum(modeldata$composition * neighborhood[i,]) / sum(neighborhood[i,]))
}
modeldata$neighborhood.composition = neighborhood.composition


k1 = 25
k2 = 150

#Function to set the optimial Tweedie theta:
bm.opt = function(theta, data, k1=25, k2=150) {
	result = list()

	data = as.data.frame(data)
	model = gam(biomass~s(neighborhood.composition, k=k1) + s(x,y,k=k2), data=data, gamma=1.4, family=Tweedie(p=theta, link='log'))
	cat(paste("theta: ", theta, '\n', sep=''))
	print(summary(model))
	print(logLik(model))

	#Get the scale (a) and the location (b)
	sqrt(abs(resid(model, type='deviance'))) -> scale
	predict(model, type='link') -> loc
	m = lm(scale~loc)
	print(summary(m))

	return(abs(coef(m)[2]))
}


#Make the one-stage model
#Set up the data, including mean composition within the first-order neighborhood:
modeldata = list(indicator=ifelse(biomass.wi[,sp]>0, 1, 0), biomass=biomass.wi[,sp], logbiomass=log(biomass.wi[,sp]), composition=composition.wi[,sp], x=composition.wi[,'x'], y=composition.wi[,'y'])
neighborhood.composition = vector()
for (i in 1:n) {
	neighborhood.composition = c(neighborhood.composition, sum(modeldata$composition * neighborhood[i,]) / sum(neighborhood[i,]))
}
modeldata$neighborhood.composition = neighborhood.composition

#Now select the optimal Tweedie parameter:
tuning = optimize(bm.opt, interval=c(1,2), data=modeldata, k1=k1, k2=k2, tol=0.008)
print("Result:")
print(tuning)
print("Good night, Irene")