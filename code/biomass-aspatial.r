#load necessary R libraries
library(mgcv)
library(statmod)
library(tweedie)
library(Matrix)
library(plotrix)



#Import the data and the heatmap code.
source("code/biomass-import.r")
taxa = c('Cherries', 'Willow', 'Walnuts', 'Hickory', 'Beech', 'Fir', 'Spruce', 'Ironwoods', 'Cedar', 'Hemlock', 'Basswood', 'Ashes', 'Elms', 'Poplar', 'Pine', 'Tamarack', 'Birches', 'Maple', 'Oaks')
#args <- commandArgs(trailingOnly = TRUE)
#indx = as.numeric(args[1])
indx = 5
sp = taxa[indx]

#Establish teh file for output
#sink(paste("output/", sp, ".txt", sep=""))
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
modeldata = list(indicator=ifelse(biomass.wi[,sp]>0, 1, 0), stemdensity=stemdens.wi[,'stem.density'], biomass=biomass.wi[,sp], logbiomass=log(biomass.wi[,sp]), composition=composition.wi[,sp], x=composition.wi[,'x'], y=composition.wi[,'y'])
neighborhood.composition = vector()
for (i in 1:n) {
    neighborhood.composition = c(neighborhood.composition, sum(modeldata$composition * neighborhood[i,]) / sum(neighborhood[i,]))
}
modeldata$neighborhood.composition = neighborhood.composition

#Function to set the optimial Tweedie theta:
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

#Now select the optimal Tweedie parameter:
tuning = optimize(bm.opt, interval=c(1,2), data=modeldata, k=150, tol=0.01)
print("Result:")
print(tuning)
print("Good night, Irene")