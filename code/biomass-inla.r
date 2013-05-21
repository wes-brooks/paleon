#load necessary R libraries
library(mgcv)
library(INLA)
library(tweedie)
library(Matrix)
library(plotrix)



#Import the data and the heatmap code.
source("code/biomass-import.r")
taxa = c('Cherries', 'Willow', 'Walnuts', 'Hickory', 'Beech', 'Fir', 'Spruce', 'Ironwoods', 'Cedar', 'Hemlock', 'Basswood', 'Ashes', 'Elms', 'Poplar', 'Pine', 'Tamarack', 'Birches', 'Maple', 'Oaks')
#args <- commandArgs(trailingOnly = TRUE)
#indx = as.numeric(args[1])
indx = 2
sp = taxa[indx]


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
neighborhood.composition = vector()
for (i in 1:n) {
    neighborhood.composition = c(neighborhood.composition, sum(modeldata$composition * neighborhood[i,]) / sum(neighborhood[i,]))
}
modeldata$neighborhood.composition = neighborhood.composition
modeldata$log.neighborhood.composition = ifelse(neighborhood.composition==0,NA,log(neighborhood.composition))

mi1 = inla(biomass.nz ~ f(loc, model='besag', graph=graphfile), data=modeldata, verbose=TRUE, family='Gamma')
mi2 = inla(biomass.nz ~ f(loc, model='besag', graph=graphfile) + f(neighborhood.composition, model='rw1'), data=modeldata, verbose=TRUE, family='Gamma')
