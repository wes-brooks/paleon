setwd("~/git/paleon")

#load necessary R libraries
library(mgcv)

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
presence.models = list()
modeldata = list(indicator=ifelse(composition[,sp]==0,-1,1), composition=composition[,sp], biomass=biomass[,sp], x=composition$x, y=composition$y)
paleon.dist = dist(cbind(modeldata$x, modeldata$y))
pd2 = as.matrix(paleon.dist)
diag(pd2) = Inf

automodeldata = data.frame()

for (i in 1:dim(pd2)[1]) {
    neighbors = which(pd2[i,] == min(pd2[i,]))
    automodeldata = rbind(automodeldata, c(modeldata$indicator[i], modeldata$x[i], modeldata$y[i], mean(modeldata$indicator[neighbors]), mean(modeldata$composition[neighbors])))
}
colnames(automodeldata) = c('indicator', 'x', 'y', 'neighborhood.indicator', 'neighborhood.composition')
automodeldata$indicator = ifelse(automodeldata$indicator==-1,0,1)

automodel1 = glm(indicator~neighborhood.indicator, data=automodeldata, family='binomial')
automodel2 = glm(indicator~neighborhood.composition, data=automodeldata, family='binomial')
