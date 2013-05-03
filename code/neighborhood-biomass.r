#Get an indicator of taxon presence:
modeldata$present = ifelse(modeldata$biomass>0,1,0)

#Consider only a fraction of the data:
N = 1000
modeldata.small = modeldata[1:N,]

#We need to put the neighborhood structure of the biomass data into a matrix:
library(Matrix)
n = dim(modeldata.small)[1]

Xmat = Matrix(rep(modeldata.small$x, times=n), n, n)
Ymat = Matrix(rep(modeldata.small$y, times=n), n, n)
D = sqrt((Xmat-t(Xmat))**2 + (Ymat-t(Ymat))**2)
diag(D) = Inf

#neighbors = Matrix(0, n, n)
P = vector()
for (i in 1:n) {
    cat(paste(i, "\n", sep=""))
    #neighbors[i,which(D[i,]<8800)] = 1
    P = c(P, mean(modeldata.small$composition[which(D[i,]<8800)]))
}

mm1 = glm(modeldata.small$present ~ P, family='binomial')