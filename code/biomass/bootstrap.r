library(MASS)
library(mgcv)
library(nlme)
library(statmod, lib.loc='R-libs')
library(tweedie, lib.loc='R-libs')
library(Matrix)

#sink(paste("output/", taxon, "-log.txt", sep=""))
cat(paste("running for: ", taxon, '\n', sep=''))

#Modeling constants:
knots = 250
powertol = 0.02

#################################################################
#Modeling - Tweedie one-stage
#################################################################
#Set up the data, including mean composition within the first-order neighborhood:
dataset = biomass.wi
modeldata = list(biomass=dataset[,taxon], x=dataset[,'x'], y=dataset[,'y'])
modeldata = as.data.frame(modeldata)

#Function to set the optimial Tweedie theta:
powertune = function(theta, data, k=150) {
    #Make the model and export it (so we dont have to reproduce it after optimization)
    model = gam(biomass~s(x,y,k=k), data=data, gamma=1.4, family=Tweedie(p=theta, link='log'))
    assign('model.out', model, envir=.GlobalEnvironment)

	#Get the scale and the location
	scale <- sqrt(abs(resid(model, type='deviance')))
	loc <- predict(model, type='link')

    #We are looking for the power parameter that gives us zero slope for the scale-location plot:
	m = lm(scale~loc)
    cat(paste("Tuning. theta: ", round(theta, 3), ", slope: ", round(abs(coef(m)[2]),4), '\n', sep=''))
	return(abs(coef(m)[2]))
}

#Locate the optimal theta. The optimization also exports the optimal model as object 'model.out'
tuning = optimize(powertune, interval=c(1,2), data=modeldata, k=knots, tol=powertol)
cat(paste("\ntheta: ", round(tuning$minimum, 3), ", slope: ", round(tuning$objective,4), '\n', sep=''))
mle <- model.out

###################
#Draw from the "posterior" distribution of biomass:
###################
#Get the matrix that is multiplied by the coefficients to compute the linear predictor:
Xp = predict(mle, type='lpmatrix')

#Initialize the parametric bootstrap with the actual data and the MLEs of beta and mu:
data.boot = modeldata
fit = fitted(mle)
beta.resampled = Matrix(mvrnorm(n=100, coef(mle), mle$Vp), nrow=100, ncol=length(coef(mle)))

#Initialize structures to hold the parametric bootstrap estimates:
smoothing.params = c(mle$sp)
theta = c(tuning$minimum)
s2 = c(mle$sig2)

#Resample from the model (this is the parametric bootstrap):
S = 9
for (i in 1:S) {
    #Regenerate the output via the parametric model:
    y = rtweedie(nrow(modeldata), mu=fit, phi=mle$sig2, power=tuning$minimum)
    data.boot$biomass = y

    #Tune the model on the regenerated data
    tuning.boot = optimize(powertune, interval=c(1,2), data=data.boot, k=knots, tol=powertol)
    sp.boot = model.out$sp
    theta.boot = tuning.boot$minimum
    
    #Run the model on the original data, using the parameters from the bootstrap:
    m.boot = gam(biomass~s(x,y,k=knots), data=modeldata,
            gamma=1.4, sp=sp.boot,
            family=Tweedie(p=theta.boot, link='log'))

    #Draw a bunch of spline coefficients from this estimate of their distribution:
    beta.resampled = rbind(beta.resampled, mvrnorm(n=100, coef(m.boot), m.boot$Vp))
    
    #Add this round of parameters to the output:
    s2 = c(s2, m.boot$sig2)
    smoothing.params = c(smoothing.params, sp.boot)
    theta = c(theta, theta.boot)
}

#Evaluate the linear predictors for each draw of the coefficients:
lp = Xp %*% t(beta.resampled)

#Write the draws of the linear predictors to disk:
write(lp, file = paste("output/logbiomass-", cluster, "-", taxon, ".csv", sep=""),
    ncolumns = ncol(lp),
    append = FALSE,
    sep = ",")

#Write the parameters to disk:
params = as.data.frame(list(s2, smoothing.params, theta))
write(params, file = paste("output/params-", cluster, "-", taxon, ".csv", sep=""),
    ncolumns = ncol(3),
    append = FALSE,
    sep = ",")
