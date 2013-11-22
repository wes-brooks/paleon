require(MASS)
require(mgcv)
require(nlme)
require(statmod, lib.loc='R-libs')
require(tweedie, lib.loc='R-libs')
#require(tweedie)

#Modeling constants:
knots = 2000
powertol = 0.02

#Initiate timing:
ptm.tot = proc.time()

#Write to output:
sink(paste("output/", taxon, "-log.txt", sep=""), append=FALSE)
cat(paste("running for: ", taxon, '\n', sep=''))
sink()

#################################################################
#Modeling - Tweedie one-stage
#################################################################
#Set up the data, including mean composition within the first-order neighborhood:
dataset = biomass
modeldata = list(biomass=dataset[,taxon], x=dataset[,'x'], y=dataset[,'y'])
modeldata = as.data.frame(modeldata)

#MLE of theta:
powertune2 = function(theta, data, k=150) {
    cat(paste("Tuning. knots=", k, "; theta: ", round(theta, 3), sep=''))
    
    #Make the model and export it (so we dont have to reproduce it after optimization)
    model = gam(biomass~s(x,y,k=k), data=data, gamma=1.4, family=Tweedie(p=theta, link='log'))
    assign('model.out', model, envir=.GlobalEnv)

    ll = -logLik(model)[1]
    cat(paste("; -logLik: ", round(ll,3), "\n", sep=''))
	return(ll)
}
#Initialize timing for this iteration:
ptm = proc.time()

#Locate the optimal theta. The optimization also exports the optimal model as object 'model.out'
sink(paste("output/", taxon, "-log.txt", sep=""), append=TRUE)
tuning = optimize(powertune2, interval=c(1,2), data=modeldata, k=knots, tol=powertol)
cat(paste("\ntheta: ", round(tuning$minimum, 3), ", slope: ", round(tuning$objective,4), '\n', sep=''))
mle <- model.out

#Write the analysis of knots to disk:
print(mgcv:::k.check(mle))
sink()

#Finalize timing for this iteration:
sink(paste("output/", taxon, "-log.txt", sep=""), append=TRUE)
print(proc.time() - ptm)
sink()

###################
#Draw from the "posterior" distribution of biomass:
###################
#Get the matrix that is multiplied by the coefficients to compute the linear predictor:
Xp = predict(mle, type='lpmatrix')

#Initialize the parametric bootstrap with the actual data and the MLEs of beta and mu:
data.boot = modeldata
fit = fitted(mle)
beta.resampled = matrix(mvrnorm(n=100, coef(mle), mle$Vp), nrow=100, ncol=length(coef(mle)))

#Evaluate the linear predictors for each draw of the coefficients and write to disk:
lp = Xp %*% t(beta.resampled)
write.table(t(lp), file=paste("output/logbiomass-", cluster, "-", taxon, ".csv", sep=""),
    append=FALSE, row.names=FALSE, col.names=FALSE, sep=',')
lp = NULL
beta.resampled = NULL

#Initialize structures to hold the parametric bootstrap estimates:
smoothing.params = c(mle$sp)
theta = c(tuning$minimum)
s2 = c(mle$sig2)

#Resample from the model (this is the parametric bootstrap):
S = 19
for (i in 1:S) {
    #Initialize timing for this iteration:
    ptm = proc.time()

    #Regenerate the output via the parametric model:
    y = rtweedie(nrow(modeldata), mu=fit, phi=mle$sig2, power=tuning$minimum)
    data.boot$biomass = y

    #Tune the model on the regenerated data
    sink(paste("output/", taxon, "-log.txt", sep=""), append=TRUE)
    tuning.boot = optimize(powertune2, interval=c(1,2), data=data.boot, k=knots, tol=powertol)
    sp.boot = model.out$sp
    theta.boot = tuning.boot$minimum
    
    #Run the model on the original data, using the parameters from the bootstrap:
    m.boot = gam(biomass~s(x,y,k=knots), data=modeldata,
            gamma=1.4, sp=sp.boot,
            family=Tweedie(p=theta.boot, link='log'))
    print(mgcv:::k.check(m.boot))
    sink()
    
    #Draw a bunch of spline coefficients from this estimate of their distribution:
    beta.resampled = mvrnorm(n=100, coef(m.boot), m.boot$Vp)
    
    #Evaluate the linear predictors for each draw of the coefficients and write to disk:
    lp = Xp %*% t(beta.resampled)
    write.table(t(lp), file=paste("output/logbiomass-", cluster, "-", taxon, ".csv", sep=""),
        append=TRUE, row.names=FALSE, col.names=FALSE, sep=',')
    lp = NULL
    beta.resampled = NULL

    #Add this round of parameters to the output:
    s2 = c(s2, m.boot$sig2)
    smoothing.params = c(smoothing.params, sp.boot)
    theta = c(theta, theta.boot)
    
    #Finalize timing for this iteration:
    sink(paste("output/", taxon, "-log.txt", sep=""), append=TRUE)
    print(proc.time() - ptm)
    sink()
}

#Write the parameters to disk:
params = as.data.frame(list(s2, smoothing.params, theta))
write.table(params, file=paste("output/params-", cluster, "-", taxon, ".csv", sep=""),
    append=FALSE, row.names=FALSE, col.names=FALSE, sep=',')

#Finalize timing:
sink(paste("output/", taxon, "-log.txt", sep=""), append=TRUE)
print(proc.time() - ptm.tot)
sink()
