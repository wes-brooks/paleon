library(mgcv)

powertol = 0.01
S = 200



t1 = vector()
t2 = vector()

#Function to set the optimial Tweedie theta:
powertune = function(theta, data, k=150) {
    #Make the model and export it (so we dont have to reproduce it after optimization)
    model = gam(y~X, data=data, family=Tweedie(p=theta, link='log'))
    assign('model.out', model, envir=.GlobalEnv)

	#Get the scale and the location
	scale <- sqrt(abs(resid(model, type='deviance')))
	loc <- predict(model, type='link')

    #We are looking for the power parameter that gives us zero slope for the scale-location plot:
	m = lm(scale~loc)
    cat(paste("Tuning. theta: ", round(theta, 3), ", slope: ", round(abs(coef(m)[2]),4), '\n', sep=''))
	return(abs(coef(m)[2]))
}


#Function to set the optimial Tweedie theta:
powertune2 = function(theta, data, k=150) {
    #Make the model and export it (so we dont have to reproduce it after optimization)
    model = gam(y~X, data=data, family=Tweedie(p=theta, link='log'))
    assign('model.out', model, envir=.GlobalEnv)

    ll = -logLik(model)[1]
    cat(paste(ll, "\n", sep=''))
	return(ll)
}

for (i in 1:S) {
    X = rnorm(100)
    B = 1
    nu = exp(X*B)
    y = rTweedie(mu=nu, p=1.8)

    modeldata = as.data.frame(cbind(X,y))

    #Locate the optimal theta. The optimization also exports the optimal model as object 'model.out'
    tuning = optimize(powertune, interval=c(1,2), data=modeldata, k=knots, tol=powertol)
    cat(paste("\ntheta: ", round(tuning$minimum, 3), ", slope: ", round(tuning$objective,4), '\n', sep=''))

    t1 = c(t1,tuning$minimum)


    tuning = optimize(powertune2, interval=c(1,2), data=modeldata, k=knots, tol=powertol)
    cat(paste("\ntheta: ", round(tuning$minimum, 3), ", logLik: ", round(tuning$objective,4), '\n', sep=''))

    t2 = c(t2,tuning$minimum)
}

