library(tweedie)
library(mgcv)

#Define parameters
n = 500
phi = 1
theta = 1.01
B = 1
N = 100

p.f = function(p) {
    m = gam(Y~X, family=Tweedie(p=p, link='log'))
    #scatter.smooth(fitted(m), sqrt(abs(resid(m))))
    
    #a = abs(resid(m, type='deviance'))
    b = predict(m, type='response')
    
    return(sum(dtweedie(Y, mu=b, power=p, phi=1/m$scale)))
    
    #m2 = lm(a~b)
    #return(abs(coef(m2)[2]))
}

results = vector()
for (i in 1:N) {
    cat(paste("Starting simulation ", i, ".\n", sep=""))

    #Simulate data
    X = 5*runif(n, min=-1, max=1)
    eps = rnorm(n)
    W = X*B + eps
    Z = exp(W)
    Y = rtweedie(n=n, phi=phi, mu=Z, power=theta)

    opt = optimize(p.f, lower=1, upper=2, tol=0.01)
    results = c(results, opt$minimum)
    cat(paste("Found a minimum at theta=", opt$minimum, ".\n", sep=""))
}