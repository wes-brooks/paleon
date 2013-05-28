library(tweedie)
library(mgcv)

#Define parameters
n = 500
phi = 3
thetas = c(1, 1.5, 2)
B = 1
N = 100
mu = 3
pdf("figures/tweedie-sim.pdf", width=12, height=6)
layout(t(matrix(1:3)))
for (theta in thetas) {
    Y = rtweedie(n=n, phi=phi, mu=mu, power=theta)
    hist(Y, breaks=20)
}
dev.off()

p.f = function(p) {
    m = gam(Y~X, family=Tweedie(p=p, link='log'))
    
    a = abs(resid(m, type='deviance'))
    b = predict(m, type='response')
    
    m2 = lm(a~b)
    return(coef(m2)[2]**2)
}

results = vector()
for (i in 1:N) {
    cat(paste("Starting simulation ", i, ".\n", sep=""))

    #Simulate data
    X = 5*runif(n, min=-1, max=1)
    eps = rnorm(n)
    W = X*B
    Z = exp(W)
    Y = rtweedie(n=n, phi=phi, mu=Z, power=theta)

    opt = optimize(p.f, lower=1, upper=2, tol=0.01)
    results = c(results, opt$minimum)
    cat(paste("Found a minimum at theta=", opt$minimum, ".\n", sep=""))
}


