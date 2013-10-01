#Imports
library(tweedie)

#Simulation settings
N = 100

#True parameters
B = 5
p = 1.2
phi = 1

#Simulate data
X = runif(N)
eta = X*B
mu = exp(eta)
Y = rtweedie(n=N, mu=mu, phi=phi, power=p)
d = list(X, Y)

#Estimate model with misspecified p
m2 = glm(Y~X, data=d, family=tweedie(var.power=1.8, link.power=0))
m3 = glm(Y~X, data=d, family=tweedie(var.power=1.5, link.power=0))
m4 = glm(Y~X, data=d, family=tweedie(var.power=1.05, link.power=0))
m = glm(Y~X, data=d, family=tweedie(var.power=1.2, link.power=0))

scatter.smooth(x=fitted(m), residuals(m, type='deviance')**2)