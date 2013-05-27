library(tweedie)
library(statmod)
library(mgcv)

#######################################
#Test a method for estimating p (based on matching deviance residuals to the variance function)
#######################################

#Establish parameters
phi = 1
b = 2
n = 1000
sigma = 0.5
p = 1.9
N = 200

result = vector()

for (i in 1:N) {
	#Simulate data
	X = rnorm(n)
	U = X*b
	Y = rtweedie(n=n, phi=phi, mu=exp(U), p=p)


	#opt.glm = function(p) {
	#   glm(Y~X, family=tweedie(var.power=p, link.power=0)) -> m    
	#   lm1 = lm(abs(resid(m))~predict(m))
	#   return(coef(lm1)[2]**2)
	#}

	#optimize(opt.glm, lower=1, upper=2)

	
	opt.gam = function(p) {
	   glm(Y~X, family=tweedie(var.power=p, link.power=0)) -> m    
	   lm1 = lm(abs(resid(m))~predict(m))
	   return(coef(lm1)[2]**2)
	}

	optimize(opt.gam, lower=1, upper=2, tol=0.005) -> o1
	cat(paste("For i=", i, ", minimum is at p=", round(o1$minimum,3), "\n", sep=""))
	result = c(result, o1$minimum)
}

png("~/Desktop/p-plot.png")
hist(result)
dev.off()