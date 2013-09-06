library(mgcv)

#MLE power parameter
p.sel.1 = function(p, X, Y) {
    deviance(gam(Y~X, family=Tweedie(link='log', p=p)))
}

#Slope-based power parameter
p.sel.2 = function(p, X, Y) {
    mod = gam(Y~X, family=Tweedie(link='log', p=p))
    sqrt(abs(resid(mod, type='deviance'))) -> scale
    predict(mod, type='link') -> loc
    m = lm(scale~loc)

    return(abs(coef(m)[2]))
}



n = 200
b = 5

#Simulation settings:
p = rep(seq(1.1,1.9,length.out=5), each=4)
phi = rep(c(rep(1,2),rep(10,2)),5)
Int = rep(c(0,5),10)

#Vectors to store the outputs:
t.mle = vector()
t.slope = vector()

#Run the simulations:
for (i in 1:200) {
    setting = ((i-1) %% 20) + 1

    pow = p[setting]
    disp = phi[setting]

    cat(paste("Run: ", i, ", setting: ", setting, ", intercept: ", Int[setting], ", power: ", pow, ", dispersion: ", disp, "\n", sep=""))

    x = rnorm(n, 0, 1)
    mu = exp(Int[setting]+x*b)
    y = rTweedie(mu=mu, p=pow, phi=disp)

    t.mle = c(t.mle, optimize(p.sel.1, interval=c(1,2), X=x, Y=y)$minimum)
    t.slope = c(t.slope, optimize(p.sel.2, interval=c(1,2), X=x, Y=y)$minimum)
}
