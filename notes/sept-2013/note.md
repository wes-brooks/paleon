#PalEON biomass modeling
###Biomass data
From the survey data, Simon has calculated the "observed" biomass for each grid cell. Biomass is calculated from the measurements made at ninety-ish corner points within the grid cell. Thus, the biomass calculation reflects the state of affairs at the corner points, which are imperfect as indicators of the biomass across the whole grid cell.

###Overview
We seek a model that reflects the distribution of biomass on each grid cell. Because nearby grid cells are expected to be more similar than distant ones, we will produce a spatial smooth of the biomass. In particular, we will try using a generalized additive model (GAM) to fit the spatial smooth.

A GAM is a method of modeling the mean of the response variable via a smooth function of the covariates, for example:

\[
\eta_i = f(s_{x,i}, s_{y,i})
\]

Where $\bm{\eta}$ is a vector of linear predictors that are related through a link function to the mean of the response variable. In this case the log link is used, giving us:

\[
\log{\mu} = \eta
\]

or

\[
\text{E}(
\]


The covariates over which the GAM is to smooth the biomass are the latitude and longitude.

###Biomass distribution
Biomass is a positive, continuous quantity that is zero in grid cells without any trees. While there are relatively few grid cells where no trees were observed at all, for each individual taxon, there are . It appears to follow a power law, meaning that the variability grows exponentially with the mean.

Some distributions that are commonly used to model continuous data following a power law are the log-normal distribution and the gamma distribution with a log link. However, neither of these distributions allows for exact zeroes in the data. For this purpose, we look to a Tweedie distribution.

###Power parameter
