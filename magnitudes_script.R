# import libraries
library(ggplot2)
library(INLA)
library(inlabru)

# set seed for reproducibility
set.seed(1)

# set GR law parameters
b_value <- 1 # b-value for GR-law
M0 <- 2.5 # magnitude of completeness
beta_value <- log(10)*b_value # beta parameter (transformation of b-value)

# simulate data
n_sim <- 100 # number of points
mags_values <- rexp(n = n_sim, rate = beta_value ) + M0 # magnitudes from exponential distribution

# plot GR law
mags_bins <- seq(M0, max(mags_values), by = 0.1) #set magnitude bins
mags_counts <- vapply(mags_bins, \(bin) sum(mags_values >= bin), 0) # calculate frequencies (number of points with magnitude greater or equal of bin)
line_value <- log10(n_sim) - b_value*(mags_bins - M0) # calculate theoretical values of counts


# visualise - I call this GR plot
ggplot() +
    geom_point(aes(x = mags_bins, y = log10(mags_counts))) +
    geom_line(aes(x = mags_bins, y = line_value))

# inlabru model fitting
# initialise dataframe
df_bru <- data.frame(mags = mags_values - M0)
# initialise component - on the left variable to be modelled, on the right linear model for log-rate (log(beta))
cmp <- mags ~ Intercept(1)
# fit the model
fit.bru <- bru(
    components = cmp, # components
    bru_obs( # model for observations
        mags ~ ., # set formula - this means use everything in the dataframe passed as data
        data = df_bru, # set dataframe of observations
        family = "exponential" # set model (exponential in this case)
    )
)

# look at summary information about model fitting
summary(fit.bru)

# retrieve posterior of beta parameter
# inla.tmarginal calculates the distribution of the transform of a distribution
beta_posterior = inla.tmarginal(exp, # function to be used
                                fit.bru$marginals.fixed$Intercept) # distribution to be trasnformed
# retrieve beta value maximum likelihood estimator
ML_estimator = 1/mean(mags_values - M0)

# set dataframe to plot vertical lines
vlines <- data.frame(
    xintercept = c(beta_value, ML_estimator),
    type = c("true", "ML")
)

# visualise posterior distribution of beta along with true value and ML estimator
ggplot(beta_posterior) +
    geom_line(aes(x,y)) +
    geom_vline(data = vlines, aes(xintercept = xintercept, color = type))

# retrieve posterior distribution of b-value
b_posterior = inla.tmarginal(\(x) exp(x)/log(10), fit.bru$marginals.fixed$Intercept)

# again dataframe for vertical lines
vlines2 <- vlines
vlines2$xintercept = vlines$xintercept/log(10)

# plot b-value distribution
ggplot(b_posterior) +
    geom_line(aes(x,y)) +
    geom_vline(data = vlines2, aes(xintercept = xintercept, color = type))


# calculate elements for GR plot
# GR plot using the posterior mean
b_value_post_mean = exp(fit.bru$summary.fixed$mean)/log(10)
line_value_post_mean <- log10(n_sim) - b_value_post_mean*(mags_bins - M0)

# GR plot using the posterior lower quantile
b_value_post_lower = inla.qmarginal(0.025,
                                    b_posterior)
line_value_post_lower <- log10(n_sim) - b_value_post_lower*(mags_bins - M0)

# GR plot using the posterior upper quantile
b_value_post_upper = inla.qmarginal(0.975,
                                    b_posterior)
line_value_post_upper <- log10(n_sim) - b_value_post_upper*(mags_bins - M0)
# GR plot using the maximum likelihood estimator
line_value_ml <- log10(n_sim) - (ML_estimator/log(10))*(mags_bins - M0)

# GR plot
pl.line1 <- ggplot() +
    geom_point(aes(x = mags_bins, y = log10(mags_counts))) +
    geom_line(aes(x = mags_bins, y = line_value, color = 'true')) +
    geom_line(aes(x = mags_bins, y = line_value_ml, color = 'ML')) +
    geom_line(aes(x = mags_bins, y = line_value_post_mean, color = 'posterior'), linetype = 2) +
    geom_line(aes(x = mags_bins, y = line_value_post_lower, color = 'posterior'), linetype = 2) +
    geom_line(aes(x = mags_bins, y = line_value_post_upper, color = 'posterior'), linetype = 2)
# visualise
pl.line1


###################################
# SAME EXAMPLE BUT WITH MORE DATA #
###################################

# sample the data
mags_values2 <- rexp(n = n_sim*10, rate = beta_value ) + M0
# calculates observed counts (used just for plotting)
mags_counts2 <- vapply(mags_bins, \(bin) sum(mags_values2 >= bin), 0)

# create dataframe for inlabru
df_bru2 <- data.frame(mags = mags_values2 - M0)
# fit model with inlabru
fit.bru2 <- bru(
    components = cmp,
    bru_obs(
        mags ~ .,
        data = df_bru2,
        family = "exponential"
    )
)

# extract posterior distribution
beta_posterior2 = inla.tmarginal(exp, fit.bru2$marginals.fixed$Intercept)
# calculate maximum likelihood estimator
ML_estimator2 = 1/mean(mags_values2 - M0)

# store posteriors in dataframes for plotting
beta_posterior = data.frame(beta_posterior)
beta_posterior2 = data.frame(beta_posterior2)

# create a column with the number of observations for color of the plot
beta_posterior$n = n_sim
beta_posterior2$n = n_sim*10
# create a column with the maximum likelihood estimate for vertical lines
beta_posterior$ML = ML_estimator
beta_posterior2$ML = ML_estimator2
# bind the dataframes by row (stack them together)
beta_posterior_bind = rbind(beta_posterior, beta_posterior2)

# beta posterior plot
ggplot(beta_posterior_bind) +
    geom_line(aes(x,y, color = factor(n))) +
    geom_vline(aes(xintercept = ML, color = factor(n))) +
    geom_vline(aes(xintercept = beta_value, color = 'true'))

# GR plot using posterior mean
b_value_post_mean2 = exp(fit.bru2$summary.fixed$mean)/log(10)
line_value_post_mean2 <- log10(n_sim*10) - b_value_post_mean2*(mags_bins - M0)

# GR plot using posterior lower quantile
b_value_post_lower2 = exp(fit.bru2$summary.fixed$`0.025quant`)/log(10)
line_value_post_lower2 <- log10(n_sim*10) - b_value_post_lower*(mags_bins - M0)

# GR plot using posterior upper quantile
b_value_post_upper2 = exp(fit.bru2$summary.fixed$`0.975quant`)/log(10)
line_value_post_upper2 <- log10(n_sim*10) - b_value_post_upper2*(mags_bins - M0)

# calculate true GR plot (the b-value is the same but log(n_sim) changed)
line_value2 <- log10(n_sim*10) - b_value*(mags_bins - M0)
# calculate  GR plot using maximum likelihood estimator
line_value2_ml <- log10(n_sim*10) - (ML_estimator2/log(10))*(mags_bins - M0)

# create GR plot
pl.line2 <- ggplot() +
    geom_point(aes(x = mags_bins, y = log10(mags_counts2))) +
    geom_line(aes(x = mags_bins, y = line_value2, color = 'true')) +
    geom_line(aes(x = mags_bins, y = line_value2_ml, color = 'ML')) +
    geom_line(aes(x = mags_bins, y = line_value_post_mean2, color = 'posterior - 1000'), linetype = 2) +
    geom_line(aes(x = mags_bins, y = line_value_post_lower2, color = 'posterior - 1000'), linetype = 2) +
    geom_line(aes(x = mags_bins, y = line_value_post_upper2, color = 'posterior - 1000'), linetype = 2)
# visualise
pl.line2

# visualise first and second GR plots together
(pl.line1|pl.line2)

#############
# REAL DATA #
#############

# import data
real_data = read.csv2(file = 'R/earthquakes.csv', sep = ',')
# isolate magnitudes
mag_values = real_data$magnitude

# creates magnitude bins for GR plot
mags_bins <- seq(min(mag_values), max(mag_values), by = 0.05) #set magnitude bins
# calculates magnitude counts for GR plot
mags_counts <- vapply(mags_bins, \(bin) sum(mag_values >= bin), 0) # calculate frequencies (number of points with magnitude greater or equal of bin)


# visualise
ggplot() +
    geom_point(aes(x = mags_bins, y = log10(mags_counts))) +
    xlim(1, 9) # set limits for x-axis
