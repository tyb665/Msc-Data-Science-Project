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

b_value_posterior = inla.tmarginal(\(x) x/log(10),
                                   beta_posterior)
plot (beta_posterior)
plot (b_value_posterior)

inla.qmarginal(c(0.025, 0.975), b_value_posterior)
inla.qmarginal(c(0.5), b_value_posterior)

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



##################################################
### Simulate data with different b-value models###
##################################################

library(ggplot2)
library(dplyr)

simulate_b_models <- function(n = 5000, M0 = 2.5, alpha = 0.05, lambda = 4000) {
  x <- seq(0, 1, length.out = n)
  i <- seq_len(n)
  
  b_models <- list(
    constant = rep(1.0, n),
    rect = ifelse(i < 0.4 * n, 0.98, ifelse(i <= 0.6 * n, 1.10, 0.98)),
    ramp = 1.05 - 0.1 * (i - 1) / (n - 1),
    sinusoidal = 1.0 + alpha * sin(2 * pi * i / lambda),
    gaussian = {
      set.seed(42)
      noise <- rnorm(n)
      # Simulate GP using low-pass filtering and smoothing
      smooth_b <- stats::filter(noise, rep(1/100, 100), sides = 2)  # moving average
      smooth_b <- as.numeric(scale(smooth_b)) * alpha + 1.0  # After standardization, add it to the mean value of 1
      # Fill the marginal NA
      smooth_b[is.na(smooth_b)] <- 1.0
      smooth_b
    }
  )
  
  # Generate earthquake magnitudes
  all_data <- lapply(names(b_models), function(model) {
    b <- b_models[[model]]
    beta <- log(10) * b
    mags <- rexp(n, rate = beta) + M0
    data.frame(x = x, i = i, b = b, beta = beta, mag = mags, model = model)
  }) %>% bind_rows()
  
  return(all_data)
}

# Simulation
sim_data <- simulate_b_models()

# Visualize b(x)
ggplot(sim_data, aes(x = x, y = b, color = model)) +
  geom_line() +
  labs(title = "b-value profiles", y = "b-value") +
  theme_minimal()

# Visualize the magnitude distribution
ggplot(sim_data, aes(x = mag, fill = model)) +
  geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
  labs(title = "Simulated magnitude distributions", x = "Magnitude") +
  theme_minimal()

  return(result)

#####################
###Inlabru fitting###
#####################

library(INLA)
library(inlabru)

fit_model <- function(data, M0 = 2.5, mesh_res = 100) {
  df <- data
  df$mag_shifted <- df$mag - M0
  mesh <- inla.mesh.1d(seq(0, 1, length.out = mesh_res), interval = c(0,1))
  
  spde <- inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(0.1, 0.5),
    prior.sigma = c(1, 0.5)  # Relax the prior restriction
  )
  
  cmp <- mag_shifted ~ Intercept(1) + beta_field(x, model = spde)
  
  fit <- bru(
    components = cmp,
    data = df,
    family = "exponential"
  )
  
  return(list(fit = fit, mesh = mesh, df = df))
}


## Use the fitting function to fit all models separately
library(dplyr)
models <- unique(sim_data$model)

fit_results <- lapply(models, function(m) {
  data_model <- dplyr::filter(sim_data, model == m)
  result <- fit_model(data_model)
  result$model <- m
  result
})
names(fit_results) <- models



## visual fitting result and the real b(x)
library(tidyr)
library(ggplot2)

plot_fitted_b <- function(result) {
  beta_mean <- result$fit$summary.random$beta_field$mean
  x_mesh <- result$mesh$loc
  b_est <- exp(beta_mean + result$fit$summary.fixed$mean[1]) / log(10)
  
  df_true <- result$df
  b_true_interp <- approx(x = df_true$x, y = df_true$b, xout = x_mesh, rule = 2)$y
  
  df_plot <- data.frame(x = x_mesh, Posterior = b_est, True = b_true_interp)
  
  df_long <- tidyr::pivot_longer(df_plot, cols = c("Posterior", "True"),
                                 names_to = "Source", values_to = "b_value")
  
  ggplot(df_long, aes(x = x, y = b_value, color = Source, linetype = Source)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("Posterior" = "blue", "True" = "red")) +
    scale_linetype_manual(values = c("Posterior" = "dashed", "True" = "solid")) +
    ylim(0.85, 1.15) +  # fixed the range of y-axis
    labs(title = paste("Model:", result$model), y = "b-value", x = "x") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "bottom"
    )
}


# Plot the fitting diagrams of each model one by one
plot_list <- lapply(fit_results, plot_fitted_b)
library(patchwork)
wrap_plots(plot_list, ncol = 2)



##Calculatie error metrics (MSE, RMSE)
calc_rmse <- function(result) {
  # Obtain intercept
  intercept <- result$fit$summary.fixed$mean[1]
  
  # Extract mesh Posterior value
  beta_mean <- result$fit$summary.random$beta_field$mean
  x_mesh <- result$mesh$loc
  b_est <- exp(beta_mean + intercept) / log(10)
  
  # Interpolate the real b value (with sorting)
  df_true <- result$df
  sorted_df <- df_true[order(df_true$x), ]
  b_true_interp <- approx(x = sorted_df$x, y = sorted_df$b, xout = x_mesh, rule = 2)$y
  
  valid <- !is.na(b_est) & !is.na(b_true_interp)
  if (sum(valid) == 0) {
    warning("No valid data points found for model: ", result$model)
    return(data.frame(model = result$model, RMSE = NA, MSE = NA))
  }
  
  rmse <- sqrt(mean((b_est[valid] - b_true_interp[valid])^2))
  mse <- mean((b_est[valid] - b_true_interp[valid])^2)
  
  data.frame(model = result$model, RMSE = rmse, MSE = mse)
}

rmse_table <- do.call(rbind, lapply(fit_results, calc_rmse))
print(rmse_table)

## Plot the histogram of RMSE and MSE
library(tidyr)
library(dplyr)

rmse_long <- rmse_table %>%
  pivot_longer(cols = c("RMSE", "MSE"),
               names_to = "Metric",
               values_to = "Value")
ggplot(rmse_long, aes(x = model, y = Value, fill = Metric)) +
  geom_col(position = "dodge") +
  labs(
    title = "RMSE and MSE for each b(x) model",
    x = "Model", y = "Error Value"
  ) + 
  scale_fill_manual(values = c("RMSE" = "#1f77b4", "MSE" = "#ff7f0e")) +
  theme_minimal()








