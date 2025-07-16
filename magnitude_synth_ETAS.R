## here we simulate observations according to the Epidemic Type Aftershock Sequence model that allows to account for temporal clustering

#' for the simulation we use the same parameters used in "Bayesian modeling of the temporal evolution of seismicity using the ETAS.inlabru package" 
#' (https://www.frontiersin.org/journals/applied-mathematics-and-statistics/articles/10.3389/fams.2023.1126759/full) in Section 3.1
#' 
#' This requires the installation of the ETAS.inlabru (https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/) R-package that can be done running the command
#'  remotes::install_github("edinburgh-seismicity-hub/ETAS.inlabru")
#'  

library(ETAS.inlabru)
library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)

# Parameters we use to generate synthetics, which we will refer to as the 'true' parameters
mu <- 0.15
K <- 0.1
alpha <- 2
c <- 0.11
p <- 1.1

# Format the true ETAS parameters for code to generate the synthetics
theta_etas <- data.frame(mu=mu, K=K, alpha=alpha, c=c, p=p)

# set right extreme of time interval (it starts at 0)
T2 <- 1500
# minimum magnitude
M0 = 2.5
# history - essenitally here we are imposing an event with magnitude 6.7 at time 1000 - you could try different options
Ht <- data.frame(ts=c(100), magnitudes=c(6.7))   # Impose a M6.7 event on day 1000 (or HT=null)


# set seed and simulate sequence
set.seed(1)
samp.etas.list <- generate_temporal_ETAS_synthetic(theta = theta_etas, #ETAS parameters 
                                                   beta.p = log(10),# beta-value (b-value*log(10)) for magnitude distribution
                                                   M0 = M0, # minimum magnitude
                                                   T1 = 0, # left extreme time interval
                                                   T2 = T2, # right extreme time interval
                                                   Ht=Ht) # history of the process - set to NULL if no history provided.
# merge in a dataframe representing the sample from ETAS model
samp.etas <- do.call(rbind, samp.etas.list)  

# time vs magnitude plot
ggplot() + 
  geom_point(aes(x = samp.etas$ts, y = samp.etas$magnitudes))

# histogram of number of observations by time
ggplot() + 
  geom_histogram(aes(x = samp.etas$ts), bins = 50)


# total number of points
n.samp = nrow(samp.etas)

# calculate magnitude counts for GR plot
mag_bins = seq(M0, (max(samp.etas$magnitudes) + 1), by = 0.1) 
mag_freqs = vapply(mag_bins, function(bin) sum(samp.etas$magnitudes >= bin), 0)

# GR plot
ggplot() + 
  geom_point(aes(x = mag_bins, y = log10(mag_freqs))) + 
  geom_line(aes(x = mag_bins, y = log10(n.samp) - 1*(mag_bins - M0)))


# set data.frame for inlabru model fitting
df_bru = data.frame(time = samp.etas$ts, 
                    mags = samp.etas$magnitudes - M0)


# set mesh points - you may want to try to extend it more than just 10
mesh_points <- seq(-10, T2 + 10, by = 20) # this sets mesh points - try others if you like
mesh1D <- fm_mesh_1d(mesh_points, degree = 2, boundary = "free")
# plot mesh
ggplot() +
  gg(mesh1D)

# create SPDE obj
# TRY DIFFERENT VALUES FOR THE RANGE PRIOR
the_spde <- inla.spde2.pcmatern(mesh1D,
                                prior.range = c(1, 0.01),
                                prior.sigma = c(1, 0.01))

# components of the model
comp <- mags ~ field(time, # this is the name of the time column in df_bru 
                     model = the_spde) + Intercept(1)

#4. fit the model
fit.bru <- bru(
  components = comp, # components
  bru_obs( # model for observations
    mags ~ field + Intercept, # set formula - this means use everything in the dataframe passed as data
    data = df_bru, # set dataframe of observations
    family = "exponential" # set model (exponential in this case)
  )
)

# look at summary informations
summary(fit.bru)

# data.frame of times for prediction
time4pred <- data.frame(time = seq(-100, T2 + 100, by = 20)) # same name used for time in df_bru

# calcualte prediction from posterior samples
pred.bru <- predict(fit.bru, # model
                    time4pred, # times at which predict
                    x ~ exp(field + Intercept)/log(10), # function to predict (b-value)
                    n.samples = 1000 #number of posterior samples
)

# plot predictions along with observations - this should be done as a double-axes plot with magnitude on the right y-axis and b-value in the left
ggplot() +
  gg(pred.bru) + 
  geom_point(aes(x = samp.etas$ts, y = samp.etas$magnitudes)) + 
  xlim(0,T2) + 
  ylim(0, 8)


##########################################################################
## CASE 2 - B-VALUE VARY WITH TIME IN CORRESPONDANCE OF BIG EARTHQUAKE ##
#########################################################################

#' This case is similar to the ones reported in Gulia & Wiemer paper on b-value temporal variations: the b-value is lower prior to a large earthquake
#' and increase after as the stress is released and large earthquakes become less likely

#' for this example we are going to use the same observarions (ETAS sample) used before and just change the magnitude
#' 

# initialise b_values (one for each observations) if observations before the large earthquake b-value = 0.8, and b-value = 1.2 if after
b_values <- c()
# populate vector with for loop
for(t in samp.etas$ts){
  if(t <= Ht$ts){
    b_values <- c(b_values, 0.8)
  }
  else{
    b_values <- c(b_values, 1.2)
  }
}

# plot b-value function - you may try different functions here
ggplot() + 
  geom_line(aes(x = samp.etas$ts, y = b_values))

# simulate new magnitude vector - for each observation we need to simulate a new magnitude

set.seed(2)
new_mags <- c() # initialise
for(idx in 1:length(samp.etas$ts)){ # populate vector with for loop
  curr_mag <- samp.etas$magnitudes[idx]
  b_val <- b_values[idx]
  if(curr_mag != 6.7){ # Avoid modifying that manually inserted major shock
    new_mag <- rexp(1, b_val*log(10)) + M0
    new_mags <- c(new_mags, new_mag)
  }
  else{
    new_mags <- c(new_mags, curr_mag)
  }
}


# time vs magnitude plots comparison
pl.mags_new <- ggplot() + 
  geom_point(aes(x = samp.etas$ts, y = new_mags)) + 
  labs(title = 'new magnitudes')

pl.mags_cont <- ggplot() + 
  geom_point(aes(x = samp.etas$ts, y = samp.etas$magnitudes)) + 
  labs(title = 'original magnitudes')

(pl.mags_cont | pl.mags_new)

  
  
# set data.frame for inlabru model fitting  
df_bru2 = data.frame(time = samp.etas$ts, 
                     mags = new_mags - M0)

# fit the model - we do not need to create again the components, the mesh, and the spde. 
# for now we try the same settings as before (we would like settings that works both for constant b-value and varying b-values)
# we are going to try a new mesh after

fit.bru2 <- bru(
  components = comp, 
  bru_obs( 
    mags ~ field + Intercept, 
    data = df_bru2, # new dataframe of observations
    family = "exponential" 
  )
)

summary(fit.bru2) # check how they changed wrt to the ones before

pred.bru2 <- predict(fit.bru2, # model
                    time4pred, # times at which predict
                    x ~ exp(field + Intercept)/log(10), # function to predict (b-value)
                    n.samples = 1000 #number of posterior samples
)

# plot posterior estimates of b-value function (not bad)
ggplot() +
  gg(pred.bru2) +
  geom_line(aes(x = samp.etas$ts, y = b_values)) + 
  xlim(0,T2) + 
  ylim(0.25, 1.75)


# here I try a different method to build the mesh - essentially the mesh points are equidistant before the large events - 
# while after the big events the distance increases 

# distance for equidistant part 
by_step = 20
# mesh points before big earthquake
mesh_points_pre <- seq(-10, Ht$ts, by = by_step)
# exponent - this regulates how fast the distance increases (you could try different options)
nn = 3
# mesh points after the big earthquake
mesh_points_after <- (seq(0^(1/nn), (T2 - Ht$ts + 10)^(1/nn), length.out = Ht$ts/by_step))^nn + Ht$ts 

# have a look - mesh points are closer at the beginning and then distance increases
plot(mesh_points_after, rep(1, length(mesh_points_after)))

# create mesh by combining mesh points into one vector
mesh1D_2 <- fm_mesh_1d(c(mesh_points_pre, mesh_points_after), degree = 2, boundary = "free")
# plot
ggplot() +
  gg(mesh1D_2)


# having a new mesh we need to create a new spde object and also the components
the_spde_2 <- inla.spde2.pcmatern(mesh1D_2,
                                prior.range = c(1, 0.01),
                                prior.sigma = c(1, 0.01))

#  components of the model
comp_2 <- mags ~ field(time, # this is the name of the time column in df_bru 
                     model = the_spde_2) + Intercept(1)

# fit the model
fit.bru3 <- bru(
  components = comp_2, # new components
  bru_obs( 
    mags ~ field + Intercept, 
    data = df_bru2, 
    family = "exponential" 
  )
)


pred.bru3 <- predict(fit.bru3, # model
                     time4pred, # times at which predict
                     x ~ exp(field + Intercept)/log(10), # function to predict (b-value)
                     n.samples = 1000 #number of posterior samples
)

pl3 <- ggplot() +
  gg(pred.bru3) +
  geom_line(aes(x = samp.etas$ts, y = b_values)) + 
  labs(title = 'equidistant mesh') +
  xlim(0,T2) + 
  ylim(0.25, 1.75)


pl2 <- ggplot() +
  gg(pred.bru2) +
  geom_line(aes(x = samp.etas$ts, y = b_values)) +
  labs(title = 'adaptive mesh') +
  xlim(0,T2) + 
  ylim(0.25, 1.75)


# compare predictions obtained with different meshes
(pl2 | pl3)

#########################
# REMARK ON THE RESULTS #
#########################

#' 1 . no much difference in terms of estimates - I think we could just consider the equidistant
#' 2. there is boundary effect at the start at the end - at the start overestimation and at the end underestimation
#'    - we can try to solve this by extending the mesh more (it was just 10 in the previous examples)
#' 3. it is varying too much 
#'    - we can try to solve this by considering a different prior for the range 
#'        


########################################################################
## REDUCING THE NUMBER OF MESH POINTS & ESTENDING THE BOUNDARIES MORE ##
########################################################################

# Now the mesh points are extended by 250 and we consider extremes to be a +-250
# Now the distance between mesh points is 50 before was 20
mesh_points_4 <- seq(-500, T2 + 500, by = 20) # this sets mesh points - try others if you like
mesh1D_4 <- fm_mesh_1d(mesh_points_4, degree = 2, boundary = "free")

ggplot() +
  gg(mesh1D_4)

# create spde
the_spde_4 <- inla.spde2.pcmatern(mesh1D_4,
                                prior.range = c(1, 0.01),
                                prior.sigma = c(1, 0.01))

# new components of the model with new spde
comp_4 <- mags ~ field(time,  
                     model = the_spde_4) + Intercept(1)

#4. fit the model
fit.bru_4 <- bru(
  components = comp_4, # new components
  bru_obs( 
    mags ~ field + Intercept, 
    data = df_bru2, 
    family = "exponential" 
  )
)


pred.bru4 <- predict(fit.bru_4, # model
                     time4pred, # times at which predict
                     x ~ exp(field + Intercept)/log(10), # function to predict (b-value)
                     n.samples = 1000 #number of posterior samples
)

pl4 <- ggplot() +
  gg(pred.bru4) +
  geom_line(aes(x = samp.etas$ts, y = b_values)) + 
  labs(title = 'extended mesh') + 
  xlim(0,1500) + 
  ylim(0.25, 1.75)


(pl4 | pl3 | pl2)

# no much difference here - so we can use the mesh with less mesh points (I would continue to extended to at least like 250 or 500)
# try different priors for the range and see how it goes :)



###1. Try different priors for the range and the standard deviation
library(dplyr)
library(tibble)

# ========================
# 1. Set search parameters
# ========================
prior_range_list <- round(c(T2/100, T2/50, T2/20, T2/10, T2/5, T2/2, T2/1.5), 1)
prior_sigma_list <- c(0.1, 0.3, 1, 3, 10)

# Initialize the result table
results <- tibble(
  range = numeric(),
  sigma = numeric(),
  DIC = numeric(),
  WAIC = numeric(),
  range_mean = numeric(),
  stdev_mean = numeric()
)

# ========================
# 2. Conduct grid search
# ========================
for (r in prior_range_list) {
  for (s in prior_sigma_list) {
    cat("Trying: prior.range =", r, "| prior.sigma =", s, "\n")
    
    # Create mesh
    mesh <- fm_mesh_1d(seq(-100, T2 + 100, by = 20), degree = 2, boundary = "free")
    
    # Create spde object
    spde <- inla.spde2.pcmatern(
      mesh,
      prior.range = c(r, 0.01),
      prior.sigma = c(s, 0.01)
    )
    
    # Define component
    comp_temp <- mags ~ field(time, model = spde) + Intercept(1)
    
    # Fit model
    fit <- tryCatch(
      bru(
        components = comp_temp,
        bru_obs(
          mags ~ field + Intercept,
          data = df_bru2,  # Use Case 2 data
          family = "exponential"
        ),
        options = list(
          control.compute = list(dic = TRUE, waic = TRUE)
        )
      ),
      error = function(e) {
        cat("Failed fitting", conditionMessage(e), "\n")
        return(NULL)
      }
    )
    
    # If the fitting is successful, extract the result
    if (!is.null(fit) && !is.null(fit$dic$dic)) {
      dic_val <- round(fit$dic$dic, 2)
      waic_val <- round(fit$waic$waic, 2)
      range_mean <- round(fit$summary.hyperpar["Range for field", "mean"], 2)
      stdev_mean <- round(fit$summary.hyperpar["Stdev for field", "mean"], 2)
      
      cat("Success| DIC:", dic_val, "| WAIC:", waic_val, "\n")
      
      results <- add_row(
        results,
        range = r,
        sigma = s,
        DIC = dic_val,
        WAIC = waic_val,
        range_mean = range_mean,
        stdev_mean = stdev_mean
      )
    } else {
      cat("The model fitting is completed but DIC/WAIC cannot be extracted. Skip it \n")
    }
  }
}

# ========================
# 3. Print results
# ========================
results_sorted <- results |> arrange(DIC)
print(results_sorted)

# Plot the result
results_heatmap <- results |>
  mutate(
    range = factor(range, levels = sort(unique(range))),
    sigma = factor(sigma, levels = sort(unique(sigma)))
  )

# Top plot:DIC heat map
p_dic <- ggplot(results_heatmap, aes(x = range, y = sigma, fill = DIC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(DIC, 1)), size = 3, color = "white") +
  scale_fill_viridis_c(option = "inferno", direction = -1) +
  labs(title = "DIC across Prior Settings",
       x = "Prior Range",
       y = "Prior Sigma",
       fill = "DIC") +
  theme_minimal()

# Bottom: WAIC heat map
p_waic <- ggplot(results_heatmap, aes(x = range, y = sigma, fill = WAIC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(WAIC, 1)), size = 3, color = "white") +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(title = "WAIC across Prior Settings",
       x = "Prior Range",
       y = "Prior Sigma",
       fill = "WAIC") +
  theme_minimal()

# Show both plot
p_dic / p_waic




###2.Try different strategies to set the mesh points

library(dplyr)
library(tibble)
library(inlabru)
library(INLA)
library(ggplot2)
library(patchwork)
library(glue)

# mesh strategies
mesh_strategies <- list(
  "Sparse" = seq(-100, T2 + 100, by = 50),
  "Dense" = seq(-100, T2 + 100, by = 10),
  "Equidistant" = seq(-10, T2 + 10, by = 20),
  "ExtendedBoundary" = seq(-500, T2 + 500, by = 20),
  "FocusOnBigQuake" = {
    big_quake_time <- Ht$ts
    mesh_before <- seq(-100, big_quake_time - 100, by = 40)
    mesh_mid <- seq(big_quake_time - 100, big_quake_time + 100, by = 10)
    mesh_after <- seq(big_quake_time + 100, T2 + 100, by = 40)
    sort(unique(c(mesh_before, mesh_mid, mesh_after)))
  }
)

# Initialize
plots <- list()
mesh_infos <- list()

# Data
data_bru <- df_bru2  

# Ergodic strategy
for (mesh_name in names(mesh_strategies)) {
  cat(glue::glue("\n Attempt the mesh stategy: {mesh_name} ...\n"))
  
  mesh_points <- mesh_strategies[[mesh_name]]
  mesh <- fm_mesh_1d(mesh_points, degree = 2, boundary = "free")
  spde <- inla.spde2.pcmatern(
    mesh,
    prior.range = c(100, 0.01),
    prior.sigma = c(0.5, 0.01)
  )
  comp <- mags ~ field(time, model = spde) + Intercept(1)
  
  # Safety fitting model
  fit <- tryCatch(
    bru(
      components = comp,
      bru_obs(mags ~ field + Intercept, data = data_bru, family = "exponential"),
      options = list(control.compute = list(dic = TRUE, waic = TRUE, hyperpar = TRUE))
    ),
    error = function(e) {
      cat(glue("Falied modeling {e$message}\n"))
      return(NULL)
    }
  )
  
  if (is.null(fit)) {
    mesh_infos[[mesh_name]] <- tibble(mesh = mesh_name, DIC = NA, WAIC = NA, range_mean = NA, stdev_mean = NA)
    plots[[mesh_name]] <- ggplot() + labs(title = paste("FAILED:", mesh_name))
    next
  }
  
  # Extract evaluation parameters (with checking)
  dic_val <- tryCatch(fit$dic$cpo$dic, error = function(e) NA)
  waic_val <- tryCatch(fit$waic$waic, error = function(e) NA)
  
  range_mean <- tryCatch(
    mean(inla.emarginal(identity, fit$marginals.hyperpar[["Range for field"]])),
    error = function(e) NA
  )
  stdev_mean <- tryCatch(
    mean(inla.emarginal(identity, fit$marginals.hyperpar[["Stdev for field"]])),
    error = function(e) NA
  )
  
  # Build tabular data (determine NA in advance)
  mesh_infos[[mesh_name]] <- tibble(
    mesh = mesh_name,
    DIC = ifelse(is.numeric(dic_val), round(dic_val, 2), NA),
    WAIC = ifelse(is.numeric(waic_val), round(waic_val, 2), NA),
    range_mean = round(range_mean, 2),
    stdev_mean = round(stdev_mean, 2)
  )
  
  # Prediction plot
  pred <- tryCatch(
    predict(fit, time4pred, ~ exp(field + Intercept) / log(10), n.samples = 1000),
    error = function(e) NULL
  )
  
  plots[[mesh_name]] <- if (!is.null(pred)) {
    ggplot() +
      gg(pred) +
      geom_line(aes(x = samp.etas$ts, y = b_values), color = "red") +
      labs(title = mesh_name) +
      xlim(0, T2) +
      ylim(0.25, 1.75)
  } else {
    ggplot() + labs(title = paste("PRED FAIL:", mesh_name))
  }
}

# Summary result
mesh_results <- bind_rows(mesh_infos) |> arrange(DIC)
print(mesh_results)

# Plot comparision
wrap_plots(plots, ncol = 2)


##3.try different simulated sequences.
library(inlabru)
library(INLA)
library(ggplot2)
library(patchwork)
library(dplyr)

# Simulation function: Supports different feature controls
simulate_sequence <- function(big_quake_time = 1000, 
                              n_events = 700,
                              clustering = "normal", 
                              T2 = 1500,
                              M0 = 2.5,
                              big_quake_mag = 6.7,
                              seed = 123) {
  set.seed(seed)
  ts <- sort(runif(n_events, 0, T2))
  
  # Simulated aggregation features (simple simulation, non-real ETAS)
  if (clustering == "high") {
    ts <- sort(c(ts, ts + runif(length(ts), 1, 5)))  # High Clustering: Add some aftershocks
  } else if (clustering == "low") {
    ts <- sort(runif(n_events, 0, T2))  # More uniform (approximately Poisson)
  }
  
  mags <- rexp(length(ts), rate = log(10)) + M0
  idx <- which.min(abs(ts - big_quake_time))
  ts[idx] <- big_quake_time
  mags[idx] <- big_quake_mag
  
  return(data.frame(ts = ts, magnitudes = mags))
}

# Model fitting function
fit_b_model <- function(df, mesh_points = seq(-500, 2000, by = 20),
                        prior_range = 100, prior_sigma = 0.5, M0 = 2.5) {
  df_bru <- data.frame(time = df$ts, mags = df$magnitudes - M0)
  
  mesh1D <- fm_mesh_1d(mesh_points, degree = 2, boundary = "free")
  the_spde <- inla.spde2.pcmatern(mesh1D,
                                  prior.range = c(prior_range, 0.01),
                                  prior.sigma = c(prior_sigma, 0.01))
  
  comp <- mags ~ field(time, model = the_spde) + Intercept(1)
  
  fit <- bru(
    components = comp,
    bru_obs(mags ~ field + Intercept, data = df_bru, family = "exponential"),
    options = list(control.compute = list(dic = TRUE, waic = TRUE, return.marginals = TRUE))
  )
  
  time4pred <- data.frame(time = seq(-100, max(df$ts) + 100, by = 20))
  pred <- predict(fit, time4pred, ~ exp(field + Intercept)/log(10), n.samples = 1000)
  
  return(pred)
}

# Different combinations of experimental Settings
experiments <- tibble::tibble(
  name = c("baseline", "early_quake", "late_quake", 
           "low_events", "high_events", "low_clustering", "high_clustering"),
  big_quake_time = c(1000, 500, 1300, 1000, 1000, 1000, 1000),
  n_events = c(700, 700, 700, 300, 1200, 700, 700),
  clustering = c("normal", "normal", "normal", "normal", "normal", "low", "high")
)

# conduct experiment
plots <- list()
for (i in seq_len(nrow(experiments))) {
  exp <- experiments[i, ]
  cat("Running scenario:", exp$name, "\n")
  
  df_sim <- simulate_sequence(
    big_quake_time = exp$big_quake_time,
    n_events = exp$n_events,
    clustering = exp$clustering
  )
  
  pred <- fit_b_model(df_sim)
  
  p <- ggplot() +
    gg(pred) +
    labs(title = paste("Scenario:", exp$name), x = "time", y = "b-value") +
    ylim(0.25, 1.75)
  
  plots[[exp$name]] <- p
}

# Show the comparison plots
patchwork::wrap_plots(plots, ncol = 2)


## 3.plus--double axis:
library(inlabru)
library(INLA)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scales)

# === Simulation function ===
simulate_sequence <- function(big_quake_time = 1000, 
                              n_events = 700,
                              clustering = "normal", 
                              T2 = 1500,
                              M0 = 2.5,
                              big_quake_mag = 6.7,
                              seed = 123) {
  set.seed(seed)
  ts <- sort(runif(n_events, 0, T2))
  
  if (clustering == "high") {
    ts <- sort(c(ts, ts + runif(length(ts), 1, 5)))
  } else if (clustering == "low") {
    ts <- sort(runif(n_events, 0, T2))
  }
  
  mags <- rexp(length(ts), rate = log(10)) + M0
  idx <- which.min(abs(ts - big_quake_time))
  ts[idx] <- big_quake_time
  mags[idx] <- big_quake_mag
  
  return(data.frame(ts = ts, magnitudes = mags))
}

# === Model fitting function ===
fit_b_model <- function(df, mesh_points = seq(-500, 2000, by = 20),
                        prior_range = 100, prior_sigma = 0.5, M0 = 2.5) {
  df_bru <- data.frame(time = df$ts, mags = df$magnitudes - M0)
  
  mesh1D <- fm_mesh_1d(mesh_points, degree = 2, boundary = "free")
  the_spde <- inla.spde2.pcmatern(mesh1D,
                                  prior.range = c(prior_range, 0.01),
                                  prior.sigma = c(prior_sigma, 0.01))
  
  comp <- mags ~ field(time, model = the_spde) + Intercept(1)
  
  fit <- bru(
    components = comp,
    bru_obs(mags ~ field + Intercept, data = df_bru, family = "exponential"),
    options = list(control.compute = list(dic = TRUE, waic = TRUE, return.marginals = TRUE))
  )
  
  time4pred <- data.frame(time = seq(-100, max(df$ts) + 100, by = 20))
  pred <- predict(fit, time4pred, ~ exp(field + Intercept)/log(10), n.samples = 1000)
  
  return(list(pred = pred, events = df))
}

# === Experimental setting ===
experiments <- tibble::tibble(
  name = c("baseline", "early_quake", "late_quake", 
           "low_events", "high_events", "low_clustering", "high_clustering"),
  big_quake_time = c(1000, 500, 1300, 1000, 1000, 1000, 1000),
  n_events = c(700, 700, 700, 300, 1200, 700, 700),
  clustering = c("normal", "normal", "normal", "normal", "normal", "low", "high")
)

# === Perform simulation, fitting and plot ===
plots <- list()
for (i in seq_len(nrow(experiments))) {
  exp <- experiments[i, ]
  cat("Running scenario:", exp$name, "\n")
  
  df_sim <- simulate_sequence(
    big_quake_time = exp$big_quake_time,
    n_events = exp$n_events,
    clustering = exp$clustering
  )
  
  result <- fit_b_model(df_sim)
  pred <- result$pred
  events <- result$events
  
  # Scale magnitude to the b-value range
  events$mag_rescaled <- rescale(events$magnitudes, to = c(0.25, 1.75))
  
  # Extraction of the main shock (magnitude ≥ 6.5)
  mainshock_df <- subset(events, magnitudes >= 6.5)
  
  # Generate a two-coordinate graph
  p <- ggplot() +
    geom_ribbon(data = pred, aes(x = time, ymin = q0.025, ymax = q0.975), fill = "grey80") +
    geom_line(data = pred, aes(x = time, y = mean), color = "black", size = 1) +
    
    # The small orange dot shows all events
    geom_point(data = events, aes(x = ts, y = mag_rescaled), 
               color = "orange", alpha = 0.5, size = 0.5) +
    
    # The main shock is marked with a big red dot
    geom_point(data = mainshock_df, aes(x = ts, y = mag_rescaled),
               color = "red", size = 2) +
    
    scale_y_continuous(
      name = "b-value", 
      limits = c(0.25, 1.75),
      sec.axis = sec_axis(~ rescale(., to = c(2.5, 6.7)), name = "Magnitude")
    ) +
    
    labs(title = paste("Scenario:", exp$name), x = "Time") +
    theme_minimal() +
    theme(
      axis.title.y.left = element_text(color = "black", size = 11),
      axis.title.y.right = element_text(color = "red", size = 11),
      axis.text.y.right = element_text(color = "red")
    )
  
  plots[[exp$name]] <- p
}

# === Show all plots ===
patchwork::wrap_plots(plots, ncol = 2)

