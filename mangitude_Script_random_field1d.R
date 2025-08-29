# set b-values
n_obs = 999


b_values = c(rep(0.75, n_obs/3), rep(1.25, n_obs/3), rep(0.75, n_obs/3))

plot(b_values, type = 'l')
M0 <- 2.5 # magnitude of completeness

# simulate magnitudes

mag_values = c()
for( b in b_values){
  beta_value <- log(10)*b # beta parameter (transformation of b-value)
  mag <- rexp(n = 1, rate = beta_value ) + M0 # magnitudes from exponential distribution
  mag_values <- c(mag_values, mag) 
}

plot(mag_values)

time_values = 1:length(mag_values)

df_bru = data.frame(time = time_values, 
                    mags = mag_values - M0)

df_bru

library(INLA)
library(inlabru)
library(ggplot2)


# 1. create the mesh
# -- have to start before the starting time of the data
# -- have to end after the end of the data
# -- TRY DIFFERENT MESHES
mesh_points <- seq(-10, n_obs + 10, by = 20) # this sets mesh points - try others if you like
mesh1D <- fm_mesh_1d(mesh_points, degree = 2, boundary = "free")

ggplot() +
  gg(mesh1D)

# 2. create SPDE obj
# TRY DIFFERENT VALUES FOR THE RANGE PRIOR
the_spde <- inla.spde2.pcmatern(mesh1D,
                                prior.range = c(100, 0.01),
                                prior.sigma = c(1, 0.01))

# 3. components of the model
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

# 5. posterior on the parameters
summary(fit.bru)

#plot(fit.bru, "Intercept")
#plot(fit.bru, "Range for field")

# 6. posterior prediction for estimating quantities
# times of the prediction
time4pred <- data.frame(time = mesh_points) # same name used for time in df_bru

pred.bru <- predict(fit.bru, # model
                     time4pred, # times at which predict
                     x ~ exp(field + Intercept)/log(10), # function to predict (b-value)
                     n.samples = 1000 #number of posterior samples
)

ggplot() +
  gg(pred.bru) + 
  geom_line(aes(x = 1:n_obs, y = b_values))




