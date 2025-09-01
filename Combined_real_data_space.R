library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(maps)

# ==== read data ====
file_path <- "D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv"

eq_data <- read_csv(file_path) %>%
  rename(mag = magnitude) %>%
  drop_na(time, mag, latitude, longitude) %>%
  mutate(
    time = as.POSIXct(time, tz = "UTC"),
    year = year(time)
  ) %>%
  filter(year >= 1990)

# ==== select magnitude ≥ m0 ====
m0 <- 2.5
eq_filtered <- eq_data %>%
  filter(mag >= m0)

# ==== Standardized spatial coordinates (for SPDE)====
eq_filtered <- eq_filtered %>%
  mutate(
    x = scale(longitude)[, 1],
    y = scale(latitude)[, 1],
    mags = mag - m0
  )

# ==== California Border data ====
states <- map_data("state")
california <- states %>% filter(region == "california")

# ==== Plot ====
ggplot() +
  # 
  geom_polygon(data = california,
               aes(x = long, y = lat, group = group, linetype = "California boundary"),
               fill = NA, color = "black", size = 0.6) +
  
  # shock points
  geom_point(data = eq_filtered,
             aes(x = longitude, y = latitude),
             alpha = 0.3, size = 0.6) +
  
  coord_fixed() +
  scale_linetype_manual(values = c("California boundary" = "solid"),
                        name = "",  
                        guide = guide_legend(override.aes = list(color = "black"))) +
  geom_point(data = eq_filtered, aes(x = longitude, y = latitude, color = mag), alpha = 0, size = 0) +

  
  # Color legend
  scale_color_gradient(
    low = "gray80", high = "black",
    name = "Magnitude"
  ) +
  
  labs(title = paste("Earthquakes ≥", m0, "since 1990"),
       x = "Longitude", y = "Latitude") +
  theme_minimal()


##########################################################################################
##Build spatial mesh grids (for SPDE models)#############################################
##########################################################################################
library(INLA)
library(inlabru)

# ==== Extract spatial coordinates ====
eq_filtered <- eq_filtered %>%
  mutate(
    x = scale(longitude)[, 1],
    y = scale(latitude)[, 1],
    mags = mag - m0
  )

coords <- dplyr::select(eq_filtered, x, y)

# ==== KDE Estimate spatial density ====
set.seed(42)
kde <- MASS::kde2d(coords$x, coords$y, n = 200)

# ==== Probabilistic sampling points are conducted based on the KDE results ====
# transfer to vector
x_vec <- rep(kde$x, each = length(kde$y))
y_vec <- rep(kde$y, times = length(kde$x))
z_vec <- as.vector(kde$z)

# Remove extremely low-density areas (avoid remote points)
threshold <- quantile(z_vec, 0.05)
valid_idx <- which(z_vec > threshold)

sample_pool <- data.frame(x = x_vec[valid_idx], y = y_vec[valid_idx], weight = z_vec[valid_idx])
sample_pool <- sample_pool[!is.na(sample_pool$weight), ]

# Sampling mesh points (it is recommended to be within 1000 points
mesh_points <- sample_n(sample_pool, size = 600, weight = weight, replace = FALSE)

# Add boundary points (range control)
boundary_x <- range(coords$x)
boundary_y <- range(coords$y)
boundary_buffer <- 0.2

boundary <- inla.nonconvex.hull(as.matrix(coords), convex = -0.05)  # 非凸边界更符合加州形状

# ==== Construct a triangular grid ====
mesh <- inla.mesh.2d(
  loc = mesh_points[, c("x", "y")],
  boundary = boundary,
  max.edge = c(0.3, 1),    # Control the length of the triangular sides (More precise)
  cutoff = 0.05            # Limit the minimum distance between points to avoid excessive density
)

# ==== Visualize Mesh ====
plot(mesh, asp = 1, main = "Spatial Mesh (KDE-based)")
points(coords$x, coords$y, col = rgb(0,0,0,0.1), pch = 16, cex = 0.3)


################################################################################
## Spatial b-value modeling: KDE-based mesh + QC + inlabru fit with CPO/PIT  ##
################################################################################

## ---- Libraries ----
suppressPackageStartupMessages({
  library(dplyr)
  library(INLA)
  library(inlabru)
  library(MASS)
  library(sp)
})

## ---- Reproducibility ----
set.seed(42)

## ---- Inputs: eq_filtered is assumed ready with columns longitude, latitude, mag ----
## Example (if you need a filter): eq_filtered <- eq_data %>% filter(mag >= m0)
## Here we assume m0 already defined in your workspace
if (!exists("eq_filtered")) stop("eq_filtered not found. Please prepare eq_filtered first.")
if (!exists("m0")) stop("m0 not found. Please define m0 (magnitude of completeness).")

## ---- Coordinates & pre-processing ----
eq_filtered <- eq_filtered %>%
  mutate(
    x    = scale(longitude)[, 1],
    y    = scale(latitude)[, 1],
    mags = mag - m0
  )

coords <- dplyr::select(eq_filtered, x, y)

## ---- KDE estimate of spatial density (for sampling mesh points) ----
kde <- MASS::kde2d(coords$x, coords$y, n = 200)

# Vectorize KDE grid
x_vec <- rep(kde$x, each = length(kde$y))
y_vec <- rep(kde$y, times = length(kde$x))
z_vec <- as.vector(kde$z)

# Remove very low-density pixels to avoid remote points
threshold  <- quantile(z_vec, 0.05, na.rm = TRUE)
valid_idx  <- which(z_vec > threshold)
sample_pool <- data.frame(
  x = x_vec[valid_idx],
  y = y_vec[valid_idx],
  weight = z_vec[valid_idx]
) %>% filter(!is.na(weight))

## ---- Sample candidate mesh points from KDE weight ----
mesh_points <- dplyr::sample_n(sample_pool, size = 600, weight = weight, replace = FALSE)

## ---- Boundary (non-convex hull fits CA-like coastlines better) ----
boundary <- inla.nonconvex.hull(as.matrix(coords), convex = -0.05)

## ---- Build 2D triangular mesh ----
mesh <- inla.mesh.2d(
  loc     = as.matrix(mesh_points[, c("x", "y")]),
  boundary = boundary,
  max.edge = c(0.3, 1.0),  # smaller -> finer mesh near data, adjust as needed
  cutoff   = 0.05          # minimum distance to avoid overly dense duplicates
)

## ---- (Optional) Visualize Mesh + points ----
plot(mesh, asp = 1, main = "Spatial Mesh (KDE-based)")
points(coords$x, coords$y, col = rgb(0,0,0,0.1), pch = 16, cex = 0.3)

## =============================================================================
##                           Quick mesh quality checks
## =============================================================================
op <- par(mfrow = c(1,2))
plot(mesh, asp = 1, main = "KDE-based Mesh (with buffer)")
points(coords$x, coords$y, col = rgb(0,0,0,0.1), pch = 16, cex = 0.3)

ta <- INLA:::inla.mesh.fem(mesh)$ta   # triangle areas
hist(ta, breaks = 30, main = "Triangle area distribution", xlab = "Area")
par(op)

cat(sprintf("\nMesh QC: triangles=%d, area range=[%.4f, %.4f], median=%.4f\n",
            length(ta), min(ta), max(ta), median(ta)))

## =============================================================================
##                           SPDE model (PC priors)
## =============================================================================
## NOTE: Set priors to your calibrated values. The following are placeholders.
## prior.range = c(r0, pr0) encodes P(range < r0) = pr0
## prior.sigma = c(s0, ps0) encodes P(sigma > s0) = ps0
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(0.5, 0.5),   # <-- replace with your tuned prior
  prior.sigma = c(1.0, 0.01)   # <-- replace with your tuned prior
)

## =============================================================================
##                           inlabru components & data
## =============================================================================
## In inlabru v2, a common pattern is to feed coordinates via SpatialPoints or
## via a mapper. Here we use SpatialPoints for clarity.
bru_data <- data.frame(
  x    = coords$x,
  y    = coords$y,
  mags = eq_filtered$mags
)
sp::coordinates(bru_data) <- ~ x + y

## Component: spatial random field + intercept
## The name 'spatial' must match the term used in formula on the RHS.
cmp <- ~ Intercept(1) + spatial(coordinates, model = spde)

## =============================================================================
##                           Fit model with CPO/PIT
## =============================================================================
## IMPORTANT:
##  - Replace 'gaussian' and the formula on the next line to match your b-value
##    observation model. For example, if modeling magnitudes directly with
##    Gaussian errors: mags ~ Intercept + spatial
##  - If you use a different likelihood (e.g., Poisson for counts/intensity),
##    change 'family' and formula accordingly.
fit <- bru(
  components = cmp,
  family     = "gaussian",                     # <-- replace to your likelihood
  formula    = mags ~ Intercept + spatial,     # <-- replace to your model formula
  data       = bru_data,
  options    = list(
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
    control.inla    = list(int.strategy = "eb")
  )
)

## ---- Basic summaries ----
cat("\nModel fit summary (fixed effects):\n")
print(fit$summary.fixed)

if (!is.null(fit$summary.hyperpar)) {
  cat("\nHyperparameters summary:\n")
  print(fit$summary.hyperpar)
}

if (!is.null(fit$dic$dic))  cat(sprintf("\nDIC:  %.3f\n", fit$dic$dic))
if (!is.null(fit$waic$waic)) cat(sprintf("WAIC: %.3f\n", fit$waic$waic))

## =============================================================================
##                           CPO / PIT diagnostics
## =============================================================================
if (!is.null(fit$cpo)) {
  # CPO failures
  if (!is.null(fit$cpo$failure)) {
    cpo_fail <- sum(fit$cpo$failure > 0, na.rm = TRUE)
    cpo_rate <- cpo_fail / length(fit$cpo$failure)
    cat(sprintf("\nCPO failures: %d (rate = %.3f)\n", cpo_fail, cpo_rate))
  }
  
  # -log(CPO) summary (extreme values may flag bad fit/outliers)
  if (!is.null(fit$cpo$cpo)) {
    neglogcpo <- -log(fit$cpo$cpo)
    cat("\nSummary of -log(CPO):\n")
    print(summary(neglogcpo))
  }
  
  # PIT histogram ~ Uniform(0,1) if calibrated
  if (!is.null(fit$cpo$pit)) {
    hist(fit$cpo$pit, breaks = 20, main = "PIT histogram", xlab = "PIT")
  }
} else {
  warning("fit$cpo is NULL; CPO/PIT not computed. Check control.compute$cpo=TRUE.")
}

## =============================================================================
##                           Predict / Project (optional)
## =============================================================================
## If you want posterior of spatial field on mesh vertices:
## s_field <- predict(fit, spatial, ~ value)
## str(s_field)
##
## Or project to a grid:
## prj <- inla.mesh.projector(mesh, xlim = range(coords$x), ylim = range(coords$y), dims = c(100, 100))
## s_mean <- inla.mesh.project(prj, fit$summary.random$spatial$mean)
## image(list(x = prj$x, y = prj$y, z = s_mean), main = "Spatial field mean")

cat("\nDone.\n")


################################################################################
##### Search the optimal prior parameter combination with KDE ##########
################################################################################
library(readr)
library(dplyr)
library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)
library(tibble)
library(MASS)
library(tidyr)

# ==== 0. staderize coordinates ====
eq_filtered <- eq_filtered %>%
  mutate(
    mags = mag - m0,
    x = scale(longitude)[,1],
    y = scale(latitude)[,1]
  )

coords <- dplyr::select(eq_filtered, x, y)

# ==== 1. KDE estimated density + contour boundary ====
kde <- kde2d(eq_filtered$x, eq_filtered$y, n = 200)
contour_level <- quantile(kde$z, 0.3)  # Adjustable density threshold
contour_lines <- contourLines(kde$x, kde$y, kde$z, levels = contour_level)
boundary_coords <- cbind(contour_lines[[1]]$x, contour_lines[[1]]$y)

boundary <- inla.nonconvex.hull(boundary_coords, convex = -0.03)

# ==== 2. construct Mesh ====
# ==== 1. KDE density estimation ====
library(MASS)
kde <- kde2d(eq_filtered$x, eq_filtered$y, n = 200)

# ==== 2. Construct a KDE-based sampling point pool ====
x_vec <- rep(kde$x, each = length(kde$y))
y_vec <- rep(kde$y, times = length(kde$x))
z_vec <- as.vector(kde$z)

# Build a dataframe and remove low-density points
sample_pool <- data.frame(x = x_vec, y = y_vec, weight = z_vec) %>%
  filter(weight > quantile(weight, 0.05))  # adjustable

# Sampling mesh points (the greater the density, the higher the weight)
set.seed(42)
mesh_points <- sample_n(sample_pool, size = 800, weight = weight)

# ==== 3. Construct non-convex boundaries (the entire epicenter area) ====
boundary <- inla.nonconvex.hull(as.matrix(dplyr::select(eq_filtered, x, y)), convex = -0.03)

# ==== 4. construct mesh ====
mesh <- inla.mesh.2d(
  loc = mesh_points[, c("x", "y")],
  boundary = boundary,
  max.edge = c(0.2, 1),  # We can adjust later.the smaller, the finer
  cutoff = 0.05
)

# ==== 5. visual ====
#plot(mesh, main = "KDE-weighted Mesh")
#points(eq_filtered$x, eq_filtered$y, col = rgb(0, 0, 0, 0.1), pch = 16, cex = 0.3)

# ==== 3. set Prior grid ====
x_range <- diff(range(eq_filtered$x))
y_range <- diff(range(eq_filtered$y))
S <- round((x_range + y_range) / 2, 2)
cat("Estimated standardized space span S ≈", S, "\n")

prior_range_list <- round(c(S/50, S/20, S/10, S/5, S/2), 2)
prior_sigma_list <- c(0.1, 0.3, 1, 3, 5, 10)

# ==== 4. Construct data (including spatial coordinates) ====
df_bru_spatial <- eq_filtered %>%
  dplyr::select(mags, x, y) %>%
  mutate(coord = as.matrix(cbind(x, y)))

# ==== 5. initialize ====
results <- tibble(
  prior_range = numeric(),
  prior_sigma = numeric(),
  DIC = numeric(),
  WAIC = numeric(),
  range_mean = numeric(),
  stdev_mean = numeric()
)

# ==== 6. grid search ====
for (r in prior_range_list) {
  for (s in prior_sigma_list) {
    cat("Trying prior.range =", r, "| prior.sigma =", s, "\n")
    
    spde <- inla.spde2.pcmatern(
      mesh,
      prior.range = c(r, 0.01),
      prior.sigma = c(s, 0.01)
    )
    
    comp <- mags ~ field(coord, model = spde) + Intercept(1)
    
    fit <- tryCatch({
      bru(
        components = comp,
        data = df_bru_spatial,
        family = "exponential",
        options = list(
          control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),  # 加入 config
          control.inla = list(int.strategy = "eb")
        )
      )
    }, error = function(e) {
      message("Fit failed: ", e$message)
      return(NULL)
    })
    
    if (!is.null(fit)) {
      results <- add_row(
        results,
        prior_range = r,
        prior_sigma = s,
        DIC = round(fit$dic$dic, 2),
        WAIC = round(fit$waic$waic, 2),
        range_mean = round(fit$summary.hyperpar["Range for field", "mean"], 3),
        stdev_mean = round(fit$summary.hyperpar["Stdev for field", "mean"], 3)
      )
    }
  }
}

# ==== 7. DIC,WAIC Stdev for field (posterior mean) heatmap ====
results_plot <- results %>%
  mutate(
    prior_range = factor(prior_range, levels = sort(unique(prior_range))),
    prior_sigma = factor(prior_sigma, levels = sort(unique(prior_sigma)))
  )

best_dic_row <- results %>% filter(DIC == min(DIC, na.rm = TRUE))
best_waic_row <- results %>% filter(WAIC == min(WAIC, na.rm = TRUE))
# DIC 
p_dic <- ggplot(results_plot, aes(x = prior_range, y = prior_sigma, fill = DIC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(DIC, 1)), size = 3.2, color = "white") +
  geom_point(data = best_dic_row,
             aes(x = factor(prior_range), y = factor(prior_sigma)),
             shape = 8, color = "cyan", size = 3.5) +
  labs(title = "DIC across Prior Settings",
       x = "Prior Range", y = "Prior Sigma") +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  theme_minimal()
#WAIC hot map
p_waic <- ggplot(results_plot, aes(x = prior_range, y = prior_sigma, fill = WAIC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(WAIC, 1)), size = 3, color = "white") +
  geom_point(data = best_waic_row,
             aes(x = factor(prior_range), y = factor(prior_sigma)),
             color = "cyan", shape = 8, size = 3) +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(title = "WAIC across Prior Settings",
       x = "Prior Range", y = "Prior Sigma") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 0),
    axis.text.y = element_text(size = 10)
  )

# === Stdev for field (posterior mean) HOTMAP ===
best_stdev_row <- results %>% filter(stdev_mean == min(stdev_mean, na.rm = TRUE))

p_stdev <- ggplot(results_plot, aes(x = prior_range, y = prior_sigma, fill = stdev_mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(stdev_mean, 2)), size = 3.2, color = "white") +
  geom_point(data = best_stdev_row,
             aes(x = factor(prior_range), y = factor(prior_sigma)),
             shape = 8, color = "cyan", size = 3.5) +
  labs(title = "Stdev for Field (Posterior Mean)",
       x = "Prior Range", y = "Prior Sigma") +
  scale_fill_viridis_c(option = "inferno", direction = -1) +
  theme_minimal()

print((p_waic / p_dic / p_stdev) + plot_layout(guides = "collect"))


################################################################################
##### Full Pipeline: SPDE Spatial b-value Modeling with KDE Mesh ##############
################################################################################

library(readr)
library(dplyr)
library(INLA)
library(inlabru)
library(purrr)
library(MASS)
library(tibble)
library(ggplot2)

# ==== 0. Load data and standardize ====
file_path <- "D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv"
eq_data <- read_csv(file_path) %>%
  rename(mag = magnitude) %>%
  drop_na(mag, longitude, latitude) %>%
  filter(!is.na(mag) & !is.na(longitude) & !is.na(latitude))

m0 <- 2.5

eq_filtered <- eq_data %>%
  filter(mag >= m0) %>%
  mutate(
    mags = mag - m0,
    x = scale(longitude)[, 1],
    y = scale(latitude)[, 1]
  )

eq_filtered$coord <- cbind(eq_filtered$x, eq_filtered$y)


# ==== 1. Build KDE-weighted mesh ====
kde <- kde2d(eq_filtered$x, eq_filtered$y, n = 200)

x_vec <- rep(kde$x, each = length(kde$y))
y_vec <- rep(kde$y, times = length(kde$x))
z_vec <- as.vector(kde$z)

sample_pool <- data.frame(x = x_vec, y = y_vec, weight = z_vec) %>%
  filter(weight > quantile(weight, 0.05))

set.seed(2024)
mesh_points <- sample_n(sample_pool, size = 800, weight = weight)

boundary <- inla.nonconvex.hull(as.matrix(eq_filtered[, c("x", "y")]), convex = -0.03)

mesh <- inla.mesh.2d(
  loc = mesh_points[, c("x", "y")],
  boundary = boundary,
  max.edge = c(0.2, 1),
  cutoff = 0.05
)

# Optional: plot mesh
#plot(mesh, main = "KDE-weighted Mesh")
#points(eq_filtered$x, eq_filtered$y, col = rgb(0, 0, 0, 0.1), pch = 16, cex = 0.3)


# ==== 2. Build SPDE model ====
spde <- inla.spde2.pcmatern(
  mesh,
  prior.range = c(0.26, 0.01),
  prior.sigma = c(0.3, 0.01)
)

# ==== 3. Fit spatial b-value model ====
# construct coord（69961 x 2 matrix）
eq_filtered$coord <- cbind(eq_filtered$x, eq_filtered$y)

# Define the component: independent Intercept, explicitly named field
components <- ~ Intercept(1) + field_spatial(coord, model = spde)

# likelihood only references the latent component name and does not write the field name incorrectly
likelihood <- like(
  formula = mags ~ Intercept + field_spatial,
  family = "exponential",
  data = eq_filtered
)

# fit(No mean.linear")
fit <- bru(
  components = components,
  likelihood = likelihood,
  options = list(
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
    control.inla = list(int.strategy = "eb")
  )
)

# summary fit
summary(fit)

# ==== 5. Construct a predictive grid (KDE non-uniform sampling)====
# Sample the prediction points from the sample_pool according to the KDE weights
# Sampling prediction point
pred_pool <- sample_n(sample_pool, size = 5000, weight = weight)

# Construct a spatial matrix
coords_pred <- as.matrix(pred_pool[, c("x", "y")])

# prediction
pred <- predict(fit,
                formula = ~ exp(field_spatial + Intercept) / log(10),
                newdata = list(coord = coords_pred),
                n.samples = 1000)


# ==== 7. visual：b-value map====
library(ggplot2)
library(viridis)
library(maps)
library(dplyr)

# —— 1) back to standard longitude and latitude —— 
lon_mean <- mean(eq_filtered$longitude, na.rm = TRUE)
lon_sd   <- sd(eq_filtered$longitude, na.rm = TRUE)
lat_mean <- mean(eq_filtered$latitude,  na.rm = TRUE)
lat_sd   <- sd(eq_filtered$latitude,  na.rm = TRUE)

pred_df <- data.frame(
  x = coords_pred[, 1],
  y = coords_pred[, 2],
  b_mean  = pred$mean,
  b_sd    = pred$sd,
  b_lower = pred$q0.025,
  b_upper = pred$q0.975
) %>%
  mutate(
    longitude = x * lon_sd + lon_mean,
    latitude  = y * lat_sd + lat_mean,
    ci_width  = b_upper - b_lower
  )

# —— 2) California Boundary Data (Administrative Boundaries/Coastline)——
usa_map <- map_data("state")
california_map <- subset(usa_map, region == "california")

# —— 3) Plot: The posterior mean of value b —— 
p_mean <- ggplot() +
  geom_tile(data = pred_df,
            aes(x = longitude, y = latitude, fill = b_mean)) +
  geom_path(data = california_map,
            aes(x = long, y = lat, group = group, linetype = "California Border"),
            color = "black", linewidth = 0.5) +
  scale_fill_viridis_c(option = "plasma", name = "b-value") +
  coord_fixed() +
  labs(
    title = "Posterior Mean of Spatial b-value",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal()

# —— 4) Plot: Posterior standard deviation (Uncertainty)——
p_sd <- ggplot() +
  geom_tile(data = pred_df,
            aes(x = longitude, y = latitude, fill = b_sd)) +
  geom_path(data = california_map,
            aes(x = long, y = lat, group = group, linetype = "California Border"),
            color = "black", linewidth = 0.5) +
  scale_fill_viridis_c(option = "magma", name = "Posterior SD") +
  coord_fixed() +
  labs(
    title = "Posterior Uncertainty of b-value (SD)",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal()

# —— 5) Plot: Width of the 95% confidence interval —— 
p_ci <- ggplot() +
  geom_tile(data = pred_df,
            aes(x = longitude, y = latitude, fill = ci_width)) +
  geom_path(data = california_map,
            aes(x = long, y = lat, group = group, linetype = "California Border"),
            color = "black", linewidth = 0.5) +
  scale_fill_viridis_c(option = "inferno", name = "95% CI Width") +
  coord_fixed() +
  labs(
    title = "Credible Interval Width for b-value",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal()

# Print
print(p_mean)
print(p_sd)
print(p_ci)
########################################################################
########################################################################
########################################################################
#########Change to map prediction, not only data points######################################################
########################################################################
# =========================
# Spatial b-value (no covariates): fit → global grid prediction → robust maps
# =========================

# 0) Libraries
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(viridis); library(maps)
  library(INLA);  library(inlabru); library(MASS)
})

set.seed(2024)

# 1) Load & prep
file_path <- "D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv"
m0 <- 2.5

eq_data <- read_csv(file_path) |>
  rename(mag = magnitude) |>
  drop_na(mag, longitude, latitude) |>
  filter(mag >= m0)

eq_filtered <- eq_data |>
  mutate(
    mags = mag - m0,                    # response for Exp(rate=β)
    x = scale(longitude)[,1],
    y = scale(latitude)[,1]
  )
eq_filtered$coord <- cbind(eq_filtered$x, eq_filtered$y)

# 2) KDE-based mesh
kde <- MASS::kde2d(eq_filtered$x, eq_filtered$y, n = 200)
x_vec <- rep(kde$x, each = length(kde$y))
y_vec <- rep(kde$y, times = length(kde$x))
z_vec <- as.vector(kde$z)

sample_pool <- data.frame(x = x_vec, y = y_vec, weight = z_vec) |>
  filter(weight > quantile(weight, 0.05, na.rm = TRUE))

mesh_points <- dplyr::sample_n(sample_pool, size = 800, weight = weight, replace = FALSE)
boundary <- INLA::inla.nonconvex.hull(as.matrix(eq_filtered[,c("x","y")]), convex = -0.03)

mesh <- INLA::inla.mesh.2d(
  loc = mesh_points[,c("x","y")],
  boundary = boundary,
  max.edge = c(0.2, 1.0),
  cutoff = 0.05
)

# 3) SPDE + components  (η = log β)
spde <- INLA::inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(0.26, 0.01),
  prior.sigma = c(0.30, 0.01)
)

components <- ~ Intercept(1) + field_spatial(coord, model = spde)

likelihood <- inlabru::like(
  formula = mags ~ Intercept + field_spatial,  # mags ~ Exp(rate = β(s))
  family  = "exponential",
  data    = eq_filtered
)

fit <- inlabru::bru(
  components = components,
  likelihood = likelihood,
  options = list(
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
    control.inla    = list(int.strategy = "eb")
  )
)
print(summary(fit))

# 4) Global regular grid (standardized x,y)
grid_res <- 220
xr <- range(eq_filtered$x); yr <- range(eq_filtered$y)
coords_pred <- as.matrix(expand.grid(
  x = seq(xr[1], xr[2], length.out = grid_res),
  y = seq(yr[1], yr[2], length.out = grid_res)
))

# 5) Predict η = log β  → b = exp(η)/ln(10)
pred_eta <- predict(
  fit,
  ~ Intercept + field_spatial,
  newdata   = list(coord = coords_pred),
  n.samples = 2000
)
ln10 <- log(10)
pred_b <- list(
  mean   = exp(pred_eta$mean) / ln10,
  sd     = pred_eta$sd * exp(pred_eta$mean) / ln10,   # delta approximation, only for visible uncertainty
  q0.025 = exp(pred_eta$q0.025) / ln10,
  q0.975 = exp(pred_eta$q0.975) / ln10
)

# 6) Back to lon/lat & plotting df
lon_mean <- mean(eq_filtered$longitude); lon_sd <- sd(eq_filtered$longitude)
lat_mean <- mean(eq_filtered$latitude ); lat_sd <- sd(eq_filtered$latitude )

pred_df <- data.frame(
  x = coords_pred[,1], y = coords_pred[,2],
  b_mean  = pred_b$mean,
  b_sd    = pred_b$sd,
  b_lower = pred_b$q0.025,
  b_upper = pred_b$q0.975
) |>
  mutate(
    longitude = x * lon_sd + lon_mean,
    latitude  = y * lat_sd + lat_mean,
    ci_width  = b_upper - b_lower
  )

# 7) Sanity checks（See if there are any outliers pulling the color axis out）
cat("Empirical b ≈", round(1/(mean(eq_filtered$mags)*ln10), 3), "\n")
cat("Map mean b ≈", round(mean(pred_df$b_mean, na.rm=TRUE), 3), "\n")
print(summary(pred_df$b_mean))   # Check the minimum/maximum/quantile

# 8) California border overlay（Do not cut to avoid stripes）
usa_map <- map_data("state")
california_map <- subset(usa_map, region == "california")

# -- Core correction: The color axis only uses the robust quantiles of b_mean (not governed by individual extreme values) --
b_limits_mean <- as.numeric(quantile(pred_df$b_mean, probs = c(0.02, 0.98), na.rm = TRUE))


p_mean <- ggplot() +
  geom_raster(data = pred_df, aes(longitude, latitude, fill = b_mean)) +
  geom_path(data = california_map, aes(long, lat, group = group), color = "black", linewidth = 0.5) +
  scale_fill_viridis_c(option = "plasma", name = "b-value",
                       limits = b_limits_mean, oob = scales::squish) +
  coord_fixed() + labs(title = "Posterior Mean of b-value", x = "Longitude", y = "Latitude") +
  theme_minimal()

# Each uncertainty graph can use its own color axis
p_sd <- ggplot() +
  geom_raster(data = pred_df, aes(longitude, latitude, fill = b_sd)) +
  geom_path(data = california_map, aes(long, lat, group = group), color = "black", linewidth = 0.5) +
  scale_fill_viridis_c(option = "magma", name = "Posterior SD") +
  coord_fixed() + labs(title = "Posterior SD of b-value", x = "Longitude", y = "Latitude") +
  theme_minimal()

p_ci <- ggplot() +
  geom_raster(data = pred_df, aes(longitude, latitude, fill = ci_width)) +
  geom_path(data = california_map, aes(long, lat, group = group), color = "black", linewidth = 0.5) +
  scale_fill_viridis_c(option = "inferno", name = "95% CI Width") +
  coord_fixed() + labs(title = "95% CI Width of b-value", x = "Longitude", y = "Latitude") +
  theme_minimal()

print(p_mean); print(p_sd); print(p_ci)



################################################################################
#Consider the covariates depths########################################
################################################################################
# ==== 0. libraries ====
library(INLA)
library(inlabru)
library(dplyr)
library(tidyr)

# ==== 1. Read and preprocess data ====
eq_data <- read.csv("D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv")

eq_filtered <- eq_data %>%
  filter(!is.na(magnitude), magnitude >= 2.5,
         !is.na(longitude), !is.na(latitude), !is.na(depth)) %>%
  mutate(
    mags        = magnitude - 2.5,                           # The magnitude minus M0 is taken as the response variable
    x           = as.numeric(scale(longitude)[,1]),          # Longitude standardization
    y           = as.numeric(scale(latitude)[,1]),           # Latitude standardization
    depth_s     = as.numeric(scale(depth)[,1]),              # Depth standardization
    coordinates = cbind(x, y)                                # Coordinates matrix
  )

# ==== 2. Construct Mesh ====
library(MASS)

# Step 1: KDE density estimation
kde <- kde2d(eq_filtered$x, eq_filtered$y, n = 200)

# Step 2: Build a sampling pool (excluding low-density areas)
x_vec <- rep(kde$x, each = length(kde$y))
y_vec <- rep(kde$y, times = length(kde$x))
z_vec <- as.vector(kde$z)

sample_pool <- data.frame(x = x_vec, y = y_vec, weight = z_vec) %>%
  filter(weight > quantile(weight, 0.05))   # Adjustable threshold

# Step 3: Weighted sampling mesh points
set.seed(2024)
mesh_points <- sample_n(sample_pool, size = 800, weight = weight)

# Step 4: Non-convex boundaries (maintaining the outline of California)
boundary <- inla.nonconvex.hull(as.matrix(eq_filtered[, c("x", "y")]), convex = -0.03)

# Step 5: Construct mesh
mesh <- inla.mesh.2d(
  loc = mesh_points[, c("x", "y")],
  boundary = boundary,
  max.edge = c(0.2, 1),
  cutoff = 0.05
)


# ==== 3. Define two SPDE models ====
spde0 <- inla.spde2.pcmatern(mesh, prior.range = c(0.26, 0.01), prior.sigma = c(0.30, 0.01))
spde1 <- inla.spde2.pcmatern(mesh, prior.range = c(0.26, 0.01), prior.sigma = c(0.20, 0.01))

# Construct the projection matrix
A0 <- inla.spde.make.A(mesh, loc = eq_filtered$coordinates)
A1 <- inla.spde.make.A(mesh, loc = eq_filtered$coordinates, weights = eq_filtered$depth_s)


# ==== 4. Define model components (mapper automatically implements variable coefficients)====
components <- ~ Intercept(1) +
  depth_lin(depth_s, model = "linear") +
  field0(coordinates, model = spde0) +
  field1(coordinates, model = spde1)

# ==== 5. Define the likelihood function ====
likelihood <- like(
  formula = mags ~ Intercept + depth_lin + field0 + field1,
  family  = "exponential",
  data    = eq_filtered
)

# ==== 6. Fitting model ====
fit_vc <- bru(
  components = components,
  likelihood = likelihood,
  options = list(
    bru_initialisation = list(A = list(A0, A1)),  # Pass in the projection matrices of the two fields
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
    control.inla = list(int.strategy = "eb")
  )
)

# ==== 7. Summary of fitting ====
summary(fit_vc)

##Prediction and visualization
# ===============================
library(INLA)
library(inlabru)
library(dplyr)
library(ggplot2)
library(viridis)
library(maps)

# ---- 0) standardized longitude and latitude ----
lon_mean <- mean(eq_filtered$longitude, na.rm = TRUE)
lon_sd   <- sd(eq_filtered$longitude, na.rm = TRUE)
lat_mean <- mean(eq_filtered$latitude,  na.rm = TRUE)
lat_sd   <- sd(eq_filtered$latitude,  na.rm = TRUE)

lonlat_back <- function(df_xy) {
  df_xy %>%
    mutate(
      longitude = x * lon_sd + lon_mean,
      latitude  = y * lat_sd + lat_mean
    )
}

usa_map <- map_data("state")
california_map <- subset(usa_map, region == "california")

plot_raster <- function(df, fill_col, title_txt, fill_name, option = "plasma") {
  ggplot() +
    geom_tile(data = df, aes(x = longitude, y = latitude, fill = .data[[fill_col]])) +
    geom_path(data = california_map,
              aes(x = long, y = lat, group = group, linetype = "California Border"),
              color = "black", linewidth = 0.5) +
    scale_linetype_manual(name = "", values = c("California Border" = "solid")) +
    guides(linetype = guide_legend(override.aes = list(color = "black", linewidth = 0.8))) +
    coord_fixed() +
    scale_fill_viridis_c(option = option, name = fill_name) +
    labs(title = title_txt, x = "Longitude", y = "Latitude") +
    theme_minimal()
}

# ---- 1) Prediction point: preferred kde sample pool, otherwise the regular grid----
if (exists("sample_pool") && all(c("x","y") %in% names(sample_pool))) {
  set.seed(2024)
  pred_points <- dplyr::sample_n(sample_pool, size = 5000, weight = weight)
  coords_pred <- as.matrix(pred_points[, c("x","y")])
} else {
  grid_res <- 200
  xr <- range(eq_filtered$x); yr <- range(eq_filtered$y)
  gx <- seq(xr[1], xr[2], length.out = grid_res)
  gy <- seq(yr[1], yr[2], length.out = grid_res)
  grid <- expand.grid(x = gx, y = gy)
  coords_pred <- as.matrix(grid[, c("x","y")])
}

# ---- 2) Generate A matrix(Note that field1 requires a depth weight)----
make_pred_As <- function(coords_pred, depth_s_vec) {
  A0_pred <- INLA::inla.spde.make.A(mesh, loc = coords_pred)                         # field0
  A1_pred <- INLA::inla.spde.make.A(mesh, loc = coords_pred, weights = depth_s_vec) # field1 * depth
  list(A0_pred = A0_pred, A1_pred = A1_pred)
}

# ---- 3) Unified Prediction Function (Key modification: The variable name in newdata must be depth_s; the formula should be '+ depth_lin')----
predict_b <- function(depth_s_vec, nsamp = 1000) {
  stopifnot(length(depth_s_vec) == nrow(coords_pred))
  As <- make_pred_As(coords_pred, depth_s_vec)
  pred <- predict(
    fit_vc,
    # depth_lin will automatically multiply the coefficient with newdata$depth_s; The multiplication of field1 is accomplished by the weights of A1
    ~ exp(Intercept + depth_lin + field0 + field1) / log(10),
    newdata = list(coordinates = coords_pred, depth_s = depth_s_vec),
    A = list(As$A0_pred, As$A1_pred),
    n.samples = nsamp
  )
  data.frame(
    x = coords_pred[,1], y = coords_pred[,2],
    b_mean  = pred$mean,
    b_sd    = pred$sd,
    b_lower = pred$q0.025,
    b_upper = pred$q0.975
  )
}

# ---- 4) Scenario A: Representativeness Depth (median) ----
depth_rep <- median(eq_filtered$depth_s, na.rm = TRUE)
pred_rep  <- predict_b(rep(depth_rep, nrow(coords_pred)), nsamp = 1000) %>% lonlat_back()
pred_rep$ci_width <- pred_rep$b_upper - pred_rep$b_lower

p_rep_mean <- plot_raster(pred_rep, "b_mean",  "Posterior Mean of b (depth = median)", "b-value", "plasma")
p_rep_sd   <- plot_raster(pred_rep, "b_sd",    "Posterior SD of b (depth = median)",   "Posterior SD", "magma")
p_rep_ci   <- plot_raster(pred_rep, "ci_width","95% CI Width of b (depth = median)",   "95% CI Width", "inferno")

print(p_rep_mean); print(p_rep_sd); print(p_rep_ci)

# ---- 5) Scenario B: Shallow vs. deep contrast + Δb ----
q_shallow <- as.numeric(quantile(eq_filtered$depth_s, 0.25, na.rm = TRUE))
q_deep    <- as.numeric(quantile(eq_filtered$depth_s, 0.75, na.rm = TRUE))  # We can change as 0.90

pred_shallow <- predict_b(rep(q_shallow, nrow(coords_pred)), nsamp = 1000) %>% lonlat_back()
pred_deep    <- predict_b(rep(q_deep,    nrow(coords_pred)), nsamp = 1000) %>% lonlat_back()

p_shallow <- plot_raster(pred_shallow, "b_mean",
                         sprintf("Mean b at shallow depth (Q25)", q_shallow),
                         "b-value", "plasma")
p_deepmap <- plot_raster(pred_deep, "b_mean",
                         sprintf("Mean b at deep depth (Q75)", q_deep),
                         "b-value", "plasma")
print(p_shallow); print(p_deepmap)

# Δb：deep − shallow
delta <- left_join(
  pred_deep, pred_shallow,
  by = c("x", "y", "longitude", "latitude"),
  suffix = c("_deep", "_shallow")
) %>%
  mutate(
    db_mean = b_mean_deep - b_mean_shallow,
    # Significance judgment of confidence intervals
    sig_neg = (b_upper_deep < b_lower_shallow),
    sig_pos = (b_lower_deep > b_upper_shallow)
  )


p_db <- plot_raster(delta, "db_mean",
                    "Δb = b(deep) − b(shallow)", "Δb", "cividis")
print(p_db)

# ---- 6) Scenario C: Depth Effect Intensity (spatial variable coefficient) map----
# Let depth_s = 1, then depth_lin + field1 is the unit slope of log β versus depth_s
depth_one  <- rep(1, nrow(coords_pred))
As_slope   <- make_pred_As(coords_pred, depth_s_vec = depth_one)
pred_slope <- predict(
  fit_vc,
  ~ depth_lin + field1,
  newdata = list(coordinates = coords_pred, depth_s = depth_one),
  A = list(As_slope$A0_pred, As_slope$A1_pred),
  n.samples = 1000
)

slope_df <- data.frame(
  x = coords_pred[,1], y = coords_pred[,2],
  slope_mean  = pred_slope$mean,
  slope_sd    = pred_slope$sd,
  slope_lower = pred_slope$q0.025,
  slope_upper = pred_slope$q0.975
) %>% lonlat_back()

p_slope    <- plot_raster(slope_df, "slope_mean",
                          "Spatially-varying depth effect on log β (β_global + δβ(s))",
                          "d logβ / d depth_s", "viridis")
p_slope_sd <- plot_raster(slope_df, "slope_sd",
                          "Uncertainty of depth effect (SD)",
                          "Posterior SD", "magma")
print(p_slope); print(p_slope_sd)

############################################################
############################################################
############map################################################
########################################################################
########################################################################
########################################################################
# ===============================
# Prediction & Visualization  (for model with depth covariate)
# ===============================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(viridis); library(maps)
})

# ---- 0) lon/lat Anti-standardization ----
lon_mean <- mean(eq_filtered$longitude, na.rm = TRUE)
lon_sd   <- sd(eq_filtered$longitude,   na.rm = TRUE)
lat_mean <- mean(eq_filtered$latitude,  na.rm = TRUE)
lat_sd   <- sd(eq_filtered$latitude,    na.rm = TRUE)

lonlat_back <- function(df_xy) {
  df_xy %>%
    mutate(
      longitude = x * lon_sd + lon_mean,
      latitude  = y * lat_sd + lat_mean
    )
}

usa_map <- map_data("state")
california_map <- subset(usa_map, region == "california")

plot_raster <- function(df, fill_col, title_txt, fill_name, option = "plasma",
                        limits = NULL, midpoint = NULL) {
  p <- ggplot() +
    geom_raster(data = df, aes(x = longitude, y = latitude, fill = .data[[fill_col]])) +
    geom_path(data = california_map,
              aes(x = long, y = lat, group = group),
              color = "black", linewidth = 0.5) +
    coord_fixed() +
    labs(title = title_txt, x = "Longitude", y = "Latitude") +
    theme_minimal()
  if (is.null(midpoint)) {
    p + scale_fill_viridis_c(option = option, name = fill_name,
                             limits = limits, oob = scales::squish)
  } else {
    p + scale_fill_gradient2(name = fill_name, midpoint = midpoint,
                             low = "#3b4cc0", mid = "white", high = "#b40426",
                             limits = limits, oob = scales::squish)
  }
}

grid_res <- 220
xr <- range(eq_filtered$x); yr <- range(eq_filtered$y)
coords_pred <- as.matrix(expand.grid(
  x = seq(xr[1], xr[2], length.out = grid_res),
  y = seq(yr[1], yr[2], length.out = grid_res)
))

# ---- 2) Generate an A matrix (field0 normal, field1 with depth_s weights)----
make_pred_As <- function(coords_pred, depth_s_vec) {
  A0_pred <- INLA::inla.spde.make.A(mesh, loc = coords_pred)                          # field0
  A1_pred <- INLA::inla.spde.make.A(mesh, loc = coords_pred, weights = depth_s_vec)   # field1 * depth
  list(A0_pred = A0_pred, A1_pred = A1_pred)
}

# ---- 3) Unified prediction function: First calculate η=logβ, and then b=exp(η)/ln10----
predict_b <- function(depth_s_vec, nsamp = 1500) {
  stopifnot(length(depth_s_vec) == nrow(coords_pred))
  As <- make_pred_As(coords_pred, depth_s_vec)
  
  # η = Intercept + depth_lin + field0 + field1
  pred_eta <- predict(
    fit_vc,
    ~ Intercept + depth_lin + field0 + field1,
    newdata = list(coordinates = coords_pred, depth_s = depth_s_vec),
    A = list(As$A0_pred, As$A1_pred),     # 顺序与 formula 中 field0, field1 对齐
    n.samples = nsamp
  )
  
  ln10 <- log(10)
  data.frame(
    x = coords_pred[,1], y = coords_pred[,2],
    b_mean  = exp(pred_eta$mean)    / ln10,
    b_sd    = pred_eta$sd * exp(pred_eta$mean) / ln10,   # delta-method 近似
    b_lower = exp(pred_eta$q0.025)  / ln10,
    b_upper = exp(pred_eta$q0.975)  / ln10
  )
}

# ---- 4) Scenario A: Median depth ----
depth_rep <- median(eq_filtered$depth_s, na.rm = TRUE)
pred_rep  <- predict_b(rep(depth_rep, nrow(coords_pred)), nsamp = 1500) %>%
  lonlat_back() %>%
  mutate(ci_width = b_upper - b_lower)

# Robust color axis: Take the 2% to 98% quantiles
b_limits <- as.numeric(quantile(pred_rep$b_mean, c(0.02, 0.98), na.rm = TRUE))

p_rep_mean <- plot_raster(pred_rep, "b_mean",
                          sprintf("Posterior Mean of b (depth = median)"),
                          "b-value", "plasma", limits = b_limits)
p_rep_sd   <- plot_raster(pred_rep, "b_sd",
                          "Posterior SD of b (depth = median)",
                          "Posterior SD", "magma")
p_rep_ci   <- plot_raster(pred_rep, "ci_width",
                          "95% CI Width of b (depth = median)",
                          "95% CI Width", "inferno")
print(p_rep_mean); print(p_rep_sd); print(p_rep_ci)

# ---- 5) Scenario B: Shallow/Deep contrast + Δb ----
q_shallow <- as.numeric(quantile(eq_filtered$depth_s, 0.25, na.rm = TRUE))
q_deep    <- as.numeric(quantile(eq_filtered$depth_s, 0.75, na.rm = TRUE))

pred_shallow <- predict_b(rep(q_shallow, nrow(coords_pred))) %>% lonlat_back()
pred_deep    <- predict_b(rep(q_deep,    nrow(coords_pred))) %>% lonlat_back()

p_shallow <- plot_raster(pred_shallow, "b_mean",
                         sprintf("Mean b at shallow depth (Q25)"),
                         "b-value", "plasma", limits = b_limits)
p_deepmap <- plot_raster(pred_deep, "b_mean",
                         sprintf("Mean b at deep depth (Q75)"),
                         "b-value", "plasma", limits = b_limits)
print(p_shallow); print(p_deepmap)

# Δb（deep − shallow），Divergent color marking & significant marking
delta <- left_join(
  pred_deep, pred_shallow,
  by = c("x","y"),
  suffix = c("_deep","_shallow")
) %>%
  mutate(
    longitude = x * lon_sd + lon_mean,
    latitude  = y * lat_sd + lat_mean,
    db_mean = b_mean_deep - b_mean_shallow,
    sig_neg = (b_upper_deep < b_lower_shallow),
    sig_pos = (b_lower_deep > b_upper_shallow)
  )

db_lim <- max(abs(quantile(delta$db_mean, c(0.02, 0.98), na.rm = TRUE)))
p_db <- plot_raster(delta, "db_mean",
                    "Δb = b(deep) − b(shallow)", "Δb",
                    option = NULL, limits = c(-db_lim, db_lim), midpoint = 0)
print(p_db)

# ---- 6) Scenario C: Slope of depth effect（d logβ / d depth_s）----
depth_one <- rep(1, nrow(coords_pred))
A1_only   <- INLA::inla.spde.make.A(mesh, loc = coords_pred, weights = depth_one)

pred_slope <- predict(
  fit_vc,
  ~ depth_lin + field1,   #  d logβ / d depth_s
  newdata = list(coordinates = coords_pred, depth_s = depth_one),
  A = list(A1_only),      # Pass only the A of field1; depth_lin does not require A
  n.samples = 1500
)

slope_df <- data.frame(
  x = coords_pred[,1], y = coords_pred[,2],
  slope_mean  = pred_slope$mean,
  slope_sd    = pred_slope$sd,
  slope_lower = pred_slope$q0.025,
  slope_upper = pred_slope$q0.975
) %>% lonlat_back()

# Pass only the A of field1; depth_lin does not require A
slim <- max(abs(quantile(slope_df$slope_mean, c(0.02, 0.98), na.rm = TRUE)))
p_slope    <- plot_raster(slope_df, "slope_mean",
                          "Depth effect on log β (d logβ / d depth_s)",
                          "d logβ / d depth_s",
                          option = NULL, limits = c(-slim, slim), midpoint = 0)
p_slope_sd <- plot_raster(slope_df, "slope_sd",
                          "Uncertainty of depth effect (SD)",
                          "Posterior SD", "magma")
print(p_slope); print(p_slope_sd)

# ---- 7) Sanity check----
ln10 <- log(10)
b_empirical <- 1 / (mean(eq_filtered$mags) * ln10)
cat("Empirical b ≈", round(b_empirical, 3), "\n")
cat("Median-depth map mean b ≈", round(mean(pred_rep$b_mean, na.rm = TRUE), 3), "\n")
