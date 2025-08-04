########## Space modeling
library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)

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

# ==== Visualization of epicentral distribution====
ggplot(eq_filtered, aes(x = longitude, y = latitude)) +
  geom_point(alpha = 0.3, size = 0.6) +
  coord_fixed() +
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
contour_level <- quantile(kde$z, 0.3)  # 可调密度阈值
contour_lines <- contourLines(kde$x, kde$y, kde$z, levels = contour_level)
boundary_coords <- cbind(contour_lines[[1]]$x, contour_lines[[1]]$y)

boundary <- inla.nonconvex.hull(boundary_coords, convex = -0.03)

# ==== 2. construct Mesh ====
# ==== 1. KDE density estimation ====
library(MASS)
kde <- kde2d(eq_filtered$x, eq_filtered$y, n = 200)

# ==== 2. 构建 KDE-based 采样点池 ====
x_vec <- rep(kde$x, each = length(kde$y))
y_vec <- rep(kde$y, times = length(kde$x))
z_vec <- as.vector(kde$z)

# Build a dataframe and remove low-density points
sample_pool <- data.frame(x = x_vec, y = y_vec, weight = z_vec) %>%
  filter(weight > quantile(weight, 0.05))  # 可调阈值

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
      message("❌ Fit failed: ", e$message)
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


# ==== 7. visual：b-value map ====
library(ggplot2)
library(viridis)

# coordinates information
pred_df <- data.frame(
  x = coords_pred[, 1],
  y = coords_pred[, 2],
  b_mean = pred$mean,
  b_sd = pred$sd,
  b_lower = pred$q0.025,
  b_upper = pred$q0.975
)

# b-value mean value map
ggplot(pred_df, aes(x = x, y = y, fill = b_mean)) +
  geom_tile() +
  coord_equal() +
  scale_fill_viridis_c(option = "plasma", name = "b-value") +
  labs(
    title = "Posterior Mean of Spatial b-value",
    x = "x (scaled longitude)", y = "y (scaled latitude)"
  ) +
  theme_minimal()

# spatial uncertainty map
ggplot(pred_df, aes(x = x, y = y, fill = b_sd)) +
  geom_tile() +
  coord_equal() +
  scale_fill_viridis_c(option = "magma", name = "Posterior SD") +
  labs(
    title = "Posterior Uncertainty of b-value (Standard Deviation)",
    x = "x (scaled longitude)", y = "y (scaled latitude)"
  ) +
  theme_minimal()

#95% Credible Interval Width map
pred_df$ci_width <- pred_df$b_upper - pred_df$b_lower

ggplot(pred_df, aes(x = x, y = y, fill = ci_width)) +
  geom_tile() +
  coord_equal() +
  scale_fill_viridis_c(option = "inferno", name = "95% CI Width") +
  labs(
    title = "Credible Interval Width for b-value",
    x = "x (scaled longitude)", y = "y (scaled latitude)"
  ) +
  theme_minimal()















