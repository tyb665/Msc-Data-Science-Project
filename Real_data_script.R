

#####################
#Plot for background###
#####################


#######
#Real data GR-law example
library(ggplot2)
library(readr)
library(dplyr)

# read data
data <- read_csv("D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv")

# Set the lower limit of the seismic magnitude (M0)
M0 <- 2.5

# Select the seismic magnitude that is non-missing and satisfies M >= M0
mags_values <- data %>% 
  filter(!is.na(magnitude), magnitude >= M0) %>% 
  pull(magnitude)

# Create magnitude bins
mags_bins <- seq(M0, max(mags_values), by = 0.1)

# Calculate the cumulative frequency N(M)
mags_counts <- vapply(mags_bins, \(bin) sum(mags_values >= bin), 0)

# Calculate log10 N(M)
log_counts <- log10(mags_counts)

# Fit the linear model
fit <- lm(log_counts ~ mags_bins)
slope <- coef(fit)[2]
b_value_est <- abs(round(slope, 3))
line_value <- predict(fit)

# Plot
ggplot() +
  geom_point(aes(x = mags_bins, y = log_counts), size = 2) +
  geom_line(aes(x = mags_bins, y = line_value), color = "red") +
  labs(
    title = "Log-Frequency Plot of Earthquake Magnitudes",
    subtitle = paste0("Linear fit slope = ", round(slope, 3), " → b ≈ ", b_value_est),
    x = "Magnitude",
    y = expression(log[10](N(M)))
  ) +
  theme_minimal(base_size = 14)

##############################
#####Gaussian process beta(t)
##############################
library(lubridate)
library(dplyr)
library(inlabru)
library(INLA)
library(ggplot2)

# Read data
data <- read.csv("D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv")

# Preprocessing
data_clean <- data %>%
  filter(magnitude >= 2.5) %>%
  mutate(
    time = ymd_hms(time),
    time_numeric = year(time) + yday(time) / 365.25  # 连续时间变量
  )

# Construct a one-dimensional temporal mesh
time_range <- range(data_clean$time_numeric)
time_mesh <- inla.mesh.1d(
  seq(time_range[1], time_range[2], length.out = 100)
)

# SPDE model
spde_time <- inla.spde2.pcmatern(
  mesh = time_mesh,
  prior.range = c(0.5, 0.9),  # Smoothness: P(range < 0.5) = 0.1
  prior.sigma = c(1, 0.01)    
)

# Define components
cmp <- ~ Intercept(1, model = "linear") +
  beta_field(time_numeric, model = spde_time)


# Define log-likelihood function
loglike <- function(beta, data, M0 = 2.5) {
  sum(log(beta) - beta * (data$magnitude - M0))
}

# Fitting (the number of events is regarded as point process)
bru_formula <- magnitude ~ Intercept + beta_field

# Construct inlabru likelihood
lik <- like(
  formula = bru_formula,
  family = "exponential",
  data = data_clean,
  control.family = list(
    control.link = list(model = "log")
  )
)

fit <- bru(
  formula = bru_formula,
  components = cmp,
  data = data_clean,
  family = "exponential",
  control.family = list(control.link = list(model = "log")),
  options = list(verbose = TRUE)
)



# ---------- Prediction: transform to b-value(t) ----------
# Construct prediction time points
pred_time <- data.frame(time_numeric = seq(time_range[1], time_range[2], length.out = 300))

# predict: return the posterior statistic of the expression we provided
# Here the log beta (t) - > beta (t) = exp (log beta (t)), again into b (t) = beta (t)/ln (10)
pred_b <- predict(
  fit,
  pred_time,
  ~ exp(Intercept + beta_field) / log(10),  
  n.samples = 1000
)

# For readability, give the columns a name (mean/sd/ quantile column names remain the default in inlabru)
names(pred_b)[names(pred_b) == "mean"]   <- "b_mean"
names(pred_b)[names(pred_b) == "sd"]     <- "b_sd"
names(pred_b)[names(pred_b) == "q0.025"] <- "b_q0.025"
names(pred_b)[names(pred_b) == "q0.975"] <- "b_q0.975"

# ---------- Plot b(t) ----------
library(ggplot2)

ggplot(pred_b, aes(x = time_numeric)) +
  geom_line(aes(y = b_mean), color = "blue") +
  geom_ribbon(aes(ymin = b_q0.025, ymax = b_q0.975), alpha = 0.2, fill = "skyblue") +
  labs(
    title = expression(Temporal~evolution~of~b(t)),
    x = "Time (year)",
    y = expression(b(t))
  ) +
  xlim(1932, 2022) +   
  theme_minimal()



##############################
#####Gaussian process beta(s)
##############################
library(dplyr)
library(INLA)
library(inlabru)
library(ggplot2)
library(maps)

# read data
data <- read.csv("D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv")

# Filter reliable data: M >= 2.5 and the longitude and latitude are not empty
data_clean <- data %>%
  filter(magnitude >= 2.5, !is.na(latitude), !is.na(longitude)) %>%
  mutate(x = longitude, y = latitude)

# Generate spatial mesh (based on the epicenter)
coords <- cbind(data_clean$x, data_clean$y)

mesh <- inla.mesh.2d(
  loc = coords,
  max.edge = c(0.5, 5),  # It can be adjusted depending on the density of the data
  cutoff = 0.2,
  offset = c(1, 2)
)

# Construct SPDE model
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(1, 0.01),    # P(range < 1 deg) = 0.01
  prior.sigma = c(1, 0.01)     # P(sigma > 1) = 0.01
)

# Define component as "spatial_field"
cmp <- ~ Intercept(1, model = "linear") +
  spatial_field(coordinates, model = spde)

# formula
formula <- magnitude ~ Intercept + spatial_field

# Define coordinate data (for matching the spatial index of SPDE)
data_clean$coordinates <- coords

# Fitting: Modeling magnitude (as proxy for log β(s))
fit <- bru(
  formula = formula,
  components = cmp,
  data = data_clean,
  family = "exponential",  
  control.family = list(control.link = list(model = "log")),
  options = list(verbose = TRUE)
)

# ---------- Prediction: transform to b-value(s) ----------
# Create a longitude and latitude grid covering the observation range (grid_res can be tuned according to performance)
x_range <- range(data_clean$x)
y_range <- range(data_clean$y)

pred_grid <- expand.grid(
  x = seq(x_range[1], x_range[2], length.out = grid_res),
  y = seq(y_range[1], y_range[2], length.out = grid_res)
)
pred_grid$coordinates <- cbind(pred_grid$x, pred_grid$y)

# Perform the transformation directly in the expression of predict：b(s) = exp(Intercept + spatial_field) / log(10)
pred_b <- predict(
  fit,
  pred_grid,
  ~ exp(Intercept + spatial_field) / log(10),
  n.samples = 1000
)

# Change to a clearer column name
names(pred_b)[names(pred_b) == "mean"]   <- "b_mean"
names(pred_b)[names(pred_b) == "sd"]     <- "b_sd"
names(pred_b)[names(pred_b) == "q0.025"] <- "b_q0.025"
names(pred_b)[names(pred_b) == "q0.975"] <- "b_q0.975"

# ---------- Plot b(s) ----------
usa_map <- map_data("state")
california_map <- subset(usa_map, region == "california")

ggplot() +
  # Heat map (posterior mean of b(s))
  geom_tile(data = pred_b, aes(x = x, y = y, fill = b_mean)) +
  
  # Shock points
  geom_point(
    data = data_clean,
    aes(x = x, y = y, shape = "Earthquake Epicenter"),
    color = "black", alpha = 0.3, size = 0.5
  ) +
  
  # The border of California
  geom_path(
    data = california_map,
    aes(x = long, y = lat, group = group, linetype = "California Border"),
    color = "black", linewidth = 0.5
  ) +
  
  # Color codes and legends
  scale_fill_viridis_c(option = "plasma", name = expression(b(s))) +
  scale_shape_manual(name = "", values = c("Earthquake Epicenter" = 16)) +
  scale_linetype_manual(name = "", values = c("California Border" = "solid")) +
  guides(
    shape = guide_legend(override.aes = list(size = 2, alpha = 1)),
    linetype = guide_legend(override.aes = list(linewidth = 1))
  ) +
  
  labs(
    title = expression("Posterior mean of " * b(s)),
    x = "Longitude", y = "Latitude"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )



#The section below created the picture for introduction
################################################################
########Distribution of magnitude on time and space########
################################################################
# =========================
# Time vs Magnitude (Intro)
# =========================
library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(scales)

df <- read_csv("D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv", show_col_types = FALSE) %>%
  mutate(time = ymd_hms(time, quiet = TRUE),
         date = as.Date(time)) %>%
  filter(!is.na(date), !is.na(magnitude))

tmin <- as.Date("1932-01-01")
tmax <- as.Date("2022-12-31")
M0  <- 2.5

# Calculate the annotation position: 92% on the right side, slightly above the dotted line
x_lab <- tmin + round(as.numeric(tmax - tmin) * 0.92)
y_lab <- M0 + 0.12

p_bin <- ggplot(df, aes(x = date, y = magnitude)) +
  geom_bin2d(binwidth = c(365.25, 0.1)) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "Count\n(log scale)",
    guide = guide_colorbar(barheight = unit(60, "pt"))
  ) +
  # White dotted line (M0)
  geom_hline(yintercept = M0, linetype = "dashed", color = "white", linewidth = 0.6) +
  annotate("label",
           x = x_lab, y = y_lab,
           label = sprintf("M0 = %.1f", M0),
           label.size = 0,
           fill = alpha("black", 0.35),
           colour = "white", size = 3.8) +
  labs(
    title = "Earthquake magnitudes over time",
    subtitle = "dashed line shows M0 = 2.5",
    x = "Year", y = "Magnitude"
  ) +
  scale_x_date(
    limits = c(tmin, tmax),
    date_breaks = "10 years",
    date_labels = "%Y",
    expand = c(0.01, 0.01)
  ) +
  coord_cartesian(ylim = range(df$magnitude, na.rm = TRUE)) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

p_bin

#################################################
# Space vs magnitude
#################################################
library(readr)
library(dplyr)
library(ggplot2)
library(viridisLite)  
library(maps)

# data
df <- read_csv("D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv", show_col_types = FALSE) |>
  filter(!is.na(latitude), !is.na(longitude), !is.na(magnitude), magnitude >= 2.5) |>
  mutate(x = longitude, y = latitude)

usa_map <- map_data("state")
ca_map  <- subset(usa_map, region == "california")

# Hexagonal density (bins can be fine-tuned as needed)
p_hex <- ggplot(df, aes(x = x, y = y)) +
  geom_hex(bins = 70) +
  scale_fill_viridis_c(trans = "log10", name = "Count\n(log scale)") +
  geom_point(data = df[sample.int(nrow(df), min(3000, nrow(df))), ],
             aes(x = x, y = y), color = "black", alpha = 0.15, size = 0.2) +
  geom_path(data = ca_map, aes(x = long, y = lat, group = group),
            color = "black", linewidth = 0.5) +
  labs(
    title = "Spatial distribution of seismicity (M ≥ 2.5)",
    x = "Longitude", y = "Latitude"
  ) +
  coord_quickmap() +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

p_hex
