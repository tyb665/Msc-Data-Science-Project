# ========================================================
# ======= Isolate different sequences with m0 ≥ 6.5 ======
# ========================================================
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)
Sys.setlocale("LC_TIME", "English")

# ==== Read combined data ====
file_path <- "D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv"
eq_data <- read_csv(file_path)

eq_data <- eq_data %>%
  rename(mag = magnitude) %>%
  mutate(time = as.POSIXct(time, tz = "UTC")) %>%
  drop_na(time, mag, latitude, longitude)

# ==== Extract all the main shocks (M ≥ 6.5) ====
mainshocks_all <- eq_data %>%
  filter(mag >= 6.5) %>%
  arrange(time)

# ==== Set the window size ====
window_days <- 60
lat_window <- 1.5
lon_window <- 1.5

filtered_mainshocks <- mainshocks_all[1, ]
for (i in 2:nrow(mainshocks_all)) {
  new_shock <- mainshocks_all[i, ]
  
  overlap <- any(
    abs(difftime(filtered_mainshocks$time, new_shock$time, units = "days")) <= window_days &
      abs(filtered_mainshocks$latitude - new_shock$latitude) <= lat_window &
      abs(filtered_mainshocks$longitude - new_shock$longitude) <= lon_window
  )
  
  if (!overlap) {
    filtered_mainshocks <- bind_rows(filtered_mainshocks, new_shock)
  }
}

# ==== Extract the sequences of all the main shocks ====
sequence_list <- list()
plot_list <- list()

for (i in seq_len(nrow(filtered_mainshocks))) {
  mainshock <- filtered_mainshocks[i, ]
  t0 <- mainshock$time
  lat0 <- mainshock$latitude
  lon0 <- mainshock$longitude
  
  seq_data <- eq_data %>%
    filter(
      time >= t0 - days(window_days),
      time <= t0 + days(window_days),
      abs(latitude - lat0) <= lat_window,
      abs(longitude - lon0) <= lon_window
    )
  # Skip if there is no data
  if (nrow(seq_data) == 0) next
  
  sequence_list[[i]] <- list(
    mainshock_time = t0,
    mainshock_mag = mainshock$mag,
    sequence_data = seq_data
  )
  
  # === Find the main shock point (with exactly the same time and magnitude)===
  mainshock_point <- seq_data %>%
    filter(time == t0 & mag == mainshock$mag)
  
  # subplots
  p <- ggplot(seq_data, aes(x = time, y = mag)) +
    geom_point(color = "orange", alpha = 0.7, size = 0.1) +
    geom_point(data = mainshock_point, aes(x = time, y = mag), color = "red", size = 1) +
    geom_vline(xintercept = as.numeric(t0), color = "red", linetype = "dashed") +
    labs(
      title = paste0("Sequence ", i, " | M", round(mainshock$mag, 1),
                     " @ ", format(t0, "%Y-%m-%d")),
      x = "Time", y = "Magnitude"
    ) +
    theme_minimal(base_size = 11)
  
  plot_list[[i]] <- p
}

# ==== Plot all sequences ====
# combined_all <- wrap_plots(plot_list, ncol = 2)
# print(combined_all)

# ==== Filter the main shocks after 1990 from the sequence_list ====
recent_indices <- which(sapply(sequence_list, function(seq) seq$mainshock_time >= as.POSIXct("1990-01-01")))
plot_list_recent <- plot_list[recent_indices]

# ==== Plot the sequence after 1990 ====
combined_recent <- wrap_plots(plot_list_recent, ncol = 2)
print(combined_recent)

##############
# Omori method
##############
##############
# Omori method (with fallback & symmetrical windows)
##############

# ==== Libraries ====
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(minpack.lm)

Sys.setlocale("LC_TIME", "English")

# ==== Load data ====
eq_data <- read_csv("D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv") %>%
  rename(mag = magnitude) %>%
  mutate(time = as.POSIXct(time, tz = "UTC")) %>%
  drop_na(time, mag, latitude, longitude)

# ==== Extract mainshocks (M ≥ 6.5, post-1990) ====
mainshocks_all <- eq_data %>%
  filter(mag >= 6.5, time >= as.POSIXct("1990-01-01")) %>%
  arrange(time)

# ==== De-duplicate mainshocks ====
window_days <- 60
lat_window <- 1.5
lon_window <- 1.5

filtered_mainshocks <- mainshocks_all[1, ]
for (i in 2:nrow(mainshocks_all)) {
  new_shock <- mainshocks_all[i, ]
  overlap <- any(
    abs(difftime(filtered_mainshocks$time, new_shock$time, units = "days")) <= window_days &
      abs(filtered_mainshocks$latitude - new_shock$latitude) <= lat_window &
      abs(filtered_mainshocks$longitude - new_shock$longitude) <= lon_window
  )
  if (!overlap) {
    filtered_mainshocks <- bind_rows(filtered_mainshocks, new_shock)
  }
}

# ==== Omori fit function with fallback ====
estimate_omori_window <- function(seq_data, t0, max_day = 120, decay_ratio = 20, fallback_day = 60, seq_id = NULL) {
  daily_counts <- seq_data %>%
    filter(time > t0) %>%
    mutate(day = as.integer(difftime(time, t0, units = "days"))) %>%
    count(day) %>%
    complete(day = 1:max_day, fill = list(n = 0))
  
  tryCatch({
    fit <- nlsLM(n ~ K / (day + c)^p,
                 data = daily_counts %>% filter(n > 0),
                 start = list(K = max(daily_counts$n), c = 0.5, p = 1.0),
                 lower = c(1, 0.01, 0.6),
                 upper = c(1e5, 5, 2.5),
                 control = list(maxiter = 1000))
    
    daily_counts$pred <- predict(fit, newdata = daily_counts)
    
    lambda_0 <- predict(fit, newdata = data.frame(day = 1))
    target_value <- lambda_0 / decay_ratio
    
    cutoff_day <- min(daily_counts$day[daily_counts$pred < target_value], na.rm = TRUE)
    
    if (!is.finite(cutoff_day)) {
      cutoff_day <- fallback_day
      cutoff_type <- "fallback"
    } else {
      cutoff_type <- "omori"
    }
    
    cat(sprintf("✅ Seq %s: K = %.1f, c = %.2f, p = %.2f | λ₀ = %.1f → cutoff = %s (%s)\n",
                seq_id, coef(fit)[["K"]], coef(fit)[["c"]], coef(fit)[["p"]],
                lambda_0,
                paste0(cutoff_day, " days"), cutoff_type))
    
    list(
      cutoff = cutoff_day,
      cutoff_type = cutoff_type,
      fit = fit,
      data = daily_counts
    )
  }, error = function(e) {
    cat(sprintf("❌ Seq %s: Fit failed - %s → fallback to %d days\n", seq_id, e$message, fallback_day))
    list(
      cutoff = fallback_day,
      cutoff_type = "fallback_fit_failed",
      fit = NULL,
      data = daily_counts
    )
  })
}

# ==== Main loop ====
sequence_list <- list()
cutoff_types <- c()

for (i in seq_len(nrow(filtered_mainshocks))) {
  mainshock <- filtered_mainshocks[i, ]
  
  # The original main shock time is used for display
  t0_orig <- mainshock$time
  mag0 <- mainshock$mag
  time_win_before <- 30
  
  # The main shock time used for window calculation (truncated throughout the day)
  t0 <- as.POSIXct(as.Date(t0_orig), tz = "UTC")
  
  full_window <- eq_data %>%
    filter(time > t0, time <= t0 + days(120))
  
  if (nrow(full_window) < 30) next
  
  omori_result <- estimate_omori_window(full_window, t0, seq_id = i)
  time_win_after <- omori_result$cutoff
  cutoff_type <- omori_result$cutoff_type
  
  seq_data <- eq_data %>%
    filter(time >= t0 - days(time_win_before),
           time <= t0 + days(time_win_after))
  
  if (nrow(seq_data) == 0) next
  
  cutoff_types <- c(cutoff_types, cutoff_type)
  
  sequence_list[[i]] <- list(
    mainshock_time = t0_orig,
    mag = mag0,
    time_win_before = time_win_before,
    time_win_after = time_win_after,
    cutoff_type = cutoff_type,
    sequence_data = seq_data
  )
}

# ==== Fallback override: force all to ±60 days if any fallback ====
if (any(cutoff_types != "omori")) {
  message("⚠️ Omori fit failed in some sequences. Forcing all sequences to ±60 days.")
  for (j in seq_along(sequence_list)) {
    t0_j <- as.POSIXct(as.Date(sequence_list[[j]]$mainshock_time), tz = "UTC")
    sequence_list[[j]]$time_win_before <- 30
    sequence_list[[j]]$time_win_after <- 60
    sequence_list[[j]]$cutoff_type <- "forced_fixed"
    sequence_list[[j]]$sequence_data <- eq_data %>%
      filter(time >= t0_j - days(30), time <= t0_j + days(60))
  }
}

# ==== Plotting ====
plot_list <- list()
for (i in seq_along(sequence_list)) {
  seq_info <- sequence_list[[i]]
  seq_data <- seq_info$sequence_data
  t0_orig <- seq_info$mainshock_time
  mag0 <- seq_info$mag
  time_win_before <- seq_info$time_win_before
  time_win_after <- seq_info$time_win_after
  cutoff_type <- seq_info$cutoff_type
  n_events <- nrow(seq_data)
  
  mainshock_point <- seq_data %>% filter(time == t0_orig & mag == mag0)
  
  subtitle <- paste0("±", time_win_before, "/", time_win_after, " days",
                     if (!is.null(cutoff_type)) paste0(" (", cutoff_type, ")") else "")
  
  p <- ggplot(seq_data, aes(x = time, y = mag)) +
    geom_point(color = "orange", alpha = 0.7, size = 0.1) +
    geom_point(data = mainshock_point, aes(x = time, y = mag), color = "red", size = 1) +
    geom_vline(xintercept = as.numeric(t0_orig), color = "red", linetype = "dashed") +
    annotate("text", x = max(seq_data$time), y = max(seq_data$mag, na.rm = TRUE) + 0.3,
             label = paste0("n = ", n_events), hjust = 1, size = 3.5) +
    labs(
      title = paste0("Seq ", i, " | M", round(mag0, 1), " @ ", format(t0_orig, "%Y-%m-%d")),
      subtitle = subtitle,
      x = "Time", y = "Magnitude"
    ) +
    theme_minimal(base_size = 11) +
    coord_cartesian(xlim = c(t0 - days(time_win_before), t0 + days(time_win_after)))

  plot_list[[i]] <- p
}

# ==== Show final plots ====
wrap_plots(plot_list, ncol = 2)



################################################################################
################################################################################
# Use the whole data of 1990-2021, m0 = 2.5, compare 3 strategies of mesh point
################################################################################
################################################################################
# ==== Libraries ====
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(INLA)
library(inlabru)
library(patchwork)

# ==== Set environment ====
Sys.setlocale("LC_TIME", "English")
Sys.setenv(TMPDIR = "D:/R_temp")  # Set a temporary directory 

# ==== Load data ====
file_path <- "D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv"
eq_data <- read_csv(file_path) %>%
  rename(mag = magnitude) %>%
  drop_na(time, mag, latitude, longitude) %>%
  mutate(time = as.POSIXct(time, tz = "UTC")) %>%
  filter(time >= as.POSIXct("1990-01-01"),
         time <= as.POSIXct("2021-12-31"))

# ==== Filter for m0 = 2.5 ====
m0 <- 2.5
eq_m <- eq_data %>%
  filter(mag >= m0) %>%
  mutate(
    day = as.numeric(difftime(time, min(time), units = "days")),
    mags = mag - m0
  )

# ==== Define 3 mesh types ====

# -- Mesh 1: Uniform spacing --
mesh1 <- fm_mesh_1d(seq(min(eq_m$day), max(eq_m$day), by = 100),
                    cutoff = 30, degree = 2, boundary = "free")

# -- Mesh 2: KDE-weighted sampling --
set.seed(42)
kde <- density(eq_m$day, bw = 200)
mesh_points2 <- sample(kde$x, size = 120, prob = kde$y, replace = FALSE)
mesh_points2 <- sort(unique(c(min(eq_m$day), mesh_points2, max(eq_m$day))))
mesh2 <- fm_mesh_1d(mesh_points2, cutoff = 30, degree = 2, boundary = "free")

# -- Mesh 3: Manually segmented mesh --
dense1 <- seq(0, 3000, by = 30)
dense2 <- seq(3000, 7000, by = 70)
dense3 <- seq(7000, max(eq_m$day), by = 150)
mesh_points3 <- sort(unique(c(dense1, dense2, dense3)))
mesh3 <- fm_mesh_1d(mesh_points3, cutoff = 30, degree = 2, boundary = "free")

# ==== Fit & Predict Function ====
fit_and_predict <- function(mesh, label) {
  spde_model <- inla.spde2.pcmatern(
    mesh,
    prior.range = c(200, 0.01),
    prior.sigma = c(0.5, 0.01)
  )
  
  comp <- mags ~ field(day, model = spde_model) + Intercept(1)
  
  fit <- bru(
    components = comp,
    data = eq_m,
    family = "exponential",
    options = list(control.compute = list(dic = TRUE, waic = TRUE))
  )
  
  pred <- predict(fit, data.frame(day = seq(0, max(eq_m$day), by = 20)),
                  ~ exp(field + Intercept) / log(10), n.samples = 1000)
  
  pred$mesh_type <- label
  return(list(pred = pred, dic = fit$dic$dic, waic = fit$waic$waic, mesh = mesh))
}

# ==== Run models ====
res1 <- fit_and_predict(mesh1, "Uniform mesh")
res2 <- fit_and_predict(mesh2, "KDE-based mesh")
res3 <- fit_and_predict(mesh3, "Segmented mesh")

# ==== Top plot: Magnitude vs Time ====
p0 <- ggplot(eq_m, aes(x = time, y = mag)) +
  geom_point(alpha = 0.4, size = 0.7) +
  labs(
    title = "Magnitude over time",
    x = "Time", y = "Magnitude"
  ) +
  theme_minimal(base_size = 12)

# ==== Helper: Build b-value plot with mesh lines ====
build_b_plot <- function(pred, mesh, title) {
  mesh_lines <- data.frame(day = mesh$loc)
  
  ggplot(pred, aes(x = day, y = mean)) +
    geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = "grey80") +
    geom_line(color = "black") +
    geom_vline(data = mesh_lines, aes(xintercept = day), color = "red", alpha = 0.3) +
    labs(
      title = title,
      x = "Days since 1990-01-01", y = "b-value"
    ) +
    coord_cartesian(ylim = c(0.3, 1.7)) +
    theme_minimal(base_size = 12)
}

# ==== Bottom 3 plots with mesh overlays ====
p1 <- build_b_plot(res1$pred, res1$mesh, "Uniform mesh")
p2 <- build_b_plot(res2$pred, res2$mesh, "KDE-based mesh")
p3 <- build_b_plot(res3$pred, res3$mesh, "Segmented mesh")

# ==== Combine 4 plots vertically ====
final_plot <- p0 / p1 / p2 / p3
print(final_plot)

# ==== (Optional) Save to file ====
# ggsave("b_value_mesh_comparison.png", final_plot, width = 10, height = 14, dpi = 300)


##the compare results of mesh strategies
# ==== Libraries ====
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(INLA)
library(inlabru)

# ==== Set environment ====
Sys.setlocale("LC_TIME", "English")
Sys.setenv(TMPDIR = "D:/R_temp")

# ==== Load data ====
file_path <- "D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv"
eq_data <- read_csv(file_path) %>%
  rename(mag = magnitude) %>%
  drop_na(time, mag, latitude, longitude) %>%
  mutate(time = as.POSIXct(time, tz = "UTC")) %>%
  filter(time >= as.POSIXct("1990-01-01"),
         time <= as.POSIXct("2021-12-31"))

# ==== Filter for m0 = 2.5 ====
m0 <- 2.5
eq_m <- eq_data %>%
  filter(mag >= m0) %>%
  mutate(
    day = as.numeric(difftime(time, min(time), units = "days")),
    mags = mag - m0
  )

# ==== Define 3 mesh types ====

# -- Mesh 1: Uniform spacing --
mesh1 <- fm_mesh_1d(seq(min(eq_m$day), max(eq_m$day), by = 100),
                    cutoff = 30, degree = 2, boundary = "free")

# -- Mesh 2: KDE-weighted sampling --
set.seed(42)
kde <- density(eq_m$day, bw = 200)
mesh_points2 <- sample(kde$x, size = 120, prob = kde$y, replace = FALSE)
mesh_points2 <- sort(unique(c(min(eq_m$day), mesh_points2, max(eq_m$day))))
mesh2 <- fm_mesh_1d(mesh_points2, cutoff = 30, degree = 2, boundary = "free")

# -- Mesh 3: Manually segmented mesh --
dense1 <- seq(0, 3000, by = 30)
dense2 <- seq(3000, 7000, by = 70)
dense3 <- seq(7000, max(eq_m$day), by = 150)
mesh_points3 <- sort(unique(c(dense1, dense2, dense3)))
mesh3 <- fm_mesh_1d(mesh_points3, cutoff = 30, degree = 2, boundary = "free")

# ==== Fit function ====
fit_and_predict <- function(mesh, label) {
  spde_model <- inla.spde2.pcmatern(
    mesh,
    prior.range = c(200, 0.01),
    prior.sigma = c(0.5, 0.01)
  )
  
  comp <- mags ~ field(day, model = spde_model) + Intercept(1)
  
  fit <- bru(
    components = comp,
    data = eq_m,
    family = "exponential",
    options = list(control.compute = list(dic = TRUE, waic = TRUE))
  )
  
  pred <- predict(fit, data.frame(day = seq(0, max(eq_m$day), by = 20)),
                  ~ exp(field + Intercept) / log(10), n.samples = 1000)
  
  pred$mesh_type <- label
  return(list(pred = pred, dic = fit$dic$dic, waic = fit$waic$waic))
}

# ==== Run models ====
res1 <- fit_and_predict(mesh1, "Uniform mesh")
res2 <- fit_and_predict(mesh2, "KDE-based mesh")
res3 <- fit_and_predict(mesh3, "Segmented mesh")

# ==== Combine for plotting ====
pred_all <- bind_rows(res1$pred, res2$pred, res3$pred)

# ==== Plot ====
ggplot(pred_all, aes(x = day, y = mean)) +
  geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = "grey80") +
  geom_line(color = "black") +
  facet_wrap(~mesh_type, ncol = 1, scales = "free_x") +
  labs(
    title = "b-value over time (m0 = 2.5) using different mesh strategies",
    x = "Days since 1990-01-01", y = "b-value"
  ) +
  coord_cartesian(ylim = c(0.3, 1.7)) +
  theme_minimal(base_size = 12)

# ==== Show model comparison ====
tibble(
  Mesh = c("Uniform", "KDE-based", "Segmented"),
  DIC = c(res1$dic, res2$dic, res3$dic),
  WAIC = c(res1$waic, res2$waic, res3$waic)
)


################################################################################
################################################################################
##### Try different m0########
################################################################################
################################################################################
# ==== Library ====
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(INLA)
options(inla.call = "remote")
library(inlabru)

# ==== Set environment ====
Sys.setlocale("LC_TIME", "English")
Sys.setenv(TMPDIR = "D:/R_temp")

# ==== Read combined data ====
file_path <- "D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv"
eq_data <- read_csv(file_path) %>%
  rename(mag = magnitude) %>%
  drop_na(time, mag, latitude, longitude) %>%
  mutate(time = as.POSIXct(time, tz = "UTC")) %>%
  filter(time >= as.POSIXct("1990-01-01"),
         time <= as.POSIXct("2021-12-31"))


# ==== m0 list ====
m0_list <- c(2.5, 2.8, 3.0, 3.2, 3.5, 4.0, 4.5)
fit_results <- list()

# ==== Loop over m0 ====
for (m0 in m0_list) {
  cat(">>> Fitting model with m0 =", m0, "\n")
  
  eq_m <- eq_data %>%
    filter(mag >= m0) %>%
    mutate(
      day = as.numeric(difftime(time, min(time), units = "days")),
      mags = mag - m0
    )
  
  if (nrow(eq_m) < 100) {
    message("Skipped m0 = ", m0, " due to insufficient data")
    next
  }
  
  eq_m <- eq_m[order(eq_m$day),]
  eq_m$day <- eq_m$day + 0.001  # avoid duplicate zero-day
  
  # Mesh
  extremes <- range(eq_m$day)
  mesh_points <- seq(extremes[1], extremes[2], by = 100)
  mesh <- fm_mesh_1d(mesh_points, cutoff = 30, degree = 2, boundary = "free")
  
  # SPDE model
  spde_model <- inla.spde2.pcmatern(
    mesh,
    prior.range = c(200, 0.01),
    prior.sigma = c(0.5, 0.01)
  )
  
  comp <- mags ~ field(day, model = spde_model) + Intercept(1)
  
  fit <- tryCatch({
    bru(
      components = comp,
      data = eq_m,
      family = "exponential",
      options = list(control.compute = list(dic = TRUE, waic = TRUE))
    )
  }, error = function(e) {
    message("Model failed for m0 = ", m0, ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(fit)) {
    pred_df <- predict(fit, data.frame(day = seq(0, max(eq_m$day), by = 20)),
                       ~ exp(field + Intercept) / log(10), n.samples = 1000)
    
    fit_results[[as.character(m0)]] <- list(
      m0 = m0,
      dic = fit$dic$dic,
      waic = fit$waic$waic,
      pred = pred_df
    )
  }
}

# ==== Summary of model fit ====
summary_table <- bind_rows(lapply(fit_results, function(res) {
  tibble(
    m0 = res$m0,
    DIC = round(res$dic, 2),
    WAIC = round(res$waic, 2)
  )
})) |> arrange(DIC)

print(summary_table)

# ==== Plotting b-value curves ====
plots <- list()
for (i in seq_along(fit_results)) {
  m0 <- names(fit_results)[i]
  pred <- fit_results[[i]]$pred
  p <- ggplot(pred, aes(x = day)) +
    geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = "grey80") +
    geom_line(aes(y = mean), color = "black", size = 1) +
    labs(title = paste("m0 =", m0), y = "b-value", x = "days since 1990-01-01") +
    ylim(0.3, 1.7) +
    scale_x_continuous(
      breaks = seq(0, max(pred$day), by = 365 * 5),
      labels = function(x) format(as.Date("1990-01-01") + x, "%Y")
    )
  plots[[m0]] <- p
}

patchwork::wrap_plots(plots, ncol = 2)


###################################
###### Try m0 = 2.5, 3.5, 4.5#####
###################################
# ==== Libraries ====
# ==== Libraries ====
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(INLA)
library(inlabru)
library(patchwork)

# ==== Set environment ====
Sys.setlocale("LC_TIME", "English")
Sys.setenv(TMPDIR = "D:/R_temp")

# ==== Read cleaned data ====
file_path <- "D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv"
eq_data <- read_csv(file_path) %>%
  rename(mag = magnitude) %>%
  drop_na(time, mag, latitude, longitude) %>%
  mutate(time = as.POSIXct(time, tz = "UTC")) %>%
  filter(time >= as.POSIXct("1990-01-01"),
         time <= as.POSIXct("2021-12-31"))

# ==== Define m0 list ====
m0_list <- c(2.5, 3.5, 4.5)

# ==== Store results ====
fit_results <- list()
plots <- list()

# ==== Loop over m0 values ====
for (m0 in m0_list) {
  cat(">>> Fitting model with m0 =", m0, "\n")
  
  # Filter catalog
  eq_m <- eq_data %>%
    filter(mag >= m0) %>%
    mutate(
      day = as.numeric(difftime(time, min(time), units = "days")),
      mags = mag - m0
    )
  
  if (nrow(eq_m) < 100) {
    message("Skipped m0 = ", m0, " due to insufficient data")
    next
  }
  
  eq_m <- eq_m[order(eq_m$day), ]
  eq_m$day <- eq_m$day + 0.001  # avoid 0-day duplication
  
  # Create mesh
  mesh_points <- seq(min(eq_m$day), max(eq_m$day), by = 100)
  mesh <- fm_mesh_1d(mesh_points, cutoff = 30, degree = 2, boundary = "free")
  
  # Define SPDE model
  spde_model <- inla.spde2.pcmatern(
    mesh,
    prior.range = c(200, 0.01),
    prior.sigma = c(0.5, 0.01)
  )
  
  # Model components
  comp <- mags ~ field(day, model = spde_model) + Intercept(1)
  
  # Fit model
  fit <- tryCatch({
    bru(
      components = comp,
      data = eq_m,
      family = "exponential",
      options = list(control.compute = list(dic = TRUE, waic = TRUE))
    )
  }, error = function(e) {
    message("❌ Model failed for m0 = ", m0, ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(fit)) {
    pred_df <- predict(fit, data.frame(day = seq(0, max(eq_m$day), by = 20)),
                       ~ exp(field + Intercept) / log(10), n.samples = 1000)
    
    # Store results
    fit_results[[as.character(m0)]] <- list(
      m0 = m0,
      dic = fit$dic$dic,
      waic = fit$waic$waic,
      pred = pred_df
    )
    
    # Plot
    # Obtain the current sample size N
    N <- nrow(filter(eq_data, mag >= m0))
    
    # Plot
    p <- ggplot(pred_df, aes(x = day)) +
      geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = "grey80") +
      geom_line(aes(y = mean), color = "black", size = 1) +
      labs(
        title = paste0("b-value over time | m0 = ", m0, " (N = ", N, ")"),
        y = "b-value", x = "Days since 1990-01-01"
      ) +
      ylim(0.3, 1.7) +
      scale_x_continuous(
        breaks = seq(0, max(pred_df$day), by = 365 * 5),
        labels = function(x) format(as.Date("1990-01-01") + x, "%Y")
      ) +
      theme_minimal()
    
    
    plots[[as.character(m0)]] <- p
  }
}

# ==== Display plots ====
wrap_plots(plots, ncol = 1)

# ==== Summary table ====
summary_table <- bind_rows(lapply(fit_results, function(res) {
  tibble(
    m0 = res$m0,
    N = nrow(filter(eq_data, mag >= res$m0)),
    DIC = round(res$dic, 2),
    WAIC = round(res$waic, 2)
  )
})) %>% arrange(m0)

print(summary_table)



################################################################################
####Try different priors for the range and the std with KDE-based mesh####################
################################################################################
library(readr)
library(dplyr)
library(tibble)
library(lubridate)
library(ggplot2)
library(inlabru)
library(INLA)
library(patchwork)

# ==== 1. Load and filter real data ====
file_path <- "D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv"
eq_data <- read_csv(file_path) %>%
  rename(mag = magnitude) %>%
  drop_na(time, mag, latitude, longitude) %>%
  mutate(time = as.POSIXct(time, tz = "UTC")) %>%
  filter(time >= as.POSIXct("1990-01-01"),
         time <= as.POSIXct("2021-12-31"))

# ==== 2. Select m0 and prepare df ====
m0 <- 2.5
T2 <- as.numeric(difftime(as.POSIXct("2021-12-31"), as.POSIXct("1990-01-01"), units = "days"))

eq_m <- eq_data %>%
  filter(mag >= m0) %>%
  mutate(
    day = as.numeric(difftime(time, min(time), units = "days")),
    mags = mag - m0
  )

df_bru_real <- eq_m %>%
  select(time = day, mags)

# ==== 3. Construct KDE-based mesh ====
set.seed(42)
kde <- density(df_bru_real$time, bw = 300)
mesh_points_kde <- sample(kde$x, size = 120, prob = kde$y, replace = FALSE)
mesh_points_kde <- sort(unique(c(min(df_bru_real$time), mesh_points_kde, max(df_bru_real$time))))
mesh_kde <- fm_mesh_1d(mesh_points_kde, degree = 2, boundary = "free")

# ==== 4. Define prior ranges and sigmas ====
prior_range_list <- round(c(T2/1000, T2/100, T2/50, T2/20, T2/10, T2/5, T2/2, T2/1.5), 1)
prior_sigma_list <- c(0.01, 0.1, 0.3, 1, 3, 10)

# ==== 5. Grid search ====
results <- tibble(
  range = numeric(),
  sigma = numeric(),
  DIC = numeric(),
  WAIC = numeric(),
  range_mean = numeric(),
  stdev_mean = numeric()
)

for (r in prior_range_list) {
  for (s in prior_sigma_list) {
    cat("Trying: prior.range =", r, "| prior.sigma =", s, "\n")
    
    spde <- inla.spde2.pcmatern(
      mesh_kde,
      prior.range = c(r, 0.01),
      prior.sigma = c(s, 0.01)
    )
    
    comp_temp <- mags ~ field(time, model = spde) + Intercept(1)
    
    fit <- tryCatch(
      bru(
        components = comp_temp,
        bru_obs(mags ~ field + Intercept, data = df_bru_real, family = "exponential"),
        options = list(control.compute = list(dic = TRUE, waic = TRUE))
      ),
      error = function(e) {
        cat("Model failed:", conditionMessage(e), "\n")
        return(NULL)
      }
    )
    
    if (!is.null(fit) && !is.null(fit$dic$dic)) {
      dic_val <- round(fit$dic$dic, 2)
      waic_val <- round(fit$waic$waic, 2)
      range_mean <- round(fit$summary.hyperpar["Range for field", "mean"], 2)
      stdev_mean <- round(fit$summary.hyperpar["Stdev for field", "mean"], 2)
      
      cat("Success | DIC:", dic_val, "| WAIC:", waic_val, "\n")
      
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
      cat("Fitting completed but no DIC/WAIC extracted. Skipped.\n")
    }
  }
}

# ==== 6. Results table ====
results_sorted <- results %>% arrange(DIC)
print(results_sorted)

# ==== 7. Plot heatmaps ====
library(ggplot2)
library(dplyr)
library(patchwork)

# === Reconstruct the ploting data ===
results_heatmap <- results %>%
  mutate(
    range = factor(range, levels = sort(unique(range))),
    sigma = factor(sigma, levels = sort(unique(sigma)))
  )

# === Optimal setting (minimum DIC) ===
best_row <- results %>% arrange(DIC) %>% slice(1)

# === DIC hot map ===
p_dic <- ggplot(results_heatmap, aes(x = range, y = sigma, fill = DIC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(DIC, 1)), size = 3, color = "white") +
  geom_point(data = best_row, aes(x = factor(range), y = factor(sigma)), 
             shape = 8, size = 3, color = "cyan") +
  scale_fill_viridis_c(option = "inferno", direction = -1) +
  labs(title = "DIC across Prior Settings",
       x = "Prior Range", y = "Prior Sigma", fill = "DIC") +
  theme_minimal()

# === WAIC hot map ===
p_waic <- ggplot(results_heatmap, aes(x = range, y = sigma, fill = WAIC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(WAIC, 1)), size = 3, color = "white") +
  geom_point(data = best_row, aes(x = factor(range), y = factor(sigma)), 
             shape = 8, size = 3, color = "cyan") +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(title = "WAIC across Prior Settings",
       x = "Prior Range", y = "Prior Sigma", fill = "WAIC") +
  theme_minimal()

# === Display the final combination diagram ===
p_dic / p_waic







## Use the optimal mesh strategy and prior parameter to fit model with b-value curve
# ==== Reuse previous data & mesh ====

# Make sure eq_m, df_bru_real, mesh_kde has been constructed

# ==== Build SPDE model with optimal priors ====
best_spde <- inla.spde2.pcmatern(
  mesh_kde,
  prior.range = c(116.9, 0.01),
  prior.sigma = c(0.1, 0.01)
)

# ==== Define model components ====
best_comp <- mags ~ field(time, model = best_spde) + Intercept(1)

# ==== Fit the model ====
fit_best <- bru(
  components = best_comp,
  bru_obs(
    mags ~ field + Intercept,
    data = df_bru_real,
    family = "exponential"
  ),
  options = list(
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
  )
)

# ==== Prediction grid ====
pred_time <- data.frame(time = seq(0, max(df_bru_real$time), by = 20))

# ==== Predict b-value over time ====
pred_best <- predict(
  fit_best,
  pred_time,
  ~ exp(field + Intercept) / log(10),
  n.samples = 1000
)

# ==== Plot the predicted b-value curve ====
ggplot() +
  geom_ribbon(data = pred_best, aes(x = time, ymin = q0.025, ymax = q0.975), fill = "grey80") +
  geom_line(data = pred_best, aes(x = time, y = mean), color = "black", size = 1) +
  labs(
    title = "Final b-value fit (m0 = 2.5, KDE mesh, optimal priors)",
    x = "Days since 1990-01-01",
    y = "b-value"
  ) +
  coord_cartesian(ylim = c(0.3, 1.7)) +
  theme_minimal()

# ==== Summary of model quality ====
tibble(
  DIC = round(fit_best$dic$dic, 2),
  WAIC = round(fit_best$waic$waic, 2),
  Range_mean = round(inla.emarginal(identity, fit_best$marginals.hyperpar[["Range for field"]]), 2),
  Stdev_mean = round(inla.emarginal(identity, fit_best$marginals.hyperpar[["Stdev for field"]]), 2)
)


################################################################################
#For each sequence (after 1990), b-value fitting and comparative analysis before and after the main shock were conducted
################################################################################
# ==== Required libraries ====
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(tibble)
library(INLA)
library(inlabru)

# ==== Read and filter earthquake catalog from 1990 onward ====
eq_data <- read_csv("D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv") %>%
  rename(mag = magnitude) %>%
  drop_na(time, mag, latitude, longitude) %>%
  mutate(time = as.POSIXct(time, tz = "UTC")) %>%
  filter(time >= as.POSIXct("1990-01-01"))

################################################################################
# Step 0: Build sequence_list from 1990 onward with spatial & temporal filtering
################################################################################

mainshocks_all <- eq_data %>%
  filter(mag >= 6.5) %>%
  arrange(time)

filtered_mainshocks <- mainshocks_all[1, ]
for (i in 2:nrow(mainshocks_all)) {
  new_shock <- mainshocks_all[i, ]
  overlap <- any(
    abs(difftime(filtered_mainshocks$time, new_shock$time, units = "days")) <= 60 &
      abs(filtered_mainshocks$latitude - new_shock$latitude) <= 1.5 &
      abs(filtered_mainshocks$longitude - new_shock$longitude) <= 1.5
  )
  if (!overlap) {
    filtered_mainshocks <- bind_rows(filtered_mainshocks, new_shock)
  }
}

sequence_list <- list()
for (i in seq_len(nrow(filtered_mainshocks))) {
  mainshock <- filtered_mainshocks[i, ]
  t0 <- mainshock$time
  lat0 <- mainshock$latitude
  lon0 <- mainshock$longitude
  
  seq_data <- eq_data %>%
    filter(
      time >= t0 - days(90),
      time <= t0 + days(90),
      abs(latitude - lat0) <= 1.5,
      abs(longitude - lon0) <= 1.5
    )
  
  if (nrow(seq_data) > 0) {
    sequence_list[[i]] <- list(
      mainshock_time = t0,
      mainshock_mag = mainshock$mag,
      sequence_data = seq_data
    )
  }
}

# Filter sequence_list again (safety) to ensure 1990+
sequence_list <- Filter(function(seq) seq$mainshock_time >= as.POSIXct("1990-01-01"), sequence_list)

################################################################################
# For each sequence: b-value fitting and comparative analysis before and after mainshock
################################################################################

# ==== Function: KDE-based b-value fitting for one segment ====
# Unification: Global time origin
time0_global <- as.POSIXct("1990-01-01", tz = "UTC")

# Fit the entire sequence (once INLA/inlabru) and cut pre/post from the same curve
time0_global <- as.POSIXct("1990-01-01", tz = "UTC")

fit_sequence_once <- function(df_segment_full, t0, m0 = 2.5,
                              range_prior = 116.9, sigma_prior = 0.1,
                              min_n = 20, bw_kde = 30) {
  # Event count statistics (by m0
  df_events <- df_segment_full %>% dplyr::filter(mag >= m0)
  
  N_before_events <- df_events %>% dplyr::filter(time <  t0) %>% nrow()
  N_after_events  <- df_events %>% dplyr::filter(time >= t0) %>% nrow()
  
  # Fit data
  df <- df_events %>%
    dplyr::mutate(
      day  = as.numeric(difftime(time, time0_global, units = "days")),
      mags = mag - m0
    )
  if (nrow(df) < min_n) return(NULL)
  
  kde <- density(df$day, bw = bw_kde)
  mesh_points <- sample(kde$x, size = min(100, length(kde$x)), prob = kde$y, replace = FALSE)
  mesh_points <- sort(unique(c(min(df$day), mesh_points, max(df$day))))
  mesh <- fm_mesh_1d(mesh_points, degree = 2, boundary = "free")
  
  spde <- inla.spde2.pcmatern(mesh,
                              prior.range = c(range_prior, 0.01),
                              prior.sigma = c(sigma_prior, 0.01))
  comp <- mags ~ field(day, model = spde) + Intercept(1)
  
  fit <- tryCatch({
    bru(components = comp, data = df, family = "exponential",
        options = list(control.compute = list(dic=TRUE, waic=TRUE, config=TRUE)))
  }, error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  # Unified Curve prediction
  day_grid <- seq(min(df$day), max(df$day), by = 1)
  pred <- predict(fit, data.frame(day = day_grid),
                  ~ exp(field + Intercept) / log(10), n.samples = 1000)
  
  # Cut pre/post based on t0
  t0_day <- as.numeric(difftime(t0, time0_global, units = "days"))
  pred$segment <- ifelse(pred$day < t0_day, "Before", "After")
  
  # Segmented indicators (using the same curve)
  summarize_seg <- function(p) dplyr::tibble(
    mean_b = mean(p$mean, na.rm = TRUE),
    ci_low = quantile(p$mean, 0.025, na.rm = TRUE),
    ci_high= quantile(p$mean, 0.975, na.rm = TRUE)
  )
  sum_before <- summarize_seg(dplyr::filter(pred, segment=="Before"))
  sum_after  <- summarize_seg(dplyr::filter(pred, segment=="After"))
  
  list(
    fit = fit,
    pred = pred,
    summary = dplyr::tibble(
      b_before = sum_before$mean_b, ci_before_low = sum_before$ci_low,  ci_before_high = sum_before$ci_high,
      b_after  = sum_after$mean_b,  ci_after_low  = sum_after$ci_low,   ci_after_high  = sum_after$ci_high,
      N_before = N_before_events,    N_after      = N_after_events,
      diff     = sum_after$mean_b - sum_before$mean_b
    )
  )
}

# ==== Function: Analyze one mainshock sequence ====

analyze_sequence_once <- function(seq, m0 = 2.5, window_days = 90,
                                  range_prior = 116.9, sigma_prior = 0.1) {
  t0  <- seq$mainshock_time
  dat <- seq$sequence_data %>%
    filter(time >= t0 - days(window_days),
           time <= t0 + days(window_days))
  res <- fit_sequence_once(dat, t0 = t0, m0 = m0,
                           range_prior = range_prior, sigma_prior = sigma_prior)
  if (is.null(res)) return(NULL)
  
  res$summary <- res$summary %>%
    mutate(mainshock_time = t0,
           mag            = seq$mainshock_mag)
  res
}

# 主循环（替换你原来的“前后各拟合一次”的循环）
summary_table <- tibble()
plot_list <- list()

for (i in seq_along(sequence_list)) {
  cat("\n▶️ Seq", i, ":", format(sequence_list[[i]]$mainshock_time, "%Y-%m-%d"), "\n")
  res <- analyze_sequence_once(sequence_list[[i]], m0 = 2.5, window_days = 90)
  if (is.null(res)) {
    cat("❌ Sequence skipped\n"); next
  }
  summary_table <- bind_rows(summary_table, res$summary)
  
  pred <- res$pred
  sumry <- res$summary
  t0_day <- as.numeric(difftime(sequence_list[[i]]$mainshock_time, time0_global, units = "days"))
  
  # Take the main shock as 0 and the left and right sides as the symmetric number of days
  pred <- pred %>% dplyr::mutate(x_rel = day - t0_day)
  pred_b <- dplyr::filter(pred, segment == "Before")
  pred_a <- dplyr::filter(pred, segment == "After")
  
  xmin <- min(pred$x_rel, na.rm = TRUE)
  xmax <- max(pred$x_rel, na.rm = TRUE)
  
  title_date <- format(sequence_list[[i]]$mainshock_time, "%Y-%m-%d")
  title_mag  <- sprintf("M%.1f", sequence_list[[i]]$mainshock_mag)
  label_text <- sprintf("%s\n%s\nb_before = %.3f (N=%d)\nb_after  = %.3f (N=%d)",
                        title_date, title_mag,
                        sumry$b_before, sumry$N_before,
                        sumry$b_after,  sumry$N_after)
  
  y_top <- max(pred$q0.975, na.rm = TRUE) + 0.05
  
  p <- ggplot() +
    # BEFORE blue
    geom_ribbon(data = pred_b, aes(x = x_rel, ymin = q0.025, ymax = q0.975),
                fill = "steelblue", alpha = 0.18) +
    geom_line(data = pred_b, aes(x = x_rel, y = mean),
              color = "steelblue", linewidth = 1) +
    geom_rect(aes(xmin = min(pred_b$x_rel, na.rm=TRUE),
                  xmax = max(pred_b$x_rel, na.rm=TRUE),
                  ymin = sumry$ci_before_low,
                  ymax = sumry$ci_before_high),
              fill = "steelblue", alpha = 0.15) +
    geom_segment(aes(x = min(pred_b$x_rel, na.rm=TRUE),
                     xend = max(pred_b$x_rel, na.rm=TRUE),
                     y = sumry$b_before, yend = sumry$b_before),
                 color = "steelblue", linewidth = 1.1, lineend = "round") +
    
    # AFTER red
    geom_ribbon(data = pred_a, aes(x = x_rel, ymin = q0.025, ymax = q0.975),
                fill = "red", alpha = 0.18) +
    geom_line(data = pred_a, aes(x = x_rel, y = mean),
              color = "red", linewidth = 1) +
    geom_rect(aes(xmin = min(pred_a$x_rel, na.rm=TRUE),
                  xmax = max(pred_a$x_rel, na.rm=TRUE),
                  ymin = sumry$ci_after_low,
                  ymax = sumry$ci_after_high),
              fill = "red", alpha = 0.15) +
    geom_segment(aes(x = min(pred_a$x_rel, na.rm=TRUE),
                     xend = max(pred_a$x_rel, na.rm=TRUE),
                     y = sumry$b_after, yend = sumry$b_after),
                 color = "red", linewidth = 1.1, lineend = "round") +
    
    # The main shock vertical line (currently at x=0)
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    
    # Add information
    annotate("text", x = xmin + 0.02*(xmax - xmin), y = y_top,
             label = label_text, hjust = 0, vjust = 1, size = 3.6) +
    
    labs(x = "Days since sequence start", y = "b-value") +
    coord_cartesian(ylim = c(0.4, max(1.8, y_top))) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
  
  plot_list[[length(plot_list) + 1]] <- p
} 
 


# Show
plots_per_page <- 6
n_pages <- ceiling(length(plot_list) / plots_per_page)
for (page in seq_len(n_pages)) {
  start <- (page - 1) * plots_per_page + 1
  end   <- min(page * plots_per_page, length(plot_list))
  page_plot <- patchwork::wrap_plots(plot_list[start:end], ncol = 2)
  print(page_plot)
  # ggsave(paste0("sequence_page_singlefit_", page, ".png"), page_plot, width = 10, height = 8, dpi = 300)
}

# Print summary
print(summary_table %>%
        select(mainshock_time, mag,
               b_before, ci_before_low, ci_before_high, N_before,
               b_after,  ci_after_low,  ci_after_high,  N_after,
               diff))



##############################
##Draw a horizontal error bar graph
##############################
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Step 1: Clean the summary_table and organize it in long format
tidy_bvalues <- summary_table %>%
  mutate(seq_id = paste0("Seq ", row_number(), " @ ", substr(mainshock, 1, 10))) %>%
  rowwise() %>%
  mutate(
    before_low = as.numeric(str_extract(ci_before, "(?<=\\[)[^,]+")),
    before_high = as.numeric(str_extract(ci_before, "(?<=, )[^\\]]+")),
    after_low = as.numeric(str_extract(ci_after, "(?<=\\[)[^,]+")),
    after_high = as.numeric(str_extract(ci_after, "(?<=, )[^\\]]+"))
  ) %>%
  select(seq_id, b_before, before_low, before_high,
         b_after, after_low, after_high, diff) %>%
  pivot_longer(
    cols = c(b_before, b_after),
    names_to = "group", values_to = "b_value"
  ) %>%
  mutate(
    group = ifelse(group == "b_before", "Before", "After"),
    ci_low = ifelse(group == "Before", before_low, after_low),
    ci_high = ifelse(group == "Before", before_high, after_high)
  ) %>%
  select(seq_id, group, b_value, ci_low, ci_high, diff)

# Step 2: Add a diff label and a y value to each group of After
diff_labels <- tidy_bvalues %>%
  filter(group == "After") %>%
  mutate(
    label = paste0("Δ=", round(diff, 2)),
    y = factor(seq_id, levels = rev(unique(tidy_bvalues$seq_id)))  # Set the y order consistent with the graph
  )

# Step 3: Plot a graph and mark the difference
ggplot(tidy_bvalues, aes(x = b_value, y = factor(seq_id, levels = rev(unique(seq_id))), color = group)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25,
                 position = position_dodge(width = 0.5)) +
  geom_text(data = diff_labels,
            aes(x = b_value + 0.05, y = y, label = label),
            inherit.aes = FALSE,
            color = "black", size = 3.5, hjust = 0) +
  scale_color_manual(values = c("Before" = "blue", "After" = "red")) +
  labs(
    title = "Comparison of b-values Before and After Mainshocks",
    x = "b-value",
    y = "Mainshock Sequence",
    color = "Group"
  ) +
  theme_minimal(base_size = 13) +
  xlim(0.5, max(tidy_bvalues$ci_high, na.rm = TRUE) + 0.4)


################################################################################
################################################################################
######b-value covariate modeling: add depth#####################################
################################################################################
################################################################################
# ==== Libraries ====
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(INLA)
library(inlabru)

# ==== Parameters ====
m0 <- 2.5
range_prior <- 116.9
sigma_prior <- 0.1

# ==== Load full catalog data (1990+) ====
file_path <- "D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv"
eq_data <- read_csv(file_path) %>%
  rename(mag = magnitude) %>%
  drop_na(time, mag, latitude, longitude, depth) %>%
  mutate(time = as.POSIXct(time, tz = "UTC")) %>%
  filter(time >= as.POSIXct("1990-01-01"))

# ==== Filter by m0 and prepare variables ====
eq_m <- eq_data %>%
  filter(mag >= m0) %>%
  mutate(
    day = as.numeric(difftime(time, min(time), units = "days")),
    mags = mag - m0,
    depth_scaled = scale(depth)[,1]  # optional: standardized depth
  )

# ==== Construct mesh on time ====
mesh <- fm_mesh_1d(seq(min(eq_m$day), max(eq_m$day), by = 100), 
                   cutoff = 30, degree = 2, boundary = "free")

# ==== Define SPDE model ====
spde_model <- inla.spde2.pcmatern(
  mesh,
  prior.range = c(range_prior, 0.01),
  prior.sigma = c(sigma_prior, 0.01)
)

# ==== Model components: add depth as linear covariate ====
comp <- mags ~ Intercept(1) + depth_scaled + field(day, model = spde_model)

# ==== Fit the model ====
fit <- bru(
  components = comp,
  data = eq_m,
  family = "exponential",
  options = list(control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))
)

# ==== Posterior summaries ====
summary(fit)

# ==== Extract posterior of beta_depth ====
depth_post <- fit$marginals.fixed$depth_scaled

# Plot posterior of depth coefficient
plot(depth_post, type = "l", xlab = expression(beta[depth]), ylab = "Density", main = "Posterior of Depth Effect")
abline(v = inla.emarginal(identity, depth_post), col = "red", lty = 2)

# Print posterior mean and 95% CI
depth_mean <- inla.emarginal(identity, depth_post)
depth_ci <- inla.qmarginal(c(0.025, 0.975), depth_post)
cat("\nPosterior beta_depth:")
cat("\n  Mean     :", round(depth_mean, 4))
cat("\n  95% CI   :", round(depth_ci, 4), "\n")

# ==== Optionally visualize b-value over time (marginalized over depth) ====
pred_time <- data.frame(day = seq(0, max(eq_m$day), by = 20), depth_scaled = 0)

pred <- predict(fit, pred_time, ~ exp(field + Intercept) / log(10), n.samples = 1000)

ggplot(pred, aes(x = day)) +
  geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = "grey80") +
  geom_line(aes(y = mean), color = "black", size = 1) +
  labs(title = "b-value over time (m0 = 2.5, with depth covariate)",
       x = "Days since 1990-01-01", y = "b-value") +
  coord_cartesian(ylim = c(0.3, 1.7)) +
  theme_minimal()

###########
#for each sequence
# ==== Libraries ====
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(tibble)
library(INLA)
library(inlabru)
library(tidyr)
library(stringr)

# ==== Function: KDE-based b-value fitting for one segment with depth ====
# Unification: Global time origin
time0_global <- as.POSIXct("1990-01-01", tz = "UTC")

# Fit the entire sequence (once INLA/inlabru) 
time0_global <- as.POSIXct("1990-01-01", tz = "UTC")

# One-time fitting：log β(t) = Intercept + field(t) + β_depth * depth_scaled
fit_sequence_once_with_depth <- function(df_segment_full, t0,
                                         m0 = 2.5,
                                         range_prior = 116.9,
                                         sigma_prior = 0.1,
                                         min_n = 20,
                                         bw_kde = 30) {
  df_events <- df_segment_full %>% filter(mag >= m0) %>% drop_na(depth)
  N_before_events <- df_events %>% filter(time <  t0) %>% nrow()
  N_after_events  <- df_events %>% filter(time >= t0) %>% nrow()
  
  df <- df_events %>% mutate(
    day          = as.numeric(difftime(time, time0_global, units = "days")),
    mags         = mag - m0,
    depth_scaled = as.numeric(scale(depth)[, 1])
  )
  if (nrow(df) < min_n) return(NULL)
  
  kde <- density(df$day, bw = bw_kde)
  mesh_points <- sample(kde$x, size = min(100, length(kde$x)), prob = kde$y, replace = FALSE)
  mesh_points <- sort(unique(c(min(df$day), mesh_points, max(df$day))))
  mesh <- fm_mesh_1d(mesh_points, degree = 2, boundary = "free")
  
  spde <- inla.spde2.pcmatern(mesh,
                              prior.range = c(range_prior, 0.01),
                              prior.sigma = c(sigma_prior, 0.01))
  
  comp <- mags ~ field(day, model = spde) + depth_scaled + Intercept(1)
  
  fit <- tryCatch({
    bru(
      components = comp,
      data = df,
      family = "exponential",
      options = list(
        control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
      )
    )
  }, error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  # Unified curve prediction (at depth = 0)
  day_grid <- seq(min(df$day), max(df$day), by = 1)
  pred <- predict(
    fit, data.frame(day = day_grid, depth_scaled = 0),
    ~ exp(field + Intercept) / log(10), n.samples = 1000
  )
  
  # Segments (the same curve)
  t0_day <- as.numeric(difftime(t0, time0_global, units = "days"))
  pred$segment <- ifelse(pred$day < t0_day, "Before", "After")
  
  summarize_seg <- function(p) tibble(
    mean_b = mean(p$mean, na.rm = TRUE),
    ci_low = quantile(p$mean, 0.025, na.rm = TRUE),
    ci_high= quantile(p$mean, 0.975, na.rm = TRUE)
  )
  sum_before <- summarize_seg(dplyr::filter(pred, segment == "Before"))
  sum_after  <- summarize_seg(dplyr::filter(pred, segment == "After"))
  
  # Extract DIC/WAIC
  dic_val  <- if (!is.null(fit$dic$dic))   as.numeric(fit$dic$dic)   else NA_real_
  waic_val <- if (!is.null(fit$waic$waic)) as.numeric(fit$waic$waic) else NA_real_
  # Optional: Average -log CPO as an additional diagnosis
  lcpo_val <- if (!is.null(fit$cpo$cpo)) mean(-log(fit$cpo$cpo), na.rm = TRUE) else NA_real_
  
  list(
    fit = fit,
    pred = pred,
    summary = tibble(
      b_before = sum_before$mean_b, ci_before_low = sum_before$ci_low,  ci_before_high = sum_before$ci_high,
      b_after  = sum_after$mean_b,  ci_after_low  = sum_after$ci_low,   ci_after_high  = sum_after$ci_high,
      N_before = N_before_events,   N_after       = N_after_events,
      diff     = sum_after$mean_b - sum_before$mean_b,
      DIC = dic_val, WAIC = waic_val, LCPO = lcpo_val
    )
  )
}

# Encapsulation: Perform a fitting (including depth) on a single sequence once
analyze_sequence_once_with_depth <- function(seq, m0 = 2.5, window_days = 90,
                                             range_prior = 116.9, sigma_prior = 0.1) {
  t0  <- seq$mainshock_time
  dat <- seq$sequence_data %>%
    filter(time >= t0 - lubridate::days(window_days),
           time <= t0 + lubridate::days(window_days))
  res <- fit_sequence_once_with_depth(dat, t0, m0, range_prior, sigma_prior)
  if (is.null(res)) return(NULL)
  res$summary <- res$summary %>% mutate(mainshock_time = t0, mag = seq$mainshock_mag)
  res
}

# ==== Function: Plot one before/after b-value curve with label ====
plot_bvalue_segment_split <- function(seq, res_before, res_after) {
  pred_b <- res_before$pred
  pred_a <- res_after$pred
  pred_b$group <- "Before"
  pred_a$group <- "After"
  
  offset <- max(pred_b$day) + 10
  pred_a$day <- pred_a$day + offset
  plot_data <- bind_rows(pred_b, pred_a)
  
  label_text <- paste0(
    format(seq$mainshock_time, "%Y-%m-%d"), "\n",
    "M", round(seq$mainshock_mag, 1), "\n",
    "b_before = ", res_before$mean, " (N=", res_before$n, ")\n",
    "b_after  = ", res_after$mean,  " (N=", res_after$n, ")"
  )
  
  ggplot(plot_data, aes(x = day)) +
    geom_ribbon(aes(ymin = q0.025, ymax = q0.975, fill = group), alpha = 0.2) +
    geom_line(aes(y = mean, color = group), size = 1.2) +
    annotate("text", x = min(plot_data$day), y = 1.7, label = label_text,
             hjust = 0, vjust = 1, fontface = "italic", size = 3.2) +
    labs(x = "Days since segment start", y = "b-value") +
    coord_cartesian(ylim = c(0.3, 1.8)) +
    theme_minimal(base_size = 11) +
    scale_color_manual(values = c("Before" = "blue", "After" = "red")) +
    scale_fill_manual(values = c("Before" = "blue", "After" = "red")) +
    theme(legend.position = "none")
}

# ==== Main execution block ====
summary_table <- tibble()
plot_list <- list()

for (i in seq_along(sequence_list)) {
  cat("\n▶️ Seq", i, ":", format(sequence_list[[i]]$mainshock_time, "%Y-%m-%d"), "\n")
  
  res <- analyze_sequence_once_with_depth(sequence_list[[i]], m0 = 2.5, window_days = 90,
                                          range_prior = 116.9, sigma_prior = 0.1)
  if (is.null(res)) { cat("❌ Sequence skipped\n"); next }
  
  summary_table <- bind_rows(summary_table, res$summary)
  
  pred  <- res$pred
  sumry <- res$summary
  t0_day <- as.numeric(difftime(sequence_list[[i]]$mainshock_time, time0_global, units = "days"))
  
  # The relative time with the main shock as zero
  pred <- pred %>% mutate(x_rel = day - t0_day)
  pred_b <- filter(pred, segment == "Before")
  pred_a <- filter(pred, segment == "After")
  
  xmin <- min(pred$x_rel, na.rm = TRUE); xmax <- max(pred$x_rel, na.rm = TRUE)
  
  title_date <- format(sequence_list[[i]]$mainshock_time, "%Y-%m-%d")
  title_mag  <- sprintf("M%.1f", sequence_list[[i]]$mainshock_mag)
  label_text <- sprintf("%s\n%s\nb_before = %.3f (N=%d)\nb_after  = %.3f (N=%d)",
                        title_date, title_mag,
                        sumry$b_before, sumry$N_before,
                        sumry$b_after,  sumry$N_after)
  
  y_top <- max(pred$q0.975, na.rm = TRUE) + 0.05
  
  p <- ggplot() +
    # BEFORE（blue）
    geom_ribbon(data = pred_b, aes(x = x_rel, ymin = q0.025, ymax = q0.975),
                fill = "steelblue", alpha = 0.18) +
    geom_line(data = pred_b, aes(x = x_rel, y = mean),
              color = "steelblue", linewidth = 1) +
    geom_rect(aes(xmin = min(pred_b$x_rel, na.rm=TRUE),
                  xmax = max(pred_b$x_rel, na.rm=TRUE),
                  ymin = sumry$ci_before_low,
                  ymax = sumry$ci_before_high),
              fill = "steelblue", alpha = 0.15) +
    geom_segment(aes(x = min(pred_b$x_rel, na.rm=TRUE),
                     xend = max(pred_b$x_rel, na.rm=TRUE),
                     y = sumry$b_before, yend = sumry$b_before),
                 color = "steelblue", linewidth = 1.1, lineend = "round") +
    
    # AFTER（red）
    geom_ribbon(data = pred_a, aes(x = x_rel, ymin = q0.025, ymax = q0.975),
                fill = "red", alpha = 0.18) +
    geom_line(data = pred_a, aes(x = x_rel, y = mean),
              color = "red", linewidth = 1) +
    geom_rect(aes(xmin = min(pred_a$x_rel, na.rm=TRUE),
                  xmax = max(pred_a$x_rel, na.rm=TRUE),
                  ymin = sumry$ci_after_low,
                  ymax = sumry$ci_after_high),
              fill = "red", alpha = 0.15) +
    geom_segment(aes(x = min(pred_a$x_rel, na.rm=TRUE),
                     xend = max(pred_a$x_rel, na.rm=TRUE),
                     y = sumry$b_after, yend = sumry$b_after),
                 color = "red", linewidth = 1.1, lineend = "round") +
    
    # The main shock vertical line (x=0)
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    
    # Information in the upper left corner
    annotate("text", x = xmin + 0.02*(xmax - xmin), y = y_top,
             label = label_text, hjust = 0, vjust = 1, size = 3.6) +
    
    labs(x = "Days since sequence start", y = "b-value") +
    coord_cartesian(ylim = c(0.4, max(1.8, y_top))) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
  
  plot_list[[length(plot_list) + 1]] <- p
}

if (length(plot_list) > 0) {
  library(patchwork)
  plots_per_page <- 6
  n_pages <- ceiling(length(plot_list) / plots_per_page)
  for (page in seq_len(n_pages)) {
    start <- (page - 1) * plots_per_page + 1
    end   <- min(page * plots_per_page, length(plot_list))
    page_plot <- wrap_plots(plot_list[start:end], ncol = 2)  # 每页2列×3行
    print(page_plot)
    
    ggsave(sprintf("sequence_depth_page_%02d.png", page),
           page_plot, width = 12, height = 9, dpi = 300)
  }
} else {
  message("No plots generated.")
}
# Print summary
print(summary_table)

# summary_table with DIC/WAIC
dic_waic_table <- summary_table %>%
  select(mainshock_time, mag, N_before, N_after, DIC, WAIC, LCPO, 
         b_before, b_after, diff) %>%
  arrange(WAIC)   # or arrange(DIC)

print(dic_waic_table)


