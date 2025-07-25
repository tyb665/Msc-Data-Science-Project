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
combined_all <- wrap_plots(plot_list, ncol = 2)
print(combined_all)

# ==== Filter the main shocks after 1990 from the sequence_list ====
recent_indices <- which(sapply(sequence_list, function(seq) seq$mainshock_time >= as.POSIXct("1990-01-01")))
plot_list_recent <- plot_list[recent_indices]

# ==== Plot the sequence after 1990 ====
combined_recent <- wrap_plots(plot_list_recent, ncol = 2)
print(combined_recent)

##############################
##### Try different m0########
##############################
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
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(INLA)
library(inlabru)
library(patchwork)

# ==== Set environment ====
Sys.setlocale("LC_TIME", "English")
Sys.setenv(TMPDIR = "D:/R_temp")

# ==== Read combined data ====
file_path <- "D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv"
eq_data <- read_csv(file_path) %>%
  rename(mag = magnitude) %>%
  mutate(time = as.POSIXct(time, tz = "UTC")) %>%
  drop_na(time, mag, latitude, longitude)

# ==== Extract mainshocks (m ≥ 6.5, post-1990) ====
mainshocks_all <- eq_data %>%
  filter(mag >= 6.5, time >= as.POSIXct("1990-01-01")) %>%
  arrange(time)

# ==== De-duplicate overlapping mainshocks ====
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

# ==== Build sequence list ====
sequence_list <- list()
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
  
  if (nrow(seq_data) > 0) {
    sequence_list[[i]] <- list(
      mainshock_time = t0,
      mainshock_mag = mainshock$mag,
      sequence_data = seq_data
    )
  }
}

# ==== Modeling for m0 = 2.5, 3.5, 4.5 ====
m0_list <- c(2.5, 3.5, 4.5)
plots_by_m0 <- list()

for (m0 in m0_list) {
  cat("=== Modeling m0 =", m0, "===\n")
  plot_list <- list()
  
  for (i in seq_along(sequence_list)) {
    seq <- sequence_list[[i]]
    t0 <- seq$mainshock_time
    seq_data <- seq$sequence_data %>%
      filter(mag >= m0) %>%
      mutate(
        day = as.numeric(difftime(time, min(time), units = "days")),
        mags = mag - m0
      )
    
    if (nrow(seq_data) < 50) next  # Skip if data too sparse
    
    # Mesh
    seq_data <- seq_data[order(seq_data$day), ]
    seq_data$day <- seq_data$day + 0.001
    extremes <- range(seq_data$day)
    mesh_points <- seq(extremes[1], extremes[2], by = 5)
    mesh <- fm_mesh_1d(mesh_points, cutoff = 2, degree = 2, boundary = "free")
    
    # SPDE model
    spde_model <- inla.spde2.pcmatern(
      mesh,
      prior.range = c(20, 0.01),
      prior.sigma = c(0.5, 0.01)
    )
    
    # Fit
    comp <- mags ~ field(day, model = spde_model) + Intercept(1)
    
    fit <- tryCatch({
      bru(
        components = comp,
        data = seq_data,
        family = "exponential",
        options = list(control.compute = list(dic = TRUE, waic = TRUE))
      )
    }, error = function(e) NULL)
    
    if (is.null(fit)) next
    
    pred <- predict(fit, data.frame(day = seq(0, max(seq_data$day), by = 0.5)),
                    ~ exp(field + Intercept) / log(10),
                    n.samples = 1000)
    
    # Plot
    p <- ggplot(pred, aes(x = day)) +
      geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = "grey80") +
      geom_line(aes(y = mean), color = "black") +
      labs(
        title = paste0("Seq ", i, " | M", round(seq$mainshock_mag, 1),
                       " @ ", format(t0, "%Y-%m-%d")),
        x = "Days since mainshock", y = "b-value"
      ) +
      theme_minimal(base_size = 10) +
      coord_cartesian(ylim = c(0.3, 1.7))
    
    plot_list[[length(plot_list) + 1]] <- p
  }
  
  plots_by_m0[[as.character(m0)]] <- plot_list
}

# ==== Display each m0's results in pages of 6 plots ====
plots_per_page <- 6

for (m0 in names(plots_by_m0)) {
  plot_list <- plots_by_m0[[m0]]
  n_pages <- ceiling(length(plot_list) / plots_per_page)
  cat("\n>>> m0 =", m0, "with", length(plot_list), "plots (", n_pages, "pages)\n")
  
  for (page in seq_len(n_pages)) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, length(plot_list))
    page_plots <- wrap_plots(plot_list[start_idx:end_idx], ncol = 2)
    print(page_plots)
  }
}


#####
# Use the whole data of 1990-2021, m0 = 2.5, compare 3 strategies of mesh point
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
prior_range_list <- round(c(T2/100, T2/50, T2/20, T2/10, T2/5, T2/2, T2/1.5), 1)
prior_sigma_list <- c(0.1, 0.3, 1, 3, 10)

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
results_heatmap <- results %>%
  mutate(
    range = factor(range, levels = sort(unique(range))),
    sigma = factor(sigma, levels = sort(unique(sigma)))
  )

p_dic <- ggplot(results_heatmap, aes(x = range, y = sigma, fill = DIC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(DIC, 1)), size = 3, color = "white") +
  scale_fill_viridis_c(option = "inferno", direction = -1) +
  labs(title = "DIC across Prior Settings (Real Data, KDE mesh)",
       x = "Prior Range", y = "Prior Sigma", fill = "DIC") +
  theme_minimal()

p_waic <- ggplot(results_heatmap, aes(x = range, y = sigma, fill = WAIC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(WAIC, 1)), size = 3, color = "white") +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(title = "WAIC across Prior Settings (Real Data, KDE mesh)",
       x = "Prior Range", y = "Prior Sigma", fill = "WAIC") +
  theme_minimal()

# Display plots
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
#Function: Construct KDE-based mesh and fit b-value
fit_bvalue_segment <- function(df_segment, m0 = 2.5, range_prior = 116.9, sigma_prior = 0.1, min_n = 20) {
  
  # 1. Data preprocessing + filtering
  df_segment <- df_segment %>%
    filter(mag >= m0) %>%
    mutate(
      day = as.numeric(difftime(time, min(time), units = "days")),
      mags = mag - m0
    )
  
  if (nrow(df_segment) < min_n) {
    cat("❌ Not enough data: ", nrow(df_segment), "\n")
    return(NULL)
  }
  
  # 2. Mesh 
  kde <- density(df_segment$day, bw = 30)
  mesh_points <- sample(kde$x, size = min(100, length(kde$x)), prob = kde$y, replace = FALSE)
  mesh_points <- sort(unique(c(min(df_segment$day), mesh_points, max(df_segment$day))))
  mesh <- fm_mesh_1d(mesh_points, degree = 2, boundary = "free")
  
  # 3. SPDE 
  spde <- inla.spde2.pcmatern(
    mesh,
    prior.range = c(range_prior, 0.01),
    prior.sigma = c(sigma_prior, 0.01)
  )
  comp <- mags ~ field(day, model = spde) + Intercept(1)
  
  # 4. fitting
  fit <- tryCatch({
    bru(
      components = comp,
      data = df_segment,
      family = "exponential",
      options = list(control.compute = list(dic = TRUE, waic = TRUE))
    )
  }, error = function(e) {
    cat("❌ Model fitting failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(fit)) {
    cat("❌ Fit is NULL\n")
    return(NULL)
  }
  
  # 5. Prediction + output
  pred <- predict(fit, data.frame(day = seq(0, max(df_segment$day), by = 1)),
                  ~ exp(field + Intercept) / log(10), n.samples = 1000)
  
  mean_b <- round(mean(pred$mean), 3)
  low_b <- round(mean(pred$q0.025), 3)
  high_b <- round(mean(pred$q0.975), 3)
  
  cat("✅ Fit success | b =", mean_b, " [", low_b, ",", high_b, "] | n =", nrow(df_segment), "\n")
  
  list(pred = pred, mean = mean_b, ci = c(low_b, high_b), n = nrow(df_segment))
}


# ==== Main loop: Model all sequences before and after execution and output graphs and tables ====
results_list <- list()
summary_table <- tibble()
plot_list <- list()

for (i in seq_along(sequence_list)) {
  seq <- sequence_list[[i]]
  t0 <- seq$mainshock_time
  cat("\n==============================\n")
  cat("▶️ Processing Sequence", i, "| Mainshock:", format(t0, "%Y-%m-%d %H:%M:%S"), "\n")
  
  res <- analyze_sequence_pair(seq, m0 = 2.5, window_days = 90)
  
  if (!is.null(res)) {
    cat("✅ Sequence", i, "successful\n")
    summary_table <- bind_rows(summary_table, res$summary)
    plot_list[[length(plot_list) + 1]] <- res$plot
  } else {
    cat("❌ Sequence", i, "skipped (either before/after failed or insufficient data)\n")
  }
}

# ==== Print the number of samples before and after each main shock ====
cat("\n========== Sample Sizes by Sequence ==========\n")
for (i in seq_along(sequence_list)) {
  seq <- sequence_list[[i]]
  t0 <- seq$mainshock_time
  seq_data <- seq$sequence_data
  
  df_before <- seq_data %>% filter(time >= t0 - days(90), time < t0)
  df_after  <- seq_data %>% filter(time > t0, time <= t0 + days(90))
  
  cat(sprintf("🟢 %s | N_before = %d | N_after = %d\n",
              format(t0, "%Y-%m-%d"),
              nrow(df_before), nrow(df_after)))
}

# ==== Print the result table and draw graphs ====
cat("\n========== Summary Table ==========\n")
print(summary_table)

wrap_plots(plot_list, ncol = 2)



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

# ==== Summary trend statistics ====
## !!We only have 6 sequences!!
cat("\n===== Summary of b-value Differences (After - Before) =====\n")

mean_diff <- mean(summary_table$diff)
median_diff <- median(summary_table$diff)
t_test <- t.test(summary_table$diff)

cat(sprintf("Mean Δb-value: %.3f\n", mean_diff))
cat(sprintf("Median Δb-value: %.3f\n", median_diff))
cat("One-sample t-test against 0:\n")
print(t_test)
