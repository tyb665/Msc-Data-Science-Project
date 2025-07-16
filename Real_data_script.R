# Isolate different sequences (m≥6.5) 
# ==== Load libraries ====
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)
Sys.setlocale("LC_TIME", "English")#Manually set the language environment

# ==== Read and clean the data ====
file_path <- "D:/project data/MSc_Data_Science_Project/earthquakes.csv"
eq_data <- read_csv(file_path)

eq_data <- eq_data %>%
  rename(mag = magnitude) %>%
  mutate(time = as.POSIXct(time, tz = "UTC")) %>%
  drop_na(time, mag, latitude, longitude)

# ==== Extract the main shocks with m≥6.5 ====
mainshocks <- eq_data %>%
  filter(mag >= 6.5) %>%
  arrange(time)

# ==== Set the windows ====
window_days <- 60
lat_window <- 1.5
lon_window <- 1.5

# ==== Extract the sequence and generate the graph ====
sequence_list <- list()
plot_list <- list()

for (i in seq_len(nrow(mainshocks))) {
  mainshock <- mainshocks[i, ]
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
  
  sequence_list[[i]] <- list(
    mainshock_time = t0,
    mainshock_mag = mainshock$mag,
    sequence_data = seq_data
  )
  
  # Subplot: Magnitude vs. Time
  p <- ggplot(seq_data, aes(x = time, y = mag)) +
    geom_point(color = "orange", alpha = 0.7) +
    geom_vline(xintercept = as.numeric(t0), color = "red", linetype = "dashed") +
    labs(
      title = paste0("Sequence ", i, " | M", round(mainshock$mag, 1),
                     " @ ", format(t0, "%Y-%m-%d")),
      x = "Time", y = "Magnitude"
    ) +
    theme_minimal(base_size = 10)
  
  plot_list[[i]] <- p
}

# ==== Combine all subplots ====
combined_plot <- wrap_plots(plot_list, ncol = 2)  

# ==== Show combined plot ====
print(combined_plot)


## Use the whole data 2008-2021, modify the mesh points according to density 
# ==== Library ====
library(readr)       
library(dplyr)        
library(tidyr)        
library(lubridate)    
library(ggplot2)      
library(patchwork)    
library(INLA)        
library(inlabru)      

# ==== Set the English time display ====
Sys.setlocale("LC_TIME", "English")

# ==== Read and clean the data ====
file_path <- "D:/project data/MSc_Data_Science_Project/earthquakes.csv"
eq_data <- read_csv(file_path)

eq_data <- eq_data %>%
  rename(mag = magnitude) %>%
  mutate(time = as.POSIXct(time, tz = "UTC")) %>%
  drop_na(time, mag, latitude, longitude)

# ==== Select time（2008–2021）====
eq_filtered <- eq_data %>%
  filter(time >= as.POSIXct("2008-01-01"),
         time <= as.POSIXct("2021-12-31"))

# ==== Convert the time to "days" ====
start_time <- min(eq_filtered$time)
eq_filtered <- eq_filtered %>%
  mutate(day = as.numeric(difftime(time, start_time, units = "days")),
         mags = mag - 2.5)  # Convert to residual magnitude


# ==== Create mesh（time） ====
mesh1D <- fm_mesh_1d(eq_filtered$day, cutoff = 30, degree = 2, boundary = "free")

# ==== SPDE ====
spde_model <- inla.spde2.pcmatern(
  mesh1D,
  prior.range = c(200, 0.01),     # Change to optimal combination later
  prior.sigma = c(0.5, 0.01)
)

# ==== Fitting ====
comp <- mags ~ field(day, model = spde_model) + Intercept(1)

fit <- bru(
  components = comp,
  data = eq_filtered,
  family = "exponential",
  options = list(control.compute = list(dic = TRUE, waic = TRUE))
)

# ==== Prediction & Plotting ====
time_pred <- data.frame(day = seq(0, max(eq_filtered$day), by = 20))

pred <- predict(fit, time_pred, ~ exp(field + Intercept) / log(10), n.samples = 1000)

ggplot() +
  geom_ribbon(data = pred, aes(x = day, ymin = q0.025, ymax = q0.975), fill = "grey80") +
  geom_line(data = pred, aes(x = day, y = mean), color = "black") +
  labs(
    title = "b-value over time (2008–2021)",
    x = "Days since 2008-01-01", y = "b-value"
  ) +
  theme_minimal()

## It seems the model did not converge, so we try different m0

#######
install.packages("INLA", repos = c(getOption("repos"),
                                   INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
system.file("bin/windows/64bit/inla.exe", package = "INLA")
#######


## Use data 2008-2021, try different m0
# ==== Library ====
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(INLA)
options(inla.call = "remote") 
library(inlabru)

# ==== Set the English month ====
# ==== Set environment ====
Sys.setlocale("LC_TIME", "English")                      
Sys.setenv(TMPDIR = "D:/R_temp")                         # Avoid INLA temporary file issues

# ==== Read and preprocess the data ====
file_path <- "D:/project data/MSc_Data_Science_Project/earthquakes.csv"
eq_data <- read_csv(file_path) |> 
  rename(mag = magnitude) |> 
  mutate(time = as.POSIXct(time, tz = "UTC")) |> 
  drop_na(time, mag, latitude, longitude) |> 
  filter(time >= as.POSIXct("2008-01-01"), time <= as.POSIXct("2021-12-31"))

# ==== list of m0 thresholds ====
m0_list <- c(2.5)#, 2.8, 3.0, 3.2)
fit_results <- list()

# ==== Loop modeling different m0 ====
for (m0 in m0_list) {
  cat(">>> Fitting model with m0 =", m0, "\n")
  
  eq_m <- eq_data %>%
    filter(mag > m0) %>%
    mutate(
      day = as.numeric(difftime(time, min(time), units = "days")),
      mags = mag - m0
    )
  eq_m<- as.data.frame(eq_m)
  eq_m$day =eq_m$day + 0.001
  eq_m <- eq_m[order(eq_m$day),]
  # mesh_points <- seq(start_day, end_day, by = 100)
  extremes = range(eq_m$day)
  mesh_points <- seq(extremes[1], extremes[2], by = 100)
  mesh <- fm_mesh_1d(mesh_points, cutoff = 30, degree = 2, boundary = "free")
  
  spde_model <- inla.spde2.pcmatern(
    mesh,
    prior.range = c(200, 0.01),## Change to the optimal combination later
    prior.sigma = c(0.5, 0.01)
  )
  
  comp <- mags ~ field(day, model = spde_model) + Intercept(1)
  
  fit <- tryCatch({
    bru(
      components = comp,
      data = eq_m,
      family = "exponential",
      options = list(control.compute = list(dic = TRUE, waic = TRUE, return.marginals = TRUE))
    )
  }, error = function(e) {
    message("Failed for m0 = ", m0, " | ", e$message)
    return(NULL)
  })
  
  if (!is.null(fit)) {
    pred_df <- predict(fit, data.frame(day = seq(0, max(eq_m$day), by = 20)),
                       ~ exp(field + Intercept) / log(10),
                       n.samples = 1000)
    
    fit_results[[as.character(m0)]] <- list(
      m0 = m0,
      dic = fit$dic$dic,
      waic = fit$waic$waic,
      pred = pred_df
    )
  }
}

# ==== Summary table ====
summary_table <- bind_rows(lapply(fit_results, function(res) {
  tibble(
    m0 = res$m0,
    DIC = round(res$dic, 2),
    WAIC = round(res$waic, 2)
  )
})) |> arrange(DIC)

print(summary_table)

# ==== Visualize the prediction curves of different M0 ====
plots <- list()
for (i in seq_along(fit_results)) {
  m0 <- names(fit_results)[i]
  pred <- fit_results[[i]]$pred
  p <- ggplot(pred, aes(x = day)) +
    geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = "grey80") +
    geom_line(aes(y = mean), color = "black", size = 1) +
    labs(title = paste("m0 =", m0), y = "b-value", x = "days since 2008-01-01") +
    ylim(0.3, 1.7)
  plots[[m0]] <- p
}

patchwork::wrap_plots(plots, ncol = 2)



