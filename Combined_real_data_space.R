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


# ==== 7. visual：b-value map（反标准化 + 加州边界）====
library(ggplot2)
library(viridis)
library(maps)
library(dplyr)

# —— 1) 反标准化回经纬度 —— 
# 注意：这里用 eq_filtered 中原始经纬度的均值和标准差来反标准化
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

# —— 2) 加州边界数据（行政边界/海岸线）——
usa_map <- map_data("state")
california_map <- subset(usa_map, region == "california")

# —— 3) 画图：b值后验均值 —— 
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

# —— 4) 画图：后验标准差（不确定性）——
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

# —— 5) 画图：95%置信区间宽度 —— 
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

# 打印或按需拼图
print(p_mean)
print(p_sd)
print(p_ci)



################################################################################
#Consider the covariates depths########################################
################################################################################
# ==== 0. 加载库 ====
library(INLA)
library(inlabru)
library(dplyr)
library(tidyr)

# ==== 1. 读取和预处理数据 ====
eq_data <- read.csv("D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv")

eq_filtered <- eq_data %>%
  filter(!is.na(magnitude), magnitude >= 2.5,
         !is.na(longitude), !is.na(latitude), !is.na(depth)) %>%
  mutate(
    mags        = magnitude - 2.5,                           # 震级减去 M0，作为响应变量
    x           = as.numeric(scale(longitude)[,1]),          # 经度标准化
    y           = as.numeric(scale(latitude)[,1]),           # 纬度标准化
    depth_s     = as.numeric(scale(depth)[,1]),              # 深度标准化
    coordinates = cbind(x, y)                                # 坐标矩阵
  )

# ==== 2. 构建 Mesh ====
library(MASS)

# Step 1: KDE 估计密度
kde <- kde2d(eq_filtered$x, eq_filtered$y, n = 200)

# Step 2: 构建采样池（剔除低密度区域）
x_vec <- rep(kde$x, each = length(kde$y))
y_vec <- rep(kde$y, times = length(kde$x))
z_vec <- as.vector(kde$z)

sample_pool <- data.frame(x = x_vec, y = y_vec, weight = z_vec) %>%
  filter(weight > quantile(weight, 0.05))   # 可调阈值

# Step 3: 加权采样 mesh 点
set.seed(2024)
mesh_points <- sample_n(sample_pool, size = 800, weight = weight)

# Step 4: 非凸边界（保持加州轮廓）
boundary <- inla.nonconvex.hull(as.matrix(eq_filtered[, c("x", "y")]), convex = -0.03)

# Step 5: 构建 mesh
mesh <- inla.mesh.2d(
  loc = mesh_points[, c("x", "y")],
  boundary = boundary,
  max.edge = c(0.2, 1),
  cutoff = 0.05
)


# ==== 3. 定义两个 SPDE 模型 ====
# 已有 mesh 的前提下
spde0 <- inla.spde2.pcmatern(mesh, prior.range = c(0.26, 0.01), prior.sigma = c(0.30, 0.01))
spde1 <- inla.spde2.pcmatern(mesh, prior.range = c(0.26, 0.01), prior.sigma = c(0.20, 0.01))

# 构建投影矩阵
A0 <- inla.spde.make.A(mesh, loc = eq_filtered$coordinates)
A1 <- inla.spde.make.A(mesh, loc = eq_filtered$coordinates, weights = eq_filtered$depth_s)


# ==== 4. 定义模型组件（mapper 自动实现变系数）====
components <- ~ Intercept(1) +
  depth_lin(depth_s, model = "linear") +
  field0(coordinates, model = spde0) +
  field1(coordinates, model = spde1)

# ==== 5. 定义似然函数 ====
likelihood <- like(
  formula = mags ~ Intercept + depth_lin + field0 + field1,
  family  = "exponential",
  data    = eq_filtered
)

# ==== 6. 拟合模型 ====
fit_vc <- bru(
  components = components,
  likelihood = likelihood,
  options = list(
    bru_initialisation = list(A = list(A0, A1)),  # ✅ 正确方式：显式传入两个场的投影矩阵
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
    control.inla = list(int.strategy = "eb")
  )
)

# ==== 7. 拟合结果摘要 ====
summary(fit_vc)

##Prediction and visualization
# ===============================
# 预测 + 可视化（含变系数与深浅对比）——修正版
# ===============================
library(INLA)
library(inlabru)
library(dplyr)
library(ggplot2)
library(viridis)
library(maps)

# ---- 0) 反标准化经纬度 + 底图 ----
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

# ---- 1) 预测点：优先用 KDE 样本池，否则规则网格 ----
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

# ---- 2) 生成 A 矩阵（注意 field1 需要 depth 权重） ----
make_pred_As <- function(coords_pred, depth_s_vec) {
  A0_pred <- INLA::inla.spde.make.A(mesh, loc = coords_pred)                         # field0
  A1_pred <- INLA::inla.spde.make.A(mesh, loc = coords_pred, weights = depth_s_vec) # field1 * depth
  list(A0_pred = A0_pred, A1_pred = A1_pred)
}

# ---- 3) 统一预测函数（关键修改：newdata 里变量名必须叫 depth_s；公式用 '+ depth_lin'） ----
predict_b <- function(depth_s_vec, nsamp = 1000) {
  stopifnot(length(depth_s_vec) == nrow(coords_pred))
  As <- make_pred_As(coords_pred, depth_s_vec)
  pred <- predict(
    fit_vc,
    # depth_lin 会自动用 newdata$depth_s 乘系数；field1 的乘法由 A1 的 weights 完成
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

# ---- 4) 情景 A：代表性深度（中位数） ----
depth_rep <- median(eq_filtered$depth_s, na.rm = TRUE)
pred_rep  <- predict_b(rep(depth_rep, nrow(coords_pred)), nsamp = 1000) %>% lonlat_back()
pred_rep$ci_width <- pred_rep$b_upper - pred_rep$b_lower

p_rep_mean <- plot_raster(pred_rep, "b_mean",  "Posterior Mean of b (depth = median)", "b-value", "plasma")
p_rep_sd   <- plot_raster(pred_rep, "b_sd",    "Posterior SD of b (depth = median)",   "Posterior SD", "magma")
p_rep_ci   <- plot_raster(pred_rep, "ci_width","95% CI Width of b (depth = median)",   "95% CI Width", "inferno")

print(p_rep_mean); print(p_rep_sd); print(p_rep_ci)

# ---- 5) 情景 B：浅层 vs 深层 对比 + Δb ----
q_shallow <- as.numeric(quantile(eq_filtered$depth_s, 0.25, na.rm = TRUE))
q_deep    <- as.numeric(quantile(eq_filtered$depth_s, 0.75, na.rm = TRUE))  # 可改 0.90

pred_shallow <- predict_b(rep(q_shallow, nrow(coords_pred)), nsamp = 1000) %>% lonlat_back()
pred_deep    <- predict_b(rep(q_deep,    nrow(coords_pred)), nsamp = 1000) %>% lonlat_back()

p_shallow <- plot_raster(pred_shallow, "b_mean",
                         sprintf("Mean b at shallow depth (Q25=%.2f)", q_shallow),
                         "b-value", "plasma")
p_deepmap <- plot_raster(pred_deep, "b_mean",
                         sprintf("Mean b at deep depth (Q75=%.2f)", q_deep),
                         "b-value", "plasma")
print(p_shallow); print(p_deepmap)

# Δb：深层 − 浅层
delta <- left_join(
  pred_deep, pred_shallow,
  by = c("x", "y", "longitude", "latitude"),
  suffix = c("_deep", "_shallow")
) %>%
  mutate(
    db_mean = b_mean_deep - b_mean_shallow,
    # 置信区间显著性判断
    sig_neg = (b_upper_deep < b_lower_shallow),
    sig_pos = (b_lower_deep > b_upper_shallow)
  )


p_db <- plot_raster(delta, "db_mean",
                    "Δb = b(deep) − b(shallow)", "Δb", "cividis")
print(p_db)

# ---- 6) 情景 C：深度效应强度（空间变系数） map ----
# 让 depth_s = 1，则 depth_lin + field1 就是 log β 对 depth_s 的单位斜率
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

##Spatial model vs. spatial + depth model