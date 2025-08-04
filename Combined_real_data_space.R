########## Space modeling
library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)

# ==== 读取数据 ====
file_path <- "D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv"

eq_data <- read_csv(file_path) %>%
  rename(mag = magnitude) %>%
  drop_na(time, mag, latitude, longitude) %>%
  mutate(
    time = as.POSIXct(time, tz = "UTC"),
    year = year(time)
  ) %>%
  filter(year >= 1990)

# ==== 筛选震级 ≥ m0 ====
m0 <- 2.5
eq_filtered <- eq_data %>%
  filter(mag >= m0)

# ==== 标准化空间坐标（用于SPDE）====
eq_filtered <- eq_filtered %>%
  mutate(
    x = scale(longitude)[, 1],
    y = scale(latitude)[, 1],
    mags = mag - m0
  )

# ==== 可视化震中分布 ====
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

# ==== 提取空间坐标 ====
eq_filtered <- eq_filtered %>%
  mutate(
    x = scale(longitude)[, 1],
    y = scale(latitude)[, 1],
    mags = mag - m0
  )

coords <- dplyr::select(eq_filtered, x, y)

# ==== KDE 估计空间密度 ====
set.seed(42)
kde <- MASS::kde2d(coords$x, coords$y, n = 200)

# ==== 根据KDE结果进行概率采样点 ====
# 转为向量
x_vec <- rep(kde$x, each = length(kde$y))
y_vec <- rep(kde$y, times = length(kde$x))
z_vec <- as.vector(kde$z)

# 去除极低密度区域（避免偏远点）
threshold <- quantile(z_vec, 0.05)
valid_idx <- which(z_vec > threshold)

sample_pool <- data.frame(x = x_vec[valid_idx], y = y_vec[valid_idx], weight = z_vec[valid_idx])
sample_pool <- sample_pool[!is.na(sample_pool$weight), ]

# 采样mesh点（建议1000点以内）
mesh_points <- sample_n(sample_pool, size = 600, weight = weight, replace = FALSE)

# 添加边界点（范围控制）
boundary_x <- range(coords$x)
boundary_y <- range(coords$y)
boundary_buffer <- 0.2

boundary <- inla.nonconvex.hull(as.matrix(coords), convex = -0.05)  # 非凸边界更符合加州形状

# ==== 构建三角网格 ====
mesh <- inla.mesh.2d(
  loc = mesh_points[, c("x", "y")],
  boundary = boundary,
  max.edge = c(0.3, 1),     # 控制三角边长度（更精细）
  cutoff = 0.05             # 限制点之间最小距离，避免过密
)

# ==== 可视化 Mesh ====
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

# ==== 0. 标准化坐标 ====
eq_filtered <- eq_filtered %>%
  mutate(
    mags = mag - m0,
    x = scale(longitude)[,1],
    y = scale(latitude)[,1]
  )

coords <- dplyr::select(eq_filtered, x, y)

# ==== 1. KDE估计密度 + 轮廓边界 ====
kde <- kde2d(eq_filtered$x, eq_filtered$y, n = 200)
contour_level <- quantile(kde$z, 0.3)  # 可调密度阈值
contour_lines <- contourLines(kde$x, kde$y, kde$z, levels = contour_level)
boundary_coords <- cbind(contour_lines[[1]]$x, contour_lines[[1]]$y)

boundary <- inla.nonconvex.hull(boundary_coords, convex = -0.03)

# ==== 2. 构建 Mesh ====
# ==== 1. KDE估计密度 ====
library(MASS)
kde <- kde2d(eq_filtered$x, eq_filtered$y, n = 200)

# ==== 2. 构建 KDE-based 采样点池 ====
x_vec <- rep(kde$x, each = length(kde$y))
y_vec <- rep(kde$y, times = length(kde$x))
z_vec <- as.vector(kde$z)

# 构建 dataframe 并去除低密度点
sample_pool <- data.frame(x = x_vec, y = y_vec, weight = z_vec) %>%
  filter(weight > quantile(weight, 0.05))  # 可调阈值

# 采样 mesh 点（密度越大权重越高）
set.seed(42)
mesh_points <- sample_n(sample_pool, size = 800, weight = weight)

# ==== 3. 构建非凸边界（整个震中区域） ====
boundary <- inla.nonconvex.hull(as.matrix(dplyr::select(eq_filtered, x, y)), convex = -0.03)

# ==== 4. 构建 mesh ====
mesh <- inla.mesh.2d(
  loc = mesh_points[, c("x", "y")],
  boundary = boundary,
  max.edge = c(0.2, 1),  # 可调（越小越细）
  cutoff = 0.05
)

# ==== 5. 可视化 ====
#plot(mesh, main = "KDE-weighted Mesh")
#points(eq_filtered$x, eq_filtered$y, col = rgb(0, 0, 0, 0.1), pch = 16, cex = 0.3)

# ==== 3. 设置 Prior 网格 ====
x_range <- diff(range(eq_filtered$x))
y_range <- diff(range(eq_filtered$y))
S <- round((x_range + y_range) / 2, 2)
cat("Estimated standardized space span S ≈", S, "\n")

prior_range_list <- round(c(S/50, S/20, S/10, S/5, S/2), 2)
prior_sigma_list <- c(0.1, 0.3, 1, 3, 5, 10)

# ==== 4. 构建数据（包含空间坐标） ====
df_bru_spatial <- eq_filtered %>%
  dplyr::select(mags, x, y) %>%
  mutate(coord = as.matrix(cbind(x, y)))

# ==== 5. 初始化结果 ====
results <- tibble(
  prior_range = numeric(),
  prior_sigma = numeric(),
  DIC = numeric(),
  WAIC = numeric(),
  range_mean = numeric(),
  stdev_mean = numeric()
)

# ==== 6. 网格搜索 ====
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

# ==== 7. 可视化 WAIC 热图 ====
results_plot <- results %>%
  mutate(
    prior_range = factor(prior_range, levels = sort(unique(prior_range))),
    prior_sigma = factor(prior_sigma, levels = sort(unique(prior_sigma)))
  )

best_dic_row <- results %>% filter(DIC == min(DIC, na.rm = TRUE))
best_waic_row <- results %>% filter(WAIC == min(WAIC, na.rm = TRUE))
# DIC 热图
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

# === Stdev for field (posterior mean) 热图 ===
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

# ==== 8. 输出最优参数 ====
print(results %>% arrange(WAIC))
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
# 正确构建 coord（69961 x 2 matrix）
eq_filtered$coord <- cbind(eq_filtered$x, eq_filtered$y)

# ✅ 定义组件：独立 Intercept，明确命名 field
components <- ~ Intercept(1) + field_spatial(coord, model = spde)

# ✅ likelihood 只引用 latent component 名，不写错字段名
likelihood <- like(
  formula = mags ~ Intercept + field_spatial,
  family = "exponential",
  data = eq_filtered
)

# ✅ 拟合模型（无 mean.linear）
fit <- bru(
  components = components,
  likelihood = likelihood,
  options = list(
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
    control.inla = list(int.strategy = "eb")
  )
)

# 查看模型摘要
summary(fit)

# ==== 5. 构建预测网格（KDE非均匀采样） ====
# 从 sample_pool 中按 KDE 权重采样预测点
# 采样预测点
pred_pool <- sample_n(sample_pool, size = 5000, weight = weight)

# 构造空间矩阵
coords_pred <- as.matrix(pred_pool[, c("x", "y")])

# 执行预测
pred <- predict(fit,
                formula = ~ exp(field_spatial + Intercept) / log(10),
                newdata = list(coord = coords_pred),
                n.samples = 1000)


# ==== 7. 可视化：b-value map ====
library(ggplot2)
library(viridis)

# 添加坐标信息
pred_df <- data.frame(
  x = coords_pred[, 1],
  y = coords_pred[, 2],
  b_mean = pred$mean,
  b_sd = pred$sd,
  b_lower = pred$q0.025,
  b_upper = pred$q0.975
)

# b-value 平均值图
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















