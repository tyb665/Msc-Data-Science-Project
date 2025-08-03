

#####################
#Plot for background###
#####################


#######
#Real data GR-law example
# 加载库
library(ggplot2)
library(readr)
library(dplyr)

# 读取真实数据
data <- read_csv("D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv")

# 设置震级下限（M0）
M0 <- 2.5

# 选取非缺失并满足 M >= M0 的震级
mags_values <- data %>% 
  filter(!is.na(magnitude), magnitude >= M0) %>% 
  pull(magnitude)

# 创建震级分箱（bins）
mags_bins <- seq(M0, max(mags_values), by = 0.1)

# 计算累计频率 N(M)
mags_counts <- vapply(mags_bins, \(bin) sum(mags_values >= bin), 0)

# 计算 log10 N(M)
log_counts <- log10(mags_counts)

# 拟合线性模型
fit <- lm(log_counts ~ mags_bins)
slope <- coef(fit)[2]
b_value_est <- abs(round(slope, 3))
line_value <- predict(fit)

# 绘图
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

# 读取数据
data <- read.csv("D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv")

# 预处理
data_clean <- data %>%
  filter(magnitude >= 2.5) %>%
  mutate(
    time = ymd_hms(time),
    time_numeric = year(time) + yday(time) / 365.25  # 连续时间变量
  )

# 构建一维时间 mesh
time_range <- range(data_clean$time_numeric)
time_mesh <- inla.mesh.1d(
  seq(time_range[1], time_range[2], length.out = 100)
)

# SPDE 模型
spde_time <- inla.spde2.pcmatern(
  mesh = time_mesh,
  prior.range = c(0.5, 0.9),  # 平滑程度：P(range < 0.5) = 0.1
  prior.sigma = c(1, 0.01)    # 振幅控制
)

# 定义组件
cmp <- ~ Intercept(1, model = "linear") +
  beta_field(time_numeric, model = spde_time)


# 自定义对数似然函数
loglike <- function(beta, data, M0 = 2.5) {
  sum(log(beta) - beta * (data$magnitude - M0))
}

# 拟合公式（事件的数量看作 point process）
bru_formula <- magnitude ~ Intercept + beta_field

# 构建 inlabru likelihood
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



# 构造 prediction time points
pred_time <- data.frame(time_numeric = seq(time_range[1], time_range[2], length.out = 300))

# 预测 log β(t)
pred <- predict(fit, pred_time, ~ Intercept + beta_field, n.samples = 1000)

library(ggplot2)

ggplot(pred, aes(x = time_numeric)) +
  geom_line(aes(y = mean), color = "blue") +
  geom_ribbon(aes(ymin = q0.025, ymax = q0.975), alpha = 0.2, fill = "skyblue") +
  labs(
    title = expression(Temporal~evolution~of~log~beta(t)),
    x = "Time (year)",
    y = expression(log~beta(t))
  ) +
  xlim(1920, 2023) +
  theme_minimal()


##############################
#####Gaussian process beta(s)
##############################
library(dplyr)
library(INLA)
library(inlabru)
library(ggplot2)
library(maps)

# 读取数据
data <- read.csv("D:/project data/MSc_Data_Science_Project/combined_earthquakes_cleaned.csv")

# 筛选可靠数据：M >= 2.5 且经纬度不为空
data_clean <- data %>%
  filter(magnitude >= 2.5, !is.na(latitude), !is.na(longitude)) %>%
  mutate(x = longitude, y = latitude)

# 生成空间 mesh（基于震中）
coords <- cbind(data_clean$x, data_clean$y)

mesh <- inla.mesh.2d(
  loc = coords,
  max.edge = c(0.5, 5),  # 视数据稠密程度而定，可调整
  cutoff = 0.2,
  offset = c(1, 2)
)

# 构建 SPDE model
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(1, 0.01),    # P(range < 1 deg) = 0.01
  prior.sigma = c(1, 0.01)     # P(sigma > 1) = 0.01
)

# 组件命名为 "spatial_field"
cmp <- ~ Intercept(1, model = "linear") +
  spatial_field(coordinates, model = spde)

# 公式
formula <- magnitude ~ Intercept + spatial_field

# 定义坐标数据（用于匹配 SPDE 的空间索引）
data_clean$coordinates <- coords

# 拟合：对 magnitude 建模（作为 proxy for log β(s)）
fit <- bru(
  formula = formula,
  components = cmp,
  data = data_clean,
  family = "exponential",  # 或其它近似分布，模拟 GR 模型对 β 的推断
  control.family = list(control.link = list(model = "log")),
  options = list(verbose = TRUE)
)

# 创建覆盖所有观测点的经纬度网格
grid_res <- 100
x_range <- range(data_clean$x)
y_range <- range(data_clean$y)

pred_grid <- expand.grid(
  x = seq(x_range[1], x_range[2], length.out = grid_res),
  y = seq(y_range[1], y_range[2], length.out = grid_res)
)

pred_grid$coordinates <- cbind(pred_grid$x, pred_grid$y)

# 预测 log β(s)
pred <- predict(fit, pred_grid, ~ Intercept + spatial_field, n.samples = 1000)

# 可视化：热图 + 地震点
# 获取美国地图数据
usa_map <- map_data("state")

# 仅保留 California
california_map <- subset(usa_map, region == "california")
ggplot() +
  # 热图
  geom_tile(data = pred, aes(x = x, y = y, fill = mean)) +
  
  # 地震点（点型图例）
  geom_point(
    data = data_clean,
    aes(x = x, y = y, shape = "Earthquake Epicenter"),
    color = "black", alpha = 0.3, size = 0.5
  ) +
  
  # 加州边界线（线型图例）
  geom_path(
    data = california_map,
    aes(x = long, y = lat, group = group, linetype = "California Border"),
    color = "black", size = 0.5
  ) +
  
  # 色标图例
  scale_fill_viridis_c(option = "plasma", name = expression(log~beta(s))) +
  
  # 设置 shape 图例
  scale_shape_manual(name = "", values = c("Earthquake Epicenter" = 16)) +
  
  # 设置 linetype 图例
  scale_linetype_manual(name = "", values = c("California Border" = "solid")) +
  
  # 明确将 shape 和 linetype 图例分离管理
  guides(
    shape = guide_legend(override.aes = list(size = 2, alpha = 1)),
    linetype = guide_legend(override.aes = list(size = 1))
  ) +
  
  labs(
    title = expression("Posterior mean of " * log~beta(s)),
    x = "Longitude", y = "Latitude"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

