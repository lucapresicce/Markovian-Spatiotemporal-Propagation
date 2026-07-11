################################################################################
# COPERNICUS EXPLORATORY DATA ANALYSIS
################################################################################

# packages ----------------------------------------------------------------

library(sf)
library(sp)
library(akima)
library(tidyr)
library(dplyr)
library(purrr)
library(Matrix)
library(slider)
library(viridis)
library(ggplot2)
library(ggplotify)
library(patchwork)
library(lubridate)
library(data.table)
library(rnaturalearth)


# data loading ------------------------------------------------------------

load("data/copernicus_data.RData")
load("data/copernicus_predictors.RData")
lon <- unique(coords[,1])
lat <- unique(coords[,2])

n <- dim(spatial_data)[1]
q <- dim(spatial_data)[2]
p <- dim(spatial_predictors)[2]
t <- dim(spatial_data)[3]

# time subsetting
spatial_data <- spatial_data[, , (t-144+1):t]
spatial_predictors <- spatial_predictors[, , (t-144+1):t]
time_points <- seq(from = as.Date("1850-01-01"), by = "month", length.out = t)[(t-144+1):t] # "2003-01-01"
t <- dim(spatial_data)[3]


# DYNBPS variogram informed grid ------------------------------------------

library(gstat)
library(spacetime)

# using only training data to inform grid
set.seed(42)
idx <- sample(1:n, 600, FALSE)

results <- list()
initial_model <- vgmST(
  stModel = "separable",
  space = vgm(psill = 1, model = "Exp", nugget = 0.1, range = 10),
  time  = vgm(psill = 1, model = "Exp", nugget = 0.1, range = 5),
  sill  = 1)
for (i in 1:q) {
  data_i <- spatial_data[idx, i, ]

  predictor_i <- spatial_predictors[idx,,]
  eps_i <- matrix(0, length(idx), t)
  
  for (j in 1:t) {
    P_t <- scale(predictor_i[,,j])
    L_t <- matrix(0, length(idx), 11)
    if (format(time_points, "%m")[j] != "12") {
      L_t[,as.numeric(format(time_points, "%m"))[j]] <- 1
    }
    eps_i[,j] <- residuals(lm(data_i[,j] ~ cbind(P_t, L_t)))
  }
  
  data_vec <- as.vector(eps_i)
  
  stfdf <- STFDF(sp = SpatialPoints(coords[idx,]),
                 time = as.POSIXct(1:dim(data_i)[2], origin = "2002-01-01"),
                 data = data.frame(value = data_vec))
  
  vgram <- variogramST(value ~ 1, data = stfdf, tlags = seq(0, 120, 40), cores = 12)
  
  fitted_model <- fit.StVariogram(object = vgram, model = initial_model, method = "L-BFGS-B",
                                  lower = c(0.0001, 0.0001, 0.0001),
                                  upper = c(100000, 100000, 100000))
  
  
  results[[ paste0("var_", i) ]] <- list("vgram" = vgram, "fitted_model" = fitted_model)
}

# save grid
tau_seq <- round(unique(sapply(results, function(a) 1-(a$fitted_model$space$psill[1]))), 4)
phi_seq <- round(unique(sapply(results, function(a) 1/a$fitted_model$space$range[2])), 4)

save(tau_seq, phi_seq, file = "data/variogram_informed_grid.Rdata")


# spatiotemporal residual variogram  --------------------------------------

results <- list()
initial_model <- vgmST(
  stModel = "separable",
  space = vgm(psill = 1, model = "Exp", nugget = 0.1, range = 10),
  time  = vgm(psill = 1, model = "Exp", nugget = 0.1, range = 5),
  sill  = 1)
for (i in 1:q) {
  data_i <- spatial_data[, i, ]

  predictor_i <- spatial_predictors
  eps_i <- matrix(0, n, t)

  for (j in 1:t) {
    P_t <- scale(predictor_i[,,j])
    L_t <- matrix(0, n, 11)
    if (format(time_points, "%m")[j] != "12") {
      L_t[,as.numeric(format(time_points, "%m"))[j]] <- 1
    }
    eps_i[,j] <- residuals(lm(data_i[,j] ~ cbind(P_t, L_t)))
  }
  
  data_vec <- as.vector(eps_i)
  
  stfdf <- STFDF(sp = SpatialPoints(coords),
                 time = as.POSIXct(1:dim(data_i)[2], origin = "2002-01-01"),
                 data = data.frame(value = data_vec))
  
  vgram <- variogramST(value ~ 1, data = stfdf, tlags = seq(0, 120, 10), cores = 12)
  
  fitted_model <- fit.StVariogram(object = vgram, model = initial_model, method = "L-BFGS-B",
                                  lower = c(0.0001, 0.0001, 0.0001),
                                  upper = c(100000, 100000, 100000))
  
  
  results[[ paste0("var_", i) ]] <- list("vgram" = vgram, "fitted_model" = fitted_model)
}

plot(results$var_1$vgram, results$var_1$fitted_model)
print(results$var_1$fitted_model)

plot(results$var_2$vgram, results$var_2$fitted_model)
print(results$var_2$fitted_model)

plot(results$var_3$vgram, results$var_3$fitted_model)
print(results$var_3$fitted_model)

plot(results$var_4$vgram, results$var_4$fitted_model)
print(results$var_3$fitted_model)

combined_plot <- wrap_plots(as.ggplot(plot(results$var_1$vgram, results$var_1$fitted_model, main = "Monthly temperature")),
                            as.ggplot(plot(results$var_2$vgram, results$var_2$fitted_model, main = "Monthly precipitation")),
                            as.ggplot(plot(results$var_3$vgram, results$var_3$fitted_model, main = "Monthly wind speed")),
                            as.ggplot(plot(results$var_4$vgram, results$var_4$fitted_model, main = "Monthly evaporation")),
                            ncol = 2)
ggsave("plots/copernicus_eda_stvariogram.png", combined_plot, width = 10.965, height = 6.63, dpi = 320, type = "cairo")


# hovmoller ---------------------------------------------------------------

lon_u <- sort(unique(lon))
lat_u <- sort(unique(lat))

all_vars_df <- map_dfr(1:q, function(v) {
  df <- as.data.frame(spatial_data[, v, ])
  colnames(df) <- as.character(time_points)
  df <- df %>%
    mutate(Longitude = coords[,1], Latitude = coords[,2], Variable = v) %>%
    pivot_longer(cols = matches("^\\d{4}-\\d{2}-\\d{2}$"),
                 names_to = "Time",
                 values_to = "Value") %>%
    mutate(Time = as.Date(Time))
  df
})

all_vars_df <- all_vars_df %>%
  mutate(Variable = recode_factor(as.character(Variable),
                                  `1` = "Monthly temperature",
                                  `2` = "Monthly precipitation",
                                  `3` = "Monthly wind speed",
                                  `4` = "Monthly evaporation"),
         Value = ifelse(is.nan(Value), NA, Value))

var_levels <- levels(all_vars_df$Variable)
complete_lat_grid <- expand.grid(Variable = var_levels, Latitude = lat_u, Time = time_points, stringsAsFactors = FALSE)
complete_lon_grid <- expand.grid(Variable = var_levels, Longitude = lon_u, Time = time_points, stringsAsFactors = FALSE)

hov_lat_raw <- all_vars_df %>%
  group_by(Variable, Latitude, Time) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  right_join(complete_lat_grid, by = c("Variable","Latitude","Time")) %>%
  arrange(Variable, Latitude, Time)

hov_lon_raw <- all_vars_df %>%
  group_by(Variable, Longitude, Time) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  right_join(complete_lon_grid, by = c("Variable","Longitude","Time")) %>%
  arrange(Variable, Longitude, Time)

time_num <- as.numeric(time_points)
dt_median <- median(diff(time_num), na.rm = TRUE)

overlap_factor <- 1.02
tile_width_days <- dt_median * overlap_factor

dlat <- median(diff(lat_u), na.rm = TRUE)
dlon <- median(diff(lon_u), na.rm = TRUE)
tile_height_lat <- dlat * overlap_factor
tile_height_lon <- dlon * overlap_factor

interp_fill <- function(x, tnum) {
  if (all(is.na(x))) return(x)
  approx_idx <- which(!is.na(x))
  if (length(approx_idx) < 2) {
    return(x)
  }
  interp_vals <- approx(x = tnum[approx_idx], y = x[approx_idx],
                        xout = tnum, rule = 2, ties = mean)$y
  return(interp_vals)
}

hov_lat_interp <- hov_lat_raw %>%
  group_by(Variable, Latitude) %>%
  arrange(Time) %>%
  mutate(Time_num = as.numeric(Time),
         Value_interp = { 
           tmp <- Value
           tmp <- interp_fill(tmp, Time_num)
           tmp }) %>%
  ungroup()

hov_lon_interp <- hov_lon_raw %>%
  group_by(Variable, Longitude) %>%
  arrange(Time) %>%
  mutate(Time_num = as.numeric(Time),
         Value_interp = { tmp <- Value; tmp <- interp_fill(tmp, Time_num); tmp }) %>%
  ungroup()

smooth_width_half <- 3
hov_lat <- hov_lat_interp %>%
  group_by(Variable, Latitude) %>%
  arrange(Time) %>%
  mutate(Smoothed = slider::slide_dbl(Value_interp, ~ mean(.x, na.rm = TRUE),
                                      .before = smooth_width_half, .after = smooth_width_half,
                                      .complete = FALSE)) %>%
  ungroup()

hov_lon <- hov_lon_interp %>%
  group_by(Variable, Longitude) %>%
  arrange(Time) %>%
  mutate(Smoothed = slider::slide_dbl(Value_interp, ~ mean(.x, na.rm = TRUE),
                                      .before = smooth_width_half, .after = smooth_width_half,
                                      .complete = FALSE)) %>%
  ungroup()

make_hov_plot_fixed <- function(varname, space = c("Latitude", "Longitude")) {
  space <- match.arg(space)
  if (space == "Latitude") {
    dfplot <- hov_lat %>% filter(Variable == varname)
    yvar <- "Latitude"
    tile_h <- tile_height_lat
  } else {
    dfplot <- hov_lon %>% filter(Variable == varname)
    yvar <- "Longitude"
    tile_h <- tile_height_lon
  }
  
  rng <- range(dfplot$Smoothed, na.rm = TRUE)
  if (!all(is.finite(rng))) rng <- c(0, 1)
  
  ggplot(dfplot, aes(x = Time, y = !!sym(yvar), fill = Smoothed)) +
    geom_tile(width = tile_width_days, height = tile_h, color = NA, show.legend = TRUE) +
    scale_fill_viridis_c(option = "plasma", limits = rng, na.value = "grey90") +
    scale_x_date(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(title = varname, x = "Time", y = yvar, fill = "") +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 14),
          legend.key.height = unit(0.9, "cm"))
}

plots_lat <- map(var_levels, ~ make_hov_plot_fixed(.x, "Latitude"))
plots_lon <- map(var_levels, ~ make_hov_plot_fixed(.x, "Longitude"))

final_figure_lat <- wrap_plots(plots_lat, ncol = 2)
final_figure_lon <- wrap_plots(plots_lon, ncol = 2)

print(final_figure_lat)
print(final_figure_lon)

ggsave("plots/copernicus_eda_hovmoller_lat.png", final_figure_lat, width = 10.965, height = 6.63, dpi = 320, type = "cairo")
ggsave("plots/copernicus_eda_hovmoller_lon.png", final_figure_lon, width = 10.965, height = 6.63, dpi = 320, type = "cairo")


# correlation analysis ----------------------------------------------------

nvars <- q
varnames <- c("T","R","W","E")

var_pairs <- combn(varnames, 2, simplify = FALSE)

compute_corr_one_time <- function(t){
  vals <- spatial_data[, , t]
  colnames(vals) <- varnames
  
  map_dfr(var_pairs, \(pair){
    c12 <- cor(vals[, pair[1]], vals[, pair[2]], use = "complete.obs")
    tibble(
      Pair = paste0("rho[list(",pair[1],",",pair[2],")]"),
      Var1 = pair[1],
      Var2 = pair[2],
      Time = time_points[t],
      Correlation = c12
    )
  })
}

corr_df <- map_dfr(seq_along(time_points), compute_corr_one_time)
corr_df <- corr_df %>%
  group_by(Pair) %>%
  arrange(Time) %>%
  mutate(Corr_smooth = slide_dbl(Correlation,
                                 mean, na.rm = TRUE,
                                 .before = 3, .after = 3)) %>%
  ungroup()

pair_levels <- unique(corr_df$Pair)
corr_df <- corr_df %>%
  mutate(Pair = factor(Pair, levels = pair_levels))

pair_labels <- setNames(
  sapply(pair_levels, function(p) p),  # keep the original text
  pair_levels
)

theme_pub <- theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0),
    axis.title = element_text(size = 15),
    axis.text  = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = "bottom",
    legend.key.height = unit(0.4, "cm"),
    legend.key.width  = unit(1.0, "cm")
  )

p_hov <- ggplot(corr_df, aes(x = Time, y = Pair, fill = Corr_smooth)) +
  geom_tile(width = diff(range(time_points)) / length(time_points),
            height = 0.9) +
  scale_fill_viridis_c(option = "plasma", limits = c(-1,1),
                       name = bquote(rho)) +
  scale_y_discrete(
    labels = function(x) parse(text = x)  # parse the strings into expressions
  ) +
  labs(title = "", x = "Time", y = "") +
  theme_pub +
  theme(legend.position = "right")

p_ts_facet <- ggplot(corr_df, aes(Time, Corr_smooth, color = Pair)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  scale_color_viridis_d(
    option = "plasma",
    name = NULL,
    labels = function(x) parse(text = x)  # legend uses expressions
  ) +
  labs(
    title = "",
    y = "",
    x = "Time"
  ) +
  facet_wrap(~ Pair, labeller = label_parsed) +  # parse facet labels as math
  theme_pub +
  theme(legend.position = "right")

# Save high-res PNG with Cairo
ggsave("plots/copernicus_eda_corr.png", wrap_plots(p_hov, p_ts_facet, nrow = 2), width = 10.965, height = 6.63, dpi = 320, type = "cairo")


# descriptive analysis ----------------------------------------------------

resp <- spatial_data
pred <- spatial_predictors
times <- time_points

resp_names <- c("Monthly temperature","Monthly precipitation","Monthly wind speed","Monthly evaporation")
pred_names <- c("Monthly cloud cover percentage","Monthly solar radiation","Monthly sea level pressure","Monthly humidity")
rownames(coords) <- seq_len(n)

resp_overall <- data.table(variable = resp_names,
                           mean = numeric(q),
                           median = numeric(q),
                           sd = numeric(q),
                           min = numeric(q),
                           max = numeric(q),
                           IQR = numeric(q),
                           n_obs = n * t)

for(j in seq_len(q)) {
  v <- as.vector(resp[, j, ])
  resp_overall[j, `:=` (mean = mean(v, na.rm=TRUE),
                        median = median(v, na.rm=TRUE),
                        sd = sd(v, na.rm=TRUE),
                        min = min(v, na.rm=TRUE),
                        max = max(v, na.rm=TRUE),
                        IQR = IQR(v, na.rm=TRUE))]
}

pred_overall <- data.table(variable = pred_names,
                           mean = numeric(p),
                           median = numeric(p),
                           sd = numeric(p),
                           min = numeric(p),
                           max = numeric(p),
                           IQR = numeric(p),
                           n_obs = n * t)
for(j in seq_len(p)) {
  v <- as.vector(pred[, j, ])
  pred_overall[j, `:=` (mean = mean(v, na.rm=TRUE),
                        median = median(v, na.rm=TRUE),
                        sd = sd(v, na.rm=TRUE),
                        min = min(v, na.rm=TRUE),
                        max = max(v, na.rm=TRUE),
                        IQR = IQR(v, na.rm=TRUE))]
}

resp_loc_mean <- matrix(NA, nrow = n, ncol = q)
resp_loc_sd   <- matrix(NA, nrow = n, ncol = q)
for(j in seq_len(q)){
  resp_loc_mean[, j] <- rowMeans(resp[, j, ], na.rm = TRUE)
  resp_loc_sd[, j]   <- apply(resp[, j, ], 1, sd, na.rm=TRUE)
}
colnames(resp_loc_mean) <- resp_names
colnames(resp_loc_sd)   <- paste0(resp_names, "_sd")

resp_spatial_summary <- data.table(variable = resp_names,
                                   loc_median = apply(resp_loc_mean, 2, median, na.rm=TRUE),
                                   loc_Q1 = apply(resp_loc_mean, 2, quantile, probs=0.25, na.rm=TRUE),
                                   loc_Q3 = apply(resp_loc_mean, 2, quantile, probs=0.75, na.rm=TRUE),
                                   loc_min = apply(resp_loc_mean, 2, min, na.rm=TRUE),
                                   loc_max = apply(resp_loc_mean, 2, max, na.rm=TRUE),
                                   loc_sd_median = apply(resp_loc_sd, 2, median, na.rm=TRUE))

resp_time_mean <- matrix(NA, nrow = t, ncol = q)
resp_time_sd   <- matrix(NA, nrow = t, ncol = q)
for(j in seq_len(q)){
  resp_time_mean[, j] <- colMeans(resp[, j, ], na.rm = TRUE)
  resp_time_sd[, j]   <- apply(resp[, j, ], 2, sd, na.rm = TRUE)
}
colnames(resp_time_mean) <- resp_names
colnames(resp_time_sd)   <- paste0(resp_names, "_sd")

resp_temporal_summary <- data.table(variable = resp_names,
                                    mean_of_time_means = apply(resp_time_mean, 2, mean, na.rm=TRUE),
                                    sd_of_time_means = apply(resp_time_mean, 2, sd, na.rm=TRUE),
                                    time_mean_min = apply(resp_time_mean, 2, min, na.rm=TRUE),
                                    time_mean_max = apply(resp_time_mean, 2, max, na.rm=TRUE))

theme_clean_img <- theme_void() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(0,0,0,0)
  )
for (j in seq_len(q)) {
  
  vname <- resp_names[j]
  v <- as.vector(resp[, j, ])
  df <- data.frame(value = v)
  
  p_hist <- ggplot(df, aes(value)) +
    geom_histogram(bins = 40, fill = "black", color = "black") +
    theme_clean_img
  
  ggsave(
    filename = sprintf("plots/hist_resp_%02d.png", j),
    plot = p_hist, dpi = 320,
    bg = "transparent"
  )
  
  p_box <- as.ggplot(~boxplot(v, col = NA, border = "black",
                              horizontal = T, outline = F, axes = F, lwd = 4))
  
  ggsave(
    filename = sprintf("plots/box_resp_%02d.png", j),
    plot = p_box, dpi = 320,
    bg = "transparent"
  )
}

latex_table <- data.table(
  Variable = resp_overall$variable,
  Mean = sprintf("%.3f", resp_overall$mean),
  SD   = sprintf("%.3f", resp_overall$sd),
  Min  = sprintf("%.3f", resp_overall$min),
  Max  = sprintf("%.3f", resp_overall$max),
  Hist = sprintf("plots/hist_resp_%02d.png", seq_len(q)),
  Box  = sprintf("plots/box_resp_%02d.png", seq_len(q))
)


