# re-install damaged packages
# detach("package:gstat", unload = TRUE)
# install.packages("gstat")
library(gstat)
# install.packages("spacetime")
library(spacetime)

# Working packages
library(mniw)
library(classInt)
library(RColorBrewer)
library(sp)
library(fields)
library(parallel)
library(foreach)
library(doParallel)
library(geoR)
library(corrplot)
library(sp)
library(mapproj)
library(patchwork)
library(akima)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rlang)
library(dplyr)
library(tidyr)
library(purrr)

# Associated R package
# devtools::install_github("lucapresicce/spFFBS") # if not previously installed
library(spFFBS)
library(spBPS)


# Working R functions -----------------------------------------------------

get_palette_by_index <- function(index) {
  if (index == 1) {
    return(list(type = "distiller", palette = "RdBu", direction = -1))
  } else if (index == 4) {
    return(list(type = "viridis", option = "D", direction = -1))
  } else if (index == 3) {
    return(list(type = "viridis", option = "E", direction = -1))
  } else if (index == 2) {
    return(list(type = "viridis", option = "C", direction = -1))
  } else {
    return(list(type = "viridis", option = "D", direction = -1))
  }
}

plot_europe_surface <- function(coords, data, var, time_label = NULL, show_title = TRUE) {
  library(akima)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(rlang)
  
  if (is.numeric(var)) {
    var_name <- names(data)[var]
    z <- data[[var]]
    var_index <- var
  } else {
    var_name <- var
    z <- data[[var]]
    var_index <- which(names(data) == var)
  }
  
  lon <- coords[, 1]
  lat <- coords[, 2]
  
  valid <- !is.na(lon) & !is.na(lat) & !is.na(z)
  lon <- lon[valid]
  lat <- lat[valid]
  z <- z[valid]
  
  interp_res <- interp(x = lon, y = lat, z = z, duplicate = "mean", extrap = FALSE)
  interp_df <- expand.grid(x = interp_res$x, y = interp_res$y)
  interp_df[[var_name]] <- as.vector(interp_res$z)
  names(interp_df)[1:2] <- c("Longitude", "Latitude")
  
  europe <- ne_countries(continent = "Europe", scale = "medium", returnclass = "sf")
  palette <- get_palette_by_index(var_index)
  
  fill_scale <- switch(
    palette$type,
    "distiller" = scale_fill_distiller(
      palette = palette$palette,
      direction = palette$direction,
      name = var_name,
      na.value = "transparent"
    ),
    "viridis" = scale_fill_viridis_c(
      option = palette$option,
      direction = palette$direction,
      name = var_name,
      na.value = "transparent"
    ),
    "manual" = scale_fill_gradientn(
      colors = palette$colors,
      name = var_name,
      na.value = "transparent"
    )
  )
  
  gg <- ggplot() +
    geom_raster(data = interp_df, aes(x = Longitude, y = Latitude, fill = .data[[var_name]]), interpolate = TRUE) +
    geom_sf(data = europe, fill = NA, color = "black", size = 0.3) +
    fill_scale +
    coord_sf(xlim = range(lon), ylim = range(lat), expand = FALSE) +
    theme_minimal()
  
  if (show_title && !is.null(time_label)) {
    gg <- gg + labs(title = time_label)
  }
  
  return(gg)
}

plot_europe_surface_timeseries <- function(coords, data_3d, var,
                                           start_date = NULL,
                                           time_unit = "days",
                                           max_plots = 9,
                                           show_titles = TRUE) {
  library(patchwork)
  library(lubridate)
  
  n_time <- dim(data_3d)[3]
  max_plots <- min(max_plots, n_time)
  raw_indices <- unique(round(seq(1, n_time, length.out = max_plots)))
  time_indices <- sort(unique(c(raw_indices, n_time)))[1:max_plots]
  
  plots <- list()
  
  for (i in seq_along(time_indices)) {
    t <- time_indices[i]
    if (t > n_time || is.na(t)) next
    data_slice <- as.data.frame(data_3d[,,t])
    
    time_label <- if (!is.null(start_date)) {
      start <- as.Date(start_date)
      label_date <- switch(
        time_unit,
        "days" = start + days(t - 1),
        "months" = start %m+% months(t - 1),
        "years" = start %m+% years(t - 1),
        stop("Invalid time_unit")
      )
      format(label_date, "%Y-%m-%d")
    } else {
      paste0("t=", t)
    }
    
    plots[[i]] <- plot_europe_surface(coords, data_slice, var, time_label, show_title = show_titles)
  }
  
  while (length(plots) < max_plots) {
    plots[[length(plots) + 1]] <- ggplot() + theme_void()
  }
  
  wrap_plots(plots, ncol = ceiling(sqrt(max_plots)))
}

plot_europe_surface_variables <- function(coords, data_3d, 
                                          time_point = 1,
                                          vars = NULL,
                                          var_names = NULL,
                                          max_plots = 9,
                                          show_titles = TRUE) {
  library(patchwork)
  
  is_3d <- length(dim(data_3d)) == 3
  if (is_3d) {
    n_vars <- dim(data_3d)[2]
    data_slice <- as.data.frame(data_3d[,,time_point])
  } else {
    n_vars <- ncol(data_3d)
    data_slice <- as.data.frame(data_3d)
  }
  
  if (is.null(vars)) {
    vars <- seq_len(n_vars)
  }
  
  vars <- vars[vars <= n_vars]
  max_plots <- min(length(vars), max_plots)
  vars <- vars[1:max_plots]
  
  plots <- list()
  
  for (i in seq_along(vars)) {
    var_idx <- vars[i]
    var_label <- if (!is.null(var_names)) var_names[var_idx] else paste0("Var ", var_idx)
    title <- if (show_titles) var_label else NULL
    plots[[i]] <- plot_europe_surface(coords, data_slice, var_idx, time_label = title, show_title = show_titles)
  }
  
  while (length(plots) < max_plots) {
    plots[[length(plots) + 1]] <- ggplot() + theme_void()
  }
  
  wrap_plots(plots, ncol = ceiling(sqrt(max_plots)))
}


# Data Loading ------------------------------------------------------------

# Load
load("data/copernicus_data.RData")
N <- dim(spatial_data)[1]
q <- dim(spatial_data)[2]
t <- dim(spatial_data)[3]

dimnames(spatial_data)[[3]] <- 1:1980

# eda  --------------------------------------------------------------------

results <- list()
initial_model <- vgmST(
  stModel = "separable",
  space = vgm(psill = 1, model = "Exp", nugget = 0.1, range = 10),
  time  = vgm(psill = 1, model = "Exp", nugget = 0.1, range = 5),
  sill  = 1)
for (i in 1:q) {
  data_i <- spatial_data[, i, (t-144+1):(t-24)]
  data_vec <- as.vector(data_i)

  stfdf <- STFDF(sp = SpatialPoints(coords),
                 time = as.POSIXct(1:dim(data_i)[2], origin = "2002-01-01"),
                 data = data.frame(value = data_vec))

  vgram <- variogramST(value ~ 1, data = stfdf, tlags = seq(0, 120, 10), cores = 4)

  fitted_model <- fit.StVariogram(object = vgram, model = initial_model, method = "L-BFGS-B",
                                  lower = c(0.0001, 0.0001, 0.0001),
                                  upper = c(100000, 100000, 100000))

  results[[paste0("var_", i)]] <- fitted_model
}

print(results)


# Setting Bayesian Predictive Stacking ------------------------------------

# Define competetive (J) models
tau_seq <- unique(sapply(results, function(a) 1-(a$space$psill[1])))
phi_seq <- unique(sapply(results, function(a) 1/a$space$range[2]))
par_grid <- spBPS::expand_grid_cpp(rev(tau_seq), rev(phi_seq))
J <- nrow(par_grid)

# Data subsetting ---------------------------------------------------------

# number of selected times
tt <- 144

# Indexing
set.seed(4-8-15-16-23-42)
idx <- sample(1:nrow(spatial_data), 600, F)
spatial_df <- spatial_data[idx,,(t-tt+1):t]
head(spatial_df)

# New data dimensions
N <- dim(spatial_df)[1]
q <- dim(spatial_df)[2]
t <- dim(spatial_df)[3]

# coordinates
lon <- coords[idx, 1]
lat <- coords[idx, 2]
coord <- coords
coords <- coords[idx,]

# maximum inter-site distance
d.max <- sqrt((max(coord[,1]) - min(coord[,1]))^2 +
                (max(coord[,2]) - min(coord[,2]))^2)
# to multiply of about 111

# sinusoidally projected coordinates (scaled to 1000km units) as explanatory variables
coords.sinusoidal <- mapproject(coords[, 1], coords[, 2],
                                projection = "sinusoidal")

radius.of.earth <- 6.371            ## 6.371 * 1000 kilometers 
coords.sinusoidal <- radius.of.earth * (cbind(coords.sinusoidal$x,
                                              coords.sinusoidal$y))

# dates
dates_vector <- seq(from = as.Date("1850-01-01"), by = "1 month", length.out = 1980)[(1980-t+1):(1980)]

# seasonal dummy variables (departure from december)
design_matrix <- coords.sinusoidal
design_matrix <- array(rep(design_matrix, t), c(N, ncol(coords.sinusoidal), t))
design_matrix2 <- array(0, c(N, 2+11, t))
for (i in 1:t) {
  L_t <- matrix(0, N, 11)
  if (format(dates_vector, "%m")[i] != "12") {
    L_t[,as.numeric(format(dates_vector, "%m"))[i]] <- 1
  }
  design_matrix2[,,i] <- cbind(design_matrix[,,i], L_t)
  }
design_matrix <- design_matrix2
p <- ncol(design_matrix)

# subset for evaluating predictions
u <- 100 # holds-out locations
h <- 12  # holds-out timepoints

set.seed(4-8-15-16-23-42)
n <- N-u
train_set <- sample(1:N, n, F)
test_set  <- (1:N)[-train_set]

time_points  <- 1:(t-h)
forecast_set <- (t-h+1):t
tmax <- length(time_points)

# subsetting data
Ys <- spatial_df[train_set,,time_points]
Xs <- design_matrix[train_set,,time_points]
crds <- coords[train_set,]
Ps <- sapply(1:(tmax+h), function(s) cbind(design_matrix[train_set,,s], diag(n)), simplify = "array")
Gs <- array(rep(diag(n+p), tmax+h), c(n+p, n+p, tmax+h))
Ds <- spBPS::arma_dist(crds)

Yu <- spatial_df[test_set,,]
Xu <- design_matrix[test_set,,]
crdu <- coords[test_set,]
Du <- spBPS::arma_dist(crdu)

# prior information
m0     <- matrix(0, n+p, q)
C0 <- rbind( cbind(diag(0.005, p), matrix(0, p, n)), cbind(matrix(0, n, p), exp(-Ds)) )
nu0 <- 3
Psi0 <- diag(q)
prior <- list("m" = m0, "C" = C0, "nu" = nu0, "Psi" = Psi0)

# free memory
gc()

# DYNBPS forward filtering backward sampling ------------------------------

out <- spFFBS::spFFBS(Y = Ys, G = Gs, P = Ps, D = Ds,
                      grid = list(tau = tau_seq, phi = phi_seq), prior = prior, L = 200,
                      do_BS = T, do_forecast = T, do_spatial = T, tnew = h,
                      spatial = list(crd = crds, crdtilde = crdu, Xtilde = Xu,
                                     t = c(which(dates_vector == "2007-09-01"),
                                           which(dates_vector == "2013-09-01"))))


# Posterior evaluation ----------------------------------------------------

# RMSE
theta_postmean <- sapply(1:tmax, function(t)apply(out$BS[[t]], c(1,2), mean), simplify = "array")
theta_upp      <- sapply(1:tmax, function(t)apply(out$BS[[t]], c(1,2), quantile, 0.975), simplify = "array")
theta_low      <- sapply(1:tmax, function(t)apply(out$BS[[t]], c(1,2), quantile, 0.025), simplify = "array")


# Forecast evaluation -----------------------------------------------------

# RMSE
Y_postmean <- apply(out$forecast$Y_pred, c(1,2,4), mean)
Y_upp      <- apply(out$forecast$Y_pred, c(1,2,4), quantile, 0.975)
Y_low      <- apply(out$forecast$Y_pred, c(1,2,4), quantile, 0.025)
sq_res <- apply((Ys - Y_postmean[,,1:tmax])^2, c(1,2), mean)
sqrt( sq_res ) |> mean()

# RMSPE
sq_res <- apply((spatial_df[train_set,,forecast_set] - Y_postmean[,,-(1:tmax)])^2, c(1,2), mean)
sqrt( sq_res ) |> mean()

# naming
colnames(Y_postmean) <- colnames(Yu)

# Plot points selcetion time-forecast
lon <- crds[, 1]
lat <- crds[, 2]
europe <- ne_countries(continent = "Europe", scale = "medium", returnclass = "sf")
smp <- c(178, 262, 341, 422)
gg_points <- ggplot() +
  geom_sf(data = europe, fill = rgb(190/255, 220/255, 180/255, 0.5), color = "black", size = 0.1) +
  geom_point(aes(x = lon[smp], y = lat[smp]), col = "red", size = 18, pch = "*") +
  geom_text(aes(x = lon[smp], y = lat[smp], label = smp), nudge_x = -2, nudge_y = 1, col = "red", size = 5.5) +
  coord_sf(xlim = range(lon)*1.15, ylim = range(lat)*1.05, expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank())
make_df <- function(j) {
  map_dfr(1:q, function(v) {
    tibble(
      time = 1:(tmax+h),
      smp = j,
      var = v,
      spatial = spatial_df[train_set,,][j, v, ],
      low = Y_low[j, v, ],
      upp = Y_upp[j, v, ],
      postmean = Y_postmean[j, v, ]
    )
  })
}
plot_df <- map_dfr(smp, make_df)
gg_lines <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = spatial), color = "black") +
  geom_ribbon(
    aes(ymin = low, ymax = upp, fill = factor(var)),
    alpha = 0.5,
    colour = "red",
    linetype = "dashed",
    show.legend = TRUE
  ) +
  geom_line(aes(y = postmean), color = "red") +
  facet_grid(smp ~ var, scales = "free") +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "variables") +
  guides(fill = guide_legend(override.aes = list(colour = NA, linetype = 0))) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text.y = element_text(),
    axis.title.y = element_blank()
  )

library(ggh4x)
custom_labels <- c("Temp","Rain","Wind","Evps")
combined_pointslines <- gg_points + gg_lines +
  facet_grid2(rows = vars(smp), cols = vars(var),
    scales = "free", independent = "all", space = "free", axes = "all",
    labeller = labeller(var = setNames(custom_labels, 1:4)))
x11(); combined_pointslines

width <- 360*3
height <- 360*1.5
pointsize <- 16
png("plots/copernicus_temporal_forecast_points.png", width = width, height = height, pointsize = pointsize, family = "sans")
combined_pointslines
dev.off()

# Spatial surface results
# temp
width <- 360*6
height <- 360*3
pointsize <- 16
png("plots/copernicus_forecast_temp.png", width = width, height = height, pointsize = pointsize, family = "sans")
wrap_plots(
  plot_europe_surface_timeseries(
    coords = crds,
    data_3d = Ys,
    var = 1,
    start_date = "2002-12-01",
    time_unit = "months",
    max_plots = 9
  ),
  plot_spacer(),
  plot_europe_surface_timeseries(
    coords = crds,
    data_3d = Y_postmean[,,1:tmax],
    var = 1,
    start_date = "2002-12-01",
    time_unit = "months",
    max_plots = 9
  ),
  ncol = 3
) +
  plot_layout(widths = c(1, 0.2, 1))
dev.off()
# rain
width <- 360*6
height <- 360*3
pointsize <- 16
png("plots/copernicus_forecast_rain.png", width = width, height = height, pointsize = pointsize, family = "sans")
wrap_plots(
  plot_europe_surface_timeseries(
    coords = crds,
    data_3d = Ys,
    var = 2,
    start_date = "2002-12-01",
    time_unit = "months",
    max_plots = 9
  ),
  plot_spacer(),
  plot_europe_surface_timeseries(
    coords = crds,
    data_3d = Y_postmean[,,1:tmax],
    var = 2,
    start_date = "2002-12-01",
    time_unit = "months",
    max_plots = 9
  ),
  ncol = 3
) +
  plot_layout(widths = c(1, 0.2, 1))
dev.off()
# wind
width <- 360*6
height <- 360*3
pointsize <- 16
png("plots/copernicus_forecast_wind.png", width = width, height = height, pointsize = pointsize, family = "sans")
wrap_plots(
  plot_europe_surface_timeseries(
    coords = crds,
    data_3d = Ys,
    var = 3,
    start_date = "2002-12-01",
    time_unit = "months",
    max_plots = 9
  ),
  plot_spacer(),
  plot_europe_surface_timeseries(
    coords = crds,
    data_3d = Y_postmean[,,1:tmax],
    var = 3,
    start_date = "2002-12-01",
    time_unit = "months",
    max_plots = 9
  ),
  ncol = 3
) +
  plot_layout(widths = c(1, 0.2, 1))
dev.off()
# evps
width <- 360*6
height <- 360*3
pointsize <- 16
png("plots/copernicus_forecast_evps.png", width = width, height = height, pointsize = pointsize, family = "sans")
wrap_plots(
  plot_europe_surface_timeseries(
    coords = crds,
    data_3d = spatial_df[train_set,,forecast_set],
    var = 4,
    start_date = "2012-12-01",
    time_unit = "months",
    max_plots = 9
  ),
  plot_spacer(),
  plot_europe_surface_timeseries(
    coords = crds,
    data_3d = Y_postmean[,,forecast_set],
    var = 4,
    start_date = "2012-12-01",
    time_unit = "months",
    max_plots = 9
  ),
  ncol = 3
) +
  plot_layout(widths = c(1, 0.2, 1))
dev.off()


# Interpolation evaluation ------------------------------------------------

# IN-SAMPLE
# Extract predictions
# Ysp_pred   <- sapply(1:L, function(l){out_SI_i[[l]][1:u,]}   , simplify = "array")
L <- 200
Ysp_pred   <- sapply(1:L, function(l){out$spatial$t57[[l]][1:u,]}   , simplify = "array")

# RMSPE
Ysp_postmean <- apply(Ysp_pred, c(1,2), mean)
Ysp_upp      <- apply(Ysp_pred, c(1,2), quantile, 0.975)
Ysp_low      <- apply(Ysp_pred, c(1,2), quantile, 0.025)
sq_res <- (Yu[,,which(dates_vector == "2007-09-01")] - Ysp_postmean)^2
sqrt( sq_res ) |> mean()

# naming
colnames(Ysp_postmean) <- colnames(Yu)

# plotting
width <- 360*4
height <- 360*2
pointsize <- 16
png("plots/copernicus_interpolation_insample.png", width = width, height = height, pointsize = pointsize, family = "sans")
wrap_plots(
  plot_europe_surface_variables(coords = crdu, data = Yu, time_point = which(dates_vector == "2007-09-01"), show_titles = F),
  plot_spacer(),
  plot_europe_surface_variables(coords = crdu, data = Ysp_postmean, show_titles = F)) + 
  plot_layout(widths = c(1, 0.2, 1)) + 
  plot_annotation(title = paste0("Spatial interpolation at time: ", dates_vector[which(dates_vector == "2007-09-01")]),theme = theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, vjust = -10)))
dev.off()

# OUT-OF-SAMPLE
# Extract predictions
Ysp_pred   <- sapply(1:L, function(l){out$spatial$t129[[l]][1:u,]}   , simplify = "array")

# RMSPE
Ysp_postmean <- apply(Ysp_pred, c(1,2), mean)
Ysp_upp      <- apply(Ysp_pred, c(1,2), quantile, 0.975)
Ysp_low      <- apply(Ysp_pred, c(1,2), quantile, 0.025)
sq_res <- (Yu[,,which(dates_vector == "2013-09-01")] - Ysp_postmean)^2
sqrt( sq_res ) |> mean()

# result visualization
# plotting
width <- 360*4
height <- 360*2
pointsize <- 16
png("plots/copernicus_interpolation_outsample.png", width = width, height = height, pointsize = pointsize, family = "sans")
wrap_plots(
  plot_europe_surface_variables(coords = crdu, data = Yu, time_point = which(dates_vector == "2013-09-01"), show_titles = F),
  plot_spacer(),
  plot_europe_surface_variables(coords = crdu, data = Ysp_postmean, show_titles = F)) + 
  plot_layout(widths = c(1, 0.2, 1)) + 
  plot_annotation(title = paste0("Spatial interpolation at time: ", dates_vector[which(dates_vector == "2013-09-01")]),theme = theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, vjust = -10)))
dev.off()

