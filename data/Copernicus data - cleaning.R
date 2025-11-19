

# data cleaning -----------------------------------------------------------

library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting

nc_data <- nc_open("D:\\presi\\Downloads\\copernicus data\\evspsbl_CMIP6_historical_mon_18500101-20141201.nc")
print(nc_data)
lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat")
t <- ncvar_get(nc_data, "time")
evs <- ncvar_get(nc_data, "evspsbl")

evs |> dim()

evs[1:2, 1:2, 1, 1] # memebr = different simulation from climate model

nc_data2 <- nc_open("D:\\presi\\Downloads\\copernicus data\\sfcwind_CMIP6_historical_mon_18500101-20141201.nc")
print(nc_data2)
lon2 <- ncvar_get(nc_data2, "lon")
lat2 <- ncvar_get(nc_data2, "lat")
t2 <- ncvar_get(nc_data2, "time")
sfcwind <- ncvar_get(nc_data2, "sfcwind")

sfcwind[1:2, 1:2, 1, 1] # memebr = different simulation from climate model

lon == lon2
lat == lat2
t == t2


expand.grid(lon, lat)


nc_data3 <- nc_open("D:\\presi\\Downloads\\copernicus data\\t_CMIP6_historical_mon_18500101-20141201.nc")
print(nc_data3)
lon3 <- ncvar_get(nc_data3, "lon")
lat3 <- ncvar_get(nc_data3, "lat")
t3 <- ncvar_get(nc_data3, "time")
tempr <- ncvar_get(nc_data3, "t")

lon == lon3
lat == lat3
t == t3


nc_data4 <- nc_open("D:\\presi\\Downloads\\copernicus data\\r_CMIP6_historical_mon_18500101-20141201.nc")
print(nc_data4)
lon4 <- ncvar_get(nc_data4, "lon")
lat4 <- ncvar_get(nc_data4, "lat")
t4 <- ncvar_get(nc_data4, "time")
rains <- ncvar_get(nc_data4, "r")


# average over members - simulations
tempr2 <- apply(tempr, c(1,2,3), mean)
dimnames(tempr2) <- list(lon, lat)
tempr2 |> dim()

rains2 <- apply(rains, c(1,2,3), mean)
dimnames(rains2) <- list(lon, lat)
rains2 |> dim()

wind2 <- apply(sfcwind, c(1,2,3), mean)
evs2 <- apply(wind2, c(1,2,3), mean)
wind2 |> dim()

evs2 <- apply(evs, c(1,2,3), mean)
dimnames(evs2) <- list(lon, lat)
evs2 |> dim()

rm(list = c("rains", "tempr", "evs", "sfcwind"))

coords <- expand.grid(lon, lat)
n <- nrow(coords)
time_points <- seq(from = as.Date("1850-01-01"), by = "month", length.out = length(t))
# coords_name <- character(n)
# for (i in 1:n) { coords_name[i] <- paste0(as.character(coords[i,1]), " , ", as.character(coords[i,2])) }
var_name <- c("Temp","Rain","Wind","Evps")

spatial_data <- array(0, c(n, 4, length(t)))
dimnames(spatial_data) <- list(NULL, var_name, time_points)
for (tt in seq(length(t))) {
  
  spatial_data[,,tt] <- cbind(#coords,
                             c(tempr2[,,tt]),
                             c(rains2[,,tt]),
                             c(wind2[,,tt]),
                             c(evs2[,,tt]))
  
}

coords <- as.matrix(coords)
save(spatial_data, coords, file = "data/copernicus_data.Rdata")

# plot data ---------------------------------------------------------------


df <- data.frame(
  Longitude = lon,                 # numeric vector
  Latitude = lat,                  # numeric vector
  Value = spatial_data[,1,1]       # variable of interest
)

library(ggplot2)

ggplot(df, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_tile() +
  coord_fixed(xlim = c(-25, 45), ylim = c(34, 72), expand = FALSE) +  # adjust limits for Europe
  scale_fill_viridis_c() +
  labs(title = "Variable Surface over Europe", fill = "Value") +
  theme_minimal()


library(rnaturalearth)
library(sf)

# Get Europe boundaries
europe <- ne_countries(continent = "Europe", scale = "medium", returnclass = "sf")

# Plot with overlay
ggplot() +
  geom_tile(data = df, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_sf(data = europe, fill = NA, color = "black", size = 0.3) +
  coord_sf(xlim = c(-10, 40), ylim = c(30, 60), expand = FALSE) +
  scale_fill_viridis_c() +
  labs(title = "Variable Surface over Europe", fill = "Value") +
  theme_minimal()


# interpolation

df <- data.frame(
  Longitude = lon,                 # numeric vector
  Latitude = lat,                  # numeric vector
  Value = spatial_df[,1,1] #spatial_data[,1,],       # variable of interest
)
# Install if needed
library(akima)

# Interpolate to grid
interp_res <- interp(
  x = df$Longitude,
  y = df$Latitude,
  z = df$Value,
  duplicate = "mean",     # Handles duplicates
  extrap = FALSE          # Set to TRUE if you want extrapolation
)

# Convert result to data frame for ggplot
interp_df <- expand.grid(
  x = interp_res$x,
  y = interp_res$y
)
interp_df$z <- as.vector(interp_res$z)
names(interp_df) <- c("Longitude", "Latitude", "Value")

ggplot(interp_df, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_raster(interpolate = TRUE) +  # smoother rendering
  scale_fill_viridis_c() +
  coord_fixed(xlim = c(-10, 40), ylim = c(30, 60), expand = FALSE) +
  labs(title = "Interpolated Surface over Europe", fill = "Value") +
  theme_minimal()

library(rnaturalearth)
library(sf)

europe <- ne_countries(continent = "Europe", scale = "medium", returnclass = "sf")

ggplot() +
  geom_raster(data = interp_df, aes(x = Longitude, y = Latitude, fill = Value), interpolate = TRUE) +
  geom_sf(data = europe, fill = NA, color = "black", size = 0.3) +
  scale_fill_viridis_c() +
  coord_sf(xlim = c(-10, 40), ylim = c(30, 60), expand = FALSE) +
  labs(title = "Interpolated Variable Surface over Europe", fill = "Value") +
  theme_minimal()


plot_europe_surface <- function(lon, lat, data, var) {
  
  library(akima)
  library(rnaturalearth)
  library(sf)
  library(ggplot2)
  library(rlang)  # For .data tidy eval
  
  # Get the actual variable vector and column name
  if (is.numeric(var)) {
    if(is.null(colnames(data))) {var_name <- as.character(var)} else {var_name <- colnames(data)[var]}
    z <- data[,var]
  } else {
    var_name <- var
    z <- data[,var]
  }
  
  # Remove NAs
  valid <- !is.na(lon) & !is.na(lat) & !is.na(z)
  lon <- lon[valid]
  lat <- lat[valid]
  z <- z[valid]
  
  # Interpolation
  interp_res <- interp(x = lon, y = lat, z = z, duplicate = "mean", extrap = FALSE)
  
  # Prepare data frame
  interp_df <- expand.grid(
    x = interp_res$x,
    y = interp_res$y
  )
  interp_df[[var_name]] <- as.vector(interp_res$z)
  names(interp_df)[1:2] <- c("Longitude", "Latitude")
  
  # Get Europe shapefile
  europe <- ne_countries(continent = "Europe", scale = "medium", returnclass = "sf")
  
  # Plot
  ggplot() +
    geom_raster(data = interp_df, aes(x = Longitude, y = Latitude, fill = .data[[var_name]]), interpolate = TRUE) +
    geom_sf(data = europe, fill = NA, color = "black", size = 0.3) +
    scale_fill_viridis_c() +
    coord_sf(xlim = range(lon), ylim = range(lat), expand = FALSE) +
    labs(title = paste0("Interpolated ", var_name, " Surface over Europe"), fill = var_name) +
    theme_bw()
}


plot_europe_surface(lon = lon, lat = lat, data = spatial_df[,,1], var = 1)





library(akima)
library(rnaturalearth)
library(sf)
library(ggplot2)
library(rlang)
library(patchwork)

# Core function: plots a single surface
plot_europe_surface <- function(coords, data, var, time_label = NULL) {
  library(akima)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(rlang)
  
  # Handle variable selection
  if (is.numeric(var)) {
    var_name <- names(data)[var]
    z <- data[[var]]
  } else {
    var_name <- var
    z <- data[[var]]
  }
  
  # Extract lon/lat
  lon <- coords[, 1]
  lat <- coords[, 2]
  
  # Filter out NAs
  valid <- !is.na(lon) & !is.na(lat) & !is.na(z)
  lon <- lon[valid]
  lat <- lat[valid]
  z <- z[valid]
  
  # Interpolate
  interp_res <- interp(x = lon, y = lat, z = z, duplicate = "mean", extrap = FALSE)
  
  # Build data frame
  interp_df <- expand.grid(x = interp_res$x, y = interp_res$y)
  interp_df[[var_name]] <- as.vector(interp_res$z)
  names(interp_df)[1:2] <- c("Longitude", "Latitude")
  
  # Get map
  europe <- ne_countries(continent = "Europe", scale = "medium", returnclass = "sf")
  
  # Plot
  ggplot() +
    geom_raster(data = interp_df, aes(x = Longitude, y = Latitude, fill = .data[[var_name]]), interpolate = TRUE) +
    geom_sf(data = europe, fill = NA, color = "black", size = 0.3) +
    scale_fill_viridis_c(name = var_name) +
    coord_sf(xlim = range(lon), ylim = range(lat), expand = FALSE) +
    # labs(title = paste("Time:", time_label)) +
    labs(title = time_label) +
    theme_minimal()
}

# Wrapper function: handles 3D array and creates 3x3 grid of plots
plot_europe_surface_timeseries <- function(coords, data_3d, var,
                                           start_date = NULL,
                                           time_unit = "days",  # NEW ARGUMENT
                                           max_plots = 9) {
  library(patchwork)
  library(lubridate)
  
  n_time <- dim(data_3d)[3]
  max_plots <- min(max_plots, n_time)
  
  # Generate evenly spaced time indices
  raw_indices <- unique(round(seq(1, n_time, length.out = max_plots)))
  time_indices <- sort(unique(c(raw_indices, n_time)))[1:max_plots]
  
  plots <- list()
  
  for (i in seq_along(time_indices)) {
    t <- time_indices[i]
    
    if (t > n_time || is.na(t)) next
    
    data_slice <- as.data.frame(data_3d[,,t])
    
    # Handle date label generation
    time_label <- if (!is.null(start_date)) {
      start <- as.Date(start_date)
      label_date <- switch(
        time_unit,
        "days" = start + days(t - 1),
        "months" = start %m+% months(t - 1),
        "years" = start %m+% years(t - 1),
        stop("Invalid time_unit: must be 'days', 'months', or 'years'")
      )
      format(label_date, "%Y-%m-%d")
    } else {
      paste0("t=", t)
    }
    
    plots[[i]] <- plot_europe_surface(coords, data_slice, var, time_label)
  }
  
  # Fill empty slots if needed
  while (length(plots) < max_plots) {
    plots[[length(plots) + 1]] <- ggplot() + theme_void()
  }
  
  wrap_plots(plots, ncol = ceiling(sqrt(max_plots)))
}

plot_europe_surface_variables <- function(coords, data_3d, 
                                          time_point = 1,
                                          vars = NULL,
                                          var_names = NULL,
                                          max_plots = 9) {
  library(patchwork)
  
  # Detect if 2D or 3D
  is_3d <- length(dim(data_3d)) == 3
  
  # Extract variable slice
  if (is_3d) {
    n_vars <- dim(data_3d)[2]
    data_slice <- as.data.frame(data_3d[,,time_point])
  } else {
    n_vars <- ncol(data_3d)
    data_slice <- as.data.frame(data_3d)
  }
  
  # If no variable subset specified, use all
  if (is.null(vars)) {
    vars <- seq_len(n_vars)
  }
  
  vars <- vars[vars <= n_vars]
  max_plots <- min(length(vars), max_plots)
  vars <- vars[1:max_plots]
  
  plots <- list()
  
  for (i in seq_along(vars)) {
    var_idx <- vars[i]
    
    # Label
    var_label <- if (!is.null(var_names)) {
      var_names[var_idx]
    } else {
      paste0("Var ", var_idx)
    }
    
    plots[[i]] <- plot_europe_surface(coords, data_slice, var_idx, var_label)
  }
  
  # Fill empty slots if needed
  while (length(plots) < max_plots) {
    plots[[length(plots) + 1]] <- ggplot() + theme_void()
  }
  
  wrap_plots(plots, ncol = ceiling(sqrt(max_plots)))
}




get_palette_by_index <- function(index) {
  if (index == 1) {
    return(list(type = "distiller", palette = "RdBu", direction = -1))
  } else if (index == 4) {
    return(list(type = "viridis", option = "D", direction = -1))
  } else if (index == 3) {
    return(list(type = "distiller", palette = "Greens", direction = -1))
  } else if (index == 2) {
    return(list(type = "manual", colors = c("purple", "white")))
  } else {
    return(list(type = "viridis", option = "D", direction = -1))
  }
}

plot_europe_surface <- function(coords, data, var, time_label = NULL) {
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
  
  ggplot() +
    geom_raster(data = interp_df, aes(x = Longitude, y = Latitude, fill = .data[[var_name]]), interpolate = TRUE) +
    geom_sf(data = europe, fill = NA, color = "black", size = 0.3) +
    fill_scale +
    coord_sf(xlim = range(lon), ylim = range(lat), expand = FALSE) +
    labs(title = time_label) +
    theme_minimal()
}




plot_europe_surface_timeseries(coords, data_3d = spatial_df, var = "Temp", max_plots = 16)

plot_europe_surface_timeseries(
  coords = coords,
  data_3d = spatial_df,
  var = "Temp",
  start_date = "2000-01-01",
  time_unit = "months",   # <--- Try "days", "months", or "years"
  max_plots = 16
)

plot_europe_surface_variables(coords = coords,
                              data_3d = spatial_df,
                              time_point = 5,
                              max_plots = 4)


plot_europe_surface_variables(coords = coords,
                              data_3d = spatial_df,
                              time_point = 5,
                              vars = c(1, 3, 4),
                              var_names = c("Temp", "Wind", "Evps"))





# Hovmöller Diagram Along Latitude/Longitude

library(tidyr)
library(dplyr)

# Create a long-format data frame for ALL variables
all_vars_df <- NULL

for (v in 1:4) {
  temp <- as.data.frame(spatial_data[ , v, ])
  colnames(temp) <- as.character(time_points)
  
  temp$Latitude <- locations$Latitude
  temp$Longitude <- locations$Longitude
  temp$Variable <- paste0("Var_", v)
  
  temp_long <- temp %>%
    pivot_longer(
      cols = -c(Latitude, Longitude, Variable),
      names_to = "Time",
      values_to = "Value"
    )
  
  all_vars_df <- bind_rows(all_vars_df, temp_long)
}

# Convert Time to Date (if needed)
all_vars_df$Time <- as.Date(all_vars_df$Time)

# Convert time
all_vars_df$Time <- as.Date(all_vars_df$Time)

# Clean any NaN from mean calculations
all_vars_df$Value[is.nan(all_vars_df$Value)] <- NA

# Recode variable names
all_vars_df$Variable <- recode(all_vars_df$Variable,
                               Var_1 = "Temperature",
                               Var_2 = "Precipitation",
                               Var_3 = "Wind",
                               Var_4 = "Evaporation")

# Group and compute smoothed mean over latitude (for Hovmöller Latitude–Time)
hov_lat <- all_vars_df %>%
  group_by(Variable, Latitude, Time) %>%
  summarize(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  group_by(Variable, Latitude) %>%
  arrange(Time) %>%
  mutate(Smoothed = zoo::rollapply(Value, width = 5, FUN = mean, fill = NA, align = "center", na.rm = TRUE)) %>%
  ungroup()

hov_lat <- hov_lat %>% filter(as.integer(Time) > "1950-12-01")
# Plot Hovmöller (Latitude vs Time) faceted by variable
ggplot(hov_lat, aes(x = Time, y = Latitude, fill = Smoothed)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", na.value = "white") +
  facet_wrap(~ Variable, scales = "free", ncol = 2) +
  labs(title = "Hovmöller Diagram (Latitude vs Time)", 
       x = "Time", y = "Latitude", fill = "Smoothed\nValue") +
  theme_minimal(base_size = 14)

hov_lat <- all_vars_df %>%
  group_by(Variable, Latitude, Time) %>%
  summarize(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  mutate(Value = ifelse(is.nan(Value), NA, Value))
# hov_lat <- hov_lat %>% filter(as.integer(Time) %% 10 == 0)
ggplot(hov_lat, aes(x = Time, y = Latitude, fill = Value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  facet_wrap(~ Variable, scales = "free", ncol = 2) +  # allows each to scale independently
  labs(title = "Hovmöller Diagram by Latitude", x = "Time", y = "Latitude", fill = "Value") +
  theme_minimal()



hov_lon <- all_vars_df %>%
  group_by(Variable, Longitude, Time) %>%
  summarize(Value = mean(Value, na.rm = TRUE), .groups = "drop")
# hov_lon <- hov_lon %>% filter(as.integer(Time) %% 10 == 0)
ggplot(hov_lon, aes(x = Time, y = Longitude, fill = Value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  facet_wrap(~ Variable, ncol = 2) +
  labs(title = "Hovmöller Diagrams by Longitude (Faceted by Variable)", fill = "Value") +
  theme_minimal()


