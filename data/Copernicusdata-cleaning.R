

# data cleaning -----------------------------------------------------------

library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting


# response variables ------------------------------------------------------

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


# gather data -------------------------------------------------------------

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

# predictor variables -----------------------------------------------------

# cloud cover percentage
nc_data <- nc_open("D:\\presi\\Downloads\\copernicus data\\clt_CMIP6_historical_mon_18500101-20141201.nc")
print(nc_data)
lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat")
t <- ncvar_get(nc_data, "time")
cloud <- ncvar_get(nc_data, "clt")

dim(cloud)
head(cloud/100)
# average over members - simulations
cloud2 <- apply(cloud, 1:3, mean)/100
dimnames(cloud2) <- list(lon, lat)
dim(cloud2)


# surface solar radiation downwards
nc_data2 <- nc_open("D:\\presi\\Downloads\\copernicus data\\rsds_CMIP6_historical_mon_18500101-20141201.nc")
print(nc_data2)
lon2 <- ncvar_get(nc_data2, "lon")
lat2 <- ncvar_get(nc_data2, "lat")
t2 <- ncvar_get(nc_data2, "time")
solar <- ncvar_get(nc_data2, "rsds")

dim(solar)
head(solar)
# average over members - simulations
solar2 <- apply(solar, 1:3, mean)
dimnames(solar2) <- list(lon, lat)
dim(solar2)


# sea level pressure
nc_data3 <- nc_open("D:\\presi\\Downloads\\copernicus data\\psl_CMIP6_historical_mon_18500101-20141201.nc")
print(nc_data3)
lon3 <- ncvar_get(nc_data3, "lon")
lat3 <- ncvar_get(nc_data3, "lat")
t3 <- ncvar_get(nc_data3, "time")
press <- ncvar_get(nc_data3, "psl")

dim(press)
head(press)
# average over members - simulations
press2 <- apply(press, 1:3, mean)
dimnames(press2) <- list(lon, lat)
dim(press2)


# near surface specific humidity
nc_data4 <- nc_open("D:\\presi\\Downloads\\copernicus data\\huss_CMIP6_historical_mon_18500101-20141201.nc")
print(nc_data4)
lon4 <- ncvar_get(nc_data4, "lon")
lat4 <- ncvar_get(nc_data4, "lat")
t4 <- ncvar_get(nc_data4, "time")
humid <- ncvar_get(nc_data4, "huss")

dim(humid)
head(humid)
# average over members - simulations
humid2 <- apply(humid, 1:3, mean)
dimnames(humid2) <- list(lon, lat)
dim(humid2)

# remove original data
rm(list = c("cloud", "solar", "press", "humid"))


# gather data -------------------------------------------------------------

coords <- expand.grid(lon, lat)
n <- nrow(coords)
time_points <- seq(from = as.Date("1850-01-01"), by = "month", length.out = length(t))
# coords_name <- character(n)
# for (i in 1:n) { coords_name[i] <- paste0(as.character(coords[i,1]), " , ", as.character(coords[i,2])) }
var_name <- c("Cloud","Radiation","Pressure","Humidity")

spatial_predictors <- array(0, c(n, 4, length(t)))
dimnames(spatial_predictors) <- list(NULL, var_name, time_points)
for (tt in seq(length(t))) {
  
  spatial_predictors[,,tt] <- cbind(#coords,
                                    c(cloud2[,,tt]),
                                    c(solar2[,,tt]),
                                    c(press2[,,tt]),
                                    c(humid2[,,tt]))
  
}

coords <- as.matrix(coords)
save(spatial_predictors, coords, file = "data/copernicus_predictors.Rdata")
