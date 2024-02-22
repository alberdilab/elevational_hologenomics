
#Extracting the information from the coordinate points

LandCover10 <- raster::raster("WORLDCOVER/ESA_WORLDCOVER_10M_2021_V200/MAP/ESA_WorldCover_10m_2021_v200_N42E000_Map/ESA_WorldCover_10m_2021_v200_N42E000_Map.tif")
LandCover12 <- raster::raster("WORLDCOVER/ESA_WORLDCOVER_10M_2021_V200/MAP/ESA_WorldCover_10m_2021_v200_N42E000_Map/ESA_WorldCover_10m_2021_v200_N42E000_InputQuality.tif")

extracted_values_new <- raster::extract(LandCover10, Coordinates)
extracted_values_new_c <- raster::extract(LandCover12, Coordinates)

LandCover11 <- raster::raster("WORLDCOVER/ESA_WORLDCOVER_10M_2021_V200/MAP/ESA_WorldCover_10m_2021_v200_N42W003_Map/ESA_WorldCover_10m_2021_v200_N42W003_Map.tif")
LandCover13 <- raster::raster("WORLDCOVER/ESA_WORLDCOVER_10M_2021_V200/MAP/ESA_WorldCover_10m_2021_v200_N42W003_Map/ESA_WorldCover_10m_2021_v200_N42W003_InputQuality.tif")

extracted_values_new_b <- raster::extract(LandCover11, Coordinates)
extracted_values_new_d<- raster::extract(LandCover13, Coordinates)

# Create a new data frame with extracted values
Coordinates_values_new <- data.frame(
  Coordinates,
  Raster_Value = extracted_values_new,
  Raster_Value1 = extracted_values_new_b,
  Raster_Value2 = extracted_values_new_c,
  Raster_Value3 = extracted_values_new_d
)

write.csv(Coordinates_values_new, "\\Coordinates_values_new.csv", row.names=FALSE)


####Viewing climatic data ####

precipitation_raster_june<- raster::raster("wc_precipitation/wc2.1_2.5m_prec_2021-06.tif")
precipitation_raster_july<- raster::raster("wc_precipitation/wc2.1_2.5m_prec_2021-07.tif")
tmin_raster_june<-raster::raster("wc-tmin/wc2.1_2.5m_tmin_2021-06.tif")
tmin_raster_july<- raster::raster("wc-tmin/wc2.1_2.5m_tmin_2021-07.tif")
tmax_raster_june<-raster::raster("wc-tmax/wc2.1_2.5m_tmax_2021-06.tif")
tmax_raster_july<- raster::raster("wc-tmax/wc2.1_2.5m_tmax_2021-07.tif")

precipitation_clim_june <- raster::extract(precipitation_raster_june, Coordinates)
precipitation_clim_july <- raster::extract(precipitation_raster_july, Coordinates)
tmin_clim_june <- raster::extract(tmin_raster_june, Coordinates)
tmin_clim_july <- raster::extract(tmin_raster_july, Coordinates)
tmax_clim_june <- raster::extract(tmax_raster_june, Coordinates)
tmax_clim_july <- raster::extract(tmax_raster_july, Coordinates)

# Create a new data frame with extracted values from june
Coordinates_values_clim_june <- data.frame(
  Coordinates,
  Precipitation = precipitation_clim_june,
  Tmin = tmin_clim_june,
  Tmax =tmax_clim_june
)

# Create a new data frame with extracted values from july
Coordinates_values_clim_july <- data.frame(
  Coordinates,
  Precipitation = precipitation_clim_july,
  Tmin = tmin_clim_july,
  Tmax =tmax_clim_july
)

write.csv(Coordinates_values_clim_june, "\\Coordinates_values_clim_june.csv", row.names=FALSE)
write.csv(Coordinates_values_clim_july, "\\Coordinates_values_clim_july.csv", row.names=FALSE)

