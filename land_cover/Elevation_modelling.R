### Elevation raster

#Training the model with the elevation extracted from another elevation raster
LandCover14 <- raster::raster("euro-dem-tif/data/eurodem/eurodem.tif")

# Define the extent representing the study area
working_polygon <- setExtent(LandCover14, extent(c(-1, 1, 42.3, 43)))

e <- as(extent(-1, 1, 42.3, 43), 'SpatialPolygons')
crs(e) <- "+proj=longlat +datum=WGS84 +no_defs"
working_polygon <- crop(LandCover14, e)

# Sample random points within the study area
extracted_values_training <- raster::sampleRandom(
  working_polygon,
  size = 100,
  na.rm = TRUE,
  cells = FALSE,
  xy = TRUE,
  sp = FALSE,
  asRaster = FALSE
)

coordinates_training<-extracted_values_training[,c(1,2)]

## France raster
LandCover15 <- raster::raster("DTM_France_20m_v3_by_Sonny.tif")
ext_interest <- extent(0, 500000, 4600000, 4900000) #select the interested raster
LandCover15_cropped <- crop(LandCover15, ext_interest)
crs_9 <- crs(LandCover9) #LandCover15 doesn't have the same crs so needs to be re-projected
LandCover15_reprojected <- projectRaster(LandCover15_cropped, crs = crs_9)
coordinates_training<-extracted_values_training[,c(1,2)]
extracted_values_france <- raster::extract(LandCover15_reprojected, coordinates_training)

##Spain raster
LandCover16 <- raster::raster("DTM_Spain_Mainland_20m_v2_by_Sonny.tif")
ext_interest_b <- extent(500000, 1000000, 4650000, 4800000) #select the interested raster
LandCover16_cropped <- crop(LandCover16, ext_interest_b)#LandCover16 doesn't have the same crs so needs to be re-projected
crs_9 <- crs(LandCover9)
LandCover16_reprojected <- projectRaster(LandCover16_cropped, crs = crs_9)
extracted_values_spain <- raster::extract(LandCover16_reprojected, coordinates_training)

Coordinates_values_raster <- data.frame(
  coordinates_training,
  france_elevation = extracted_values_france,
  spain_elevation = extracted_values_spain
)

write.csv(Coordinates_values_raster, "\\Coordinates_values_raster.csv", row.names=FALSE) #merge elevation values obtain in each of the raster objects. Add the elevation to the previos file with the rest of the metadata

# Modelling
formula_index_training <- formula(Elevation ~ Tmin+Tmax+Precipitation+World_cover_2021_code+Copernicus_cover_2018_code+Copernicus_tree_cover_2018_code)

model_index_training <- lm(formula_index_training, data = Coordinates_values_training)

summary(model_index)

# Predict Elevation using the linear mixed model
Coordinates_values_training$Predicted_Elevation <- predict(model_index_training, newdata = Coordinates_values_training)

# Assess the model by calculating the R^2 and plotting the predicted elevation values
r_squared <- cor(Coordinates_values_training$Elevation, Coordinates_values_training$Predicted_Elevation)^2
print(r_squared)

plot(Coordinates_values_training$Elevation, Coordinates_values_training$Predicted_Elevation, 
     xlab = "Observed Elevation", ylab = "Predicted Elevation")+
  abline(0, 1, col = "red")