##Bioclimatic index

#Create mixed model formula
formula_index <- formula(Elevation ~ Tmin+Tmax+Precipitation+World_cover_2021_class+Copernicus_cover_2018_class+Copernicus_tree_cover_2018_code+(1|Transect))

#Fit the mixed model
model_index <- lmer(formula_index, data = Pyrenees_metadata_all)

#Print the model summary
summary(model_index)

# Predict Elevation using the linear mixed model
Pyrenees_metadata_all$Predicted_Elevation <- predict(model_index, newdata = Pyrenees_metadata_all)

# Assess the model by calculating the R^2 and plotting the predicted elevation values
r_squared <- cor(Pyrenees_metadata_all$Elevation, Pyrenees_metadata_all$Predicted_Elevation)^2
print(r_squared)

plot(Pyrenees_metadata_all$Elevation, Pyrenees_metadata_all$Predicted_Elevation, 
     xlab = "Observed Elevation", ylab = "Predicted Elevation")+
    abline(0, 1, col = "red")

write.csv(Pyrenees_metadata_all, "\\Pyrenees_metadata_all.csv", row.names=FALSE)
