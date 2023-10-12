library(arcgisbinding)
library(here)
library(tidyverse)
library(raster)
library(rgdal)
library(sf)
library(maxlike)

# Initialize R-ArcGIS Bridge
arc.check_product()

# Set file path of species location data
location_points_filepath <- here::here('AfricanBuffalo', 'AfricanBuffalo.gdb', 'African_Buffalo_Locations')

# Convert arc shape to spatial data frame 
location_points <- arc.select(arc.open(location_points_filepath), sr = 3857)

# Extract coordinates
shape_info <- arc.shape(location_points)
coordinates <- cbind(shape_info$x, shape_info$y)

# Set file path of environmental parameters raster data
env_rast_filepath <- here::here('AfricanBuffalo', 'Focal_Statistics_Results.gdb')

# Read environmental parameters raster data
env_rast <- arc.open(env_rast_filepath)

# Extract environmental parameters raster names
env_rast_list <- env_rast@children$RasterDataset

# Raster stack building ----
# Initialize raster ID
raster_id <- 0

# Define empty raster list
raster_list <- list()

# Define empty raster names list
raster_name_list <- list()

# Populate environmental parameters raster data in list
for (i in env_rast_list){
  
  # Set file path of raster
  raster_filepath <- here::here('AfricanBuffalo', 'Focal_Statistics_Results.gdb', i)
  
  # Read raster as arc raster
  raster <- arc.raster(arc.open(raster_filepath))
  
  # Convert arc raster to raster
  raster_terra <- as.raster(raster)
  
  # Compute sum of pixels in raster
  pixel_count <- cellStats(raster_terra > 0, 'sum')
  
  # Define conditional statement,
  # where add raster to list if variation exists
  # compared to the number of location points for buffalo
  
  # If there are 100 buffalo, then at least 100 cells 
  # should contain non-zero values
  
  if(pixel_count >= dim(location_points)[1]){
    
    # Initialize raster ID count 
    raster_id <- raster_id + 1
    
    # Create standardized (normalized) raster
    raster_norm <- (raster_terra - cellStats(raster_terra, 'mean') / cellStats(raster_terra, 'sd'))
    
    # Add to raster list
    raster_list[[raster_id]] <- raster_norm
    
    # Add to raster names list
    raster_name_list[raster_id] <- i
    
  }
}

# Convert list to raster stack
raster_stack <- stack(raster_list)
names(raster_stack) <- raster_name_list

# MaxEnt model building ----
# Define mathematical expression
math_exp_str <- paste(raster_name_list, collapse = '+')
math_exp <- as.formula(paste('~', math_exp_str))

# Fit model to location points and environmental parameters
maxent_model <- maxlike(math_exp, raster_stack, coordinates, 
                        link = 'logit', hessian = TRUE, 
                        removeDuplicates = TRUE,
                        savedata = TRUE)

# Plot confidence interval for model coefficients
confint(maxent_model)

# Predict buffalo occurrence probability for study area
suitability_map <- predict(maxent_model)

# View buffalo occurrence probability histogram
hist(suitability_map)

# View buffalo occurrence probability raster
plot(suitability_map)
# Add location points to view
points(coordinates[, 1], coordinates[, 2])

# Set file path of export
export_filepath <- here::here('AfricanBuffalo', 'AfricanBuffalo.gdb', 'suitability')

# Export buffalo occurrence probability raster
arc.write(export_filepath, suitability_map)