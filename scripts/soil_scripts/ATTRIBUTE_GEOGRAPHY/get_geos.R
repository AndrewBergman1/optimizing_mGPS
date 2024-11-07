#################################################

# Author: Andrew Bergman
# Course: BINP37, Optimizing mGPS 

################################################

# Description:
#   The script utilizes sf, sp and rnaturalearth to fetch geographic data for a dataset tagged with longitude and latitude.


# Imports 
library(sf)
library(sp)
library(rnaturalearth)
library(rnaturalearthdata)





# Takes longitude, latitude, country and continent. Country/Continent are the names of the columns to be created.
get_geography <- function(data, 
                          longitude_col = "Longitude", 
                          latitude_col = "Latitude",
                          country_col = "Country",
                          continent_col = "Continent") {
  
  # Load the country boundaries with continent and country information
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  
  # Check if necessary columns exist in the 'world' dataset
  required_world_columns <- c("continent", "name_en")
  missing_world_columns <- setdiff(required_world_columns, names(world))
  if (length(missing_world_columns) > 0) {
    stop(paste("The 'world' dataset is missing the following required columns:", 
               paste(missing_world_columns, collapse = ", ")))
  }
  
  # Validate input data has the required longitude and latitude columns
  if (!all(c(longitude_col, latitude_col) %in% names(data))) {
    stop(paste("Input data must contain", longitude_col, "and", latitude_col, "columns."))
  }
  
  # Convert the input data to an sf object using Longitude and Latitude
  data_sf <- st_as_sf(data, coords = c(longitude_col, latitude_col), crs = 4326, remove = FALSE)
  
  # Perform a spatial join to attribute continent and country to each point
  # Using st_intersects for robustness
  # Select only the necessary columns to optimize performance
  world_subset <- world[, c("continent", "name_en")]
  
  data_with_geography <- st_join(data_sf, world_subset, join = st_intersects)
  
  # Assign the continent and country to the original data
  data[[continent_col]] <- make.names(data_with_geography$continent)
  data[[country_col]] <- make.names(data_with_geography$name_en)
  
  # Ensure that the new columns are character vectors
  data[[continent_col]] <- as.character(data[[continent_col]])
  data[[country_col]] <- as.character(data[[country_col]])
  
  # Apply make.names to sanitize the Continent and Country names
  data[[continent_col]] <- make.names(data[[continent_col]])
  data[[country_col]] <- make.names(data[[country_col]])
  
  # Optionally, replace dots with underscores or another character if preferred
  # data[[continent_col]] <- gsub("\\.", "_", data[[continent_col]])
  # data[[country_col]] <- gsub("\\.", "_", data[[country_col]])
  
  # Optionally, remove the spatial geometry to return a regular dataframe
  data <- st_drop_geometry(data)
  
  return(data)
}


# Read dataset 
data <- read.csv("~/Master/15hp_project/mGPS-master/mGPS-master/soil_scripts/serious_time/2_inteprolation/iterations/2x_interpol/interpolated_dataset.csv")


# Attribute countries
data <- get_geography(data)

# Write dataset w/ countries
write_csv(data, "~/Master/15hp_project/mGPS-master/mGPS-master/soil_scripts/serious_time/3_attribute_geo/iterations/2x_interpol/interpolated_dataset.csv")
