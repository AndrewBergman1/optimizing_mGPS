#################################################

# Author: Andrew Bergman
# Course: BINP37, Optimizing mGPS 

################################################

# Description:
#   The script finds the optimal number of clusters (2-4) and divides countries into that many clusters. This is upstream of running a clustered mGPS.



# Load necessary libraries
library(tidyverse)
library(cluster)
library(ClusterR)  
library(tidyr)

# Function to perform k-means clustering
kmeans_cluster <- function(data, k) {
  clusters <- kmeans(data, k, nstart = 20) 
  return(clusters)
}

# Function to find the optimal number of clusters using silhouette analysis
find_optimal_k <- function(data) {
  # Ensure the data has more than one unique point to avoid k-means issues
  if (nrow(data) < 2) {
    print("Not enough data for clustering, returning k = 1")
    return(1)
  }
  
  tryCatch({
    sil_width <- numeric(4)  # Initialize a numeric vector with 4 elements
    
    # Evaluate 2 to 4 clusters
    for (k in 2:4) {
      clusters <- kmeans_cluster(data, k)
      sil <- silhouette(clusters$cluster, dist(data))
      sil_width[k] <- mean(sil[, 3])  # Store the average silhouette width for each k
    }
    
    best_k <- which.max(sil_width[2:4]) + 1  # Get the best k (adjust by +1)
    
    print(paste("Best number of clusters for current country is:", best_k))
    return(best_k)
  }, error = function(e) {
    print(paste("Error during clustering: ", e$message))
    return(1)  # Default to 1 cluster in case of error
  })
}

# Function to add clustering information to the data with unique country-specific cluster labels
add_clustering_information <- function(data, clusters, country) {
  # Create unique cluster labels by prepending the country name to the cluster number
  data$Cluster <- paste(country, clusters$cluster, sep = "_")
  return(data)
}

# Set seed for reproducibility
set.seed(123)

# Read the CSV file
data <- read.csv("~/Master/15hp_project/mGPS-master/mGPS-master/soil_scripts/serious_time/3_attribute_geo/round_2/original_dataset.csv")

countries <- unique(data$Country)

# Initialize an empty data frame to store all country-specific clustered data
clustered_data <- data.frame()

# Loop through each country
for (country in countries) {
  print(paste("Processing", country))
  
  # Subset data for the specific country
  country_specific_data <- data[data$Country == country, ]
  
  # Separate non-coordinate columns and coordinates
  non_coordinates <- country_specific_data[!names(country_specific_data) %in% c("Longitude", "Latitude")]
  coordinates <- country_specific_data[names(country_specific_data) %in% c("Longitude", "Latitude")]
  
  # Find the optimal number of clusters
  best_k <- find_optimal_k(coordinates)
  
  # Perform k-means clustering
  clusters <- kmeans_cluster(coordinates, best_k)
  
  # Add clustering information with unique country-specific cluster names
  coordinates_with_clusters <- add_clustering_information(coordinates, clusters, country)
  
  # Combine the non-coordinate columns back with the clustered data
  country_clustered_data <- cbind(non_coordinates, coordinates_with_clusters)
  
  # Append the country's clustered data to the final dataframe
  clustered_data <- rbind(clustered_data, country_clustered_data)
}

write.csv(clustered_data, "~/Master/15hp_project/mGPS-master/mGPS-master/soil_scripts/serious_time/3.5_cluster_data/round_2/original_dataset.csviginal_dataset.csv")


