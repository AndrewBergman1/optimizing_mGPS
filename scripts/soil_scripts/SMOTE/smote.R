#################################################

# Author: Andrew Bergman
# Course: BINP37, Optimizing mGPS 

################################################

# Description:
#   Perform SMOTE on a dataset.


# SMOTE the soil dataset in accordance to Yali's efforts
set.seed(1)
library(tidyverse)
library(UBL)
library(leaflet)
library(caret)
library(FactoMineR)
library(ggfortify)
library(factoextra)
library(scatterplot3d)
library(rgeoda)
library(geodaData)
library(cluster)
library(sf)
library(sp)
library(rnaturalearth)
library(rnaturalearthdata)
library(Boruta)
library(rsample)

# load data
load_data <- function(meta_path, taxa_path) {
  set.seed(1)
  meta_data <- read_csv("~/Master/15hp_project/mGPS-master/mGPS-master/Data/Data/Soil/Dataset_01_22_2018_enviro.csv")
  abundance_data <- read_csv("~/Master/15hp_project/mGPS-master/mGPS-master/Data/Data/Soil/Dataset_01_22_2018_taxa.csv")
  soil_data <- merge(meta_data, abundance_data, on="ID_Environmental")
  otu_columns <- grep("^OTU", names(soil_data))
  soil_otu_data <- cbind(soil_data[, as.numeric(otu_columns)], soil_data[c("Continent", "Longitude", "Latitude")])
  soil_otu_data$Continent <- as.factor(soil_otu_data$Continent)
  return(soil_otu_data)
}

#normalise data
data_normalise <- function(df) {
  return(df / rowSums(df))
}

#pre-process data
preprocess_data <- function(spatial_data) {
  set.seed(1)
  # Create world object and attribute countries (continents are already attributed)
  world <- ne_countries(scale = "medium", returnclass = "sf")
  world <- st_transform(world, crs = 4326)  # Ensure CRS is WGS 84
  
  # Convert spatial_data to sf object
  combined_sf <- st_as_sf(spatial_data, coords = c("Longitude", "Latitude"), crs = 4326)
  
  # Perform spatial join between combined_sf and world to get country information
  combined_sf <- st_join(combined_sf, world[, c("name_en", "continent")], join = st_intersects)
  
  # Extract coordinates back into Longitude and Latitude columns
  coords <- st_coordinates(combined_sf)
  combined_sf$Longitude <- coords[, 1] 
  combined_sf$Latitude <- coords[, 2]   
  
  # Convert the sf object back to a data frame
  combined_sf <- data.frame(combined_sf)
  
  # Rename and factorize country
  combined_sf$Country <- as.factor(combined_sf$name_en)
  combined_sf$Continent <- as.factor(combined_sf$continent)
  
  # Drop unnecessary columns
  combined_sf <- combined_sf[ , !names(combined_sf) %in% c("name_en", "continent", "geometry")]
  
  # make the factor names valid for mGPS application further on
  combined_sf$Country <- make.names(combined_sf$Country)
  combined_sf$Continent <- make.names(combined_sf$Continent)
  combined_sf$Country <- as.factor(combined_sf$Country)
  combined_sf$Continent <- as.factor(combined_sf$Continent)
  
  combined_sf <- combined_sf[complete.cases(combined_sf[, c("Continent", "Country")]), ]
  
  country_counts <- table(combined_sf$Country)
  combined_sf <- combined_sf[combined_sf$Country %in% names(country_counts[country_counts > 3]), ]
  combined_sf$Country <- droplevels(combined_sf$Country)
  otu_columns <- combined_sf[, grepl("^OTU", names(combined_sf))]
  normalized_otu_columns <- data_normalise(otu_columns)
  combined_sf[, grepl("^OTU", names(combined_sf))] <- normalized_otu_columns  
  return(combined_sf)
}
#perform smote
smote_data <- function(data, cluster_data = TRUE, smote_type, class_label) {
  set.seed(1)
  
  # Ensure that the class_label is a factor
  data[[class_label]] <- as.factor(data[[class_label]])
  
  # Select only numeric columns
  data_numeric <- data[, sapply(data, is.numeric)]  # Select only numeric columns
  
  # Add class_label back for SMOTE processing
  data_numeric[[class_label]] <- data[[class_label]]
  
  # Apply SMOTE to balance the 'class_label' class
  smoted_soil_data <- tryCatch({
    UBL::SmoteClassif(
      form = as.formula(paste(class_label, "~ .")), 
      dat = data_numeric,
      C.perc = smote_type,                                          
      k = 3,  
      repl = FALSE,
      dist = "Euclidean",
      p = 2
    )
  }, error = function(e) {
    message("Error during SMOTE: ", e)
    return(NULL)
  })
  
  # Check if SMOTE returned a valid result
  if (is.null(smoted_soil_data)) {
    return(data)  # If SMOTE fails, return the original data
  }
  
  # Use the get_continents function to assign continents based on Long/Lat for SMOTEd data
  smoted_soil_data <- get_continents(smoted_soil_data)
  
  # Combine the original and SMOTEd data
  combined_data <- smoted_soil_data
  combined_data$Continent <- as.factor(combined_data$Continent)
  
  return(combined_data)
}
#plot continental fractions of dataset
plot_continental_fractions <- function(smoted_otu_data, cluster_data = T, plot = F) {
  
  set.seed(1)
  
  if (plot == T) {
    if (cluster_data == T) {
      # Plot fractions by cluster
      continental_fractions <- data.frame(table(smoted_otu_data$Cluster))
      total_samples <- sum(continental_fractions$Freq)
      
      frac <- ggplot(continental_fractions, aes(x = "", y = Freq, fill = Var1)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        geom_text(aes(label = paste0(Var1, " (", Freq, ")")), 
                  position = position_stack(vjust = 0.5)) +
        labs(
          title = "Cluster-based Continental Fractions", 
          x = NULL, 
          y = "Frequency",
          caption = paste("Total Samples:", total_samples)
        ) +
        theme_minimal()
      
      print(frac)
    } else {
      # Plot fractions by country
      continental_fractions <- data.frame(table(smoted_otu_data$Country))
      total_samples <- sum(continental_fractions$Freq)
      
      frac <- ggplot(continental_fractions, aes(x = "", y = Freq, fill = Var1)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        geom_text(aes(label = paste0(Var1, " (", Freq, ")")), 
                  position = position_stack(vjust = 0.5)) +
        labs(
          title = "Country-based Continental Fractions", 
          x = NULL, 
          y = "Frequency",
          caption = paste("Total Samples:", total_samples)
        ) +
        theme_minimal()
      
      print(frac)
    }
  }
}

#generate map plot showing all samples
plot_all_coordinates <- function(soil_otu_data, smoted_otu_data, plot = F) {
  
  set.seed(1)
  
  if (plot == T) {
    map <- leaflet() %>%
      addTiles() %>%
      # Original Data (larger red circle markers)
      addCircleMarkers(
        lng = soil_otu_data$Longitude, 
        lat = soil_otu_data$Latitude, 
        radius = 5,          # Larger red markers
        color = "red", 
        fillColor = "red", 
        fillOpacity = 1,
        stroke = FALSE
      ) %>%
      # SMOTEd Data (smaller blue circle markers)
      addCircleMarkers(
        lng = smoted_otu_data$Longitude, 
        lat = smoted_otu_data$Latitude, 
        radius = 3,          # Smaller blue markers
        color = "blue", 
        fillColor = "blue", 
        fillOpacity = 0.5, 
        stroke = FALSE
      ) %>%
      addLegend(
        position = "bottomright", 
        colors = c("red", "blue"), 
        labels = c("Original Data", "SMOTEd Data"), 
        opacity = 1
      )
    
    map
  }
}

# Generate pca plot for samples' features
plot_pca <- function(data, label = "Cluster", plot = F) {
  
  set.seed(1)
  
  # Select only numeric columns for PCA, excluding non-numeric ones like 'Continent', 'Country', etc.
  numeric_data <- data[, sapply(data, is.numeric) & !names(data) %in% c("Longitude", "Latitude")]
  
  # Perform PCA on numeric columns
  pca <- prcomp(numeric_data, scale. = TRUE)
  
  # Automatically plot the PCA results using autoplot
  auto <- autoplot(pca, data = data, colour = label, frame = T) +
    labs(title = paste("PCA plot of ", label)) +
    theme_minimal() 
  
  # Create a data frame for PCA results
  PC <- as.data.frame(pca$x)
  
  # Dynamically add the label column
  PC$label <- data[[label]]
  
  # Calculate the percentage of variance explained by each component
  explained_variance <- pca$sdev^2 / sum(pca$sdev^2) * 100
  variance_df <- data.frame(PC = paste0("PC", 1:length(explained_variance)),
                            Variance = explained_variance)
  
  # Plot eigenvalues (scree plot)
  eig <- fviz_eig(pca, addlabels = TRUE, ylim = c(0, 70))
  
  if (plot == T) {
    print(auto)
    print(eig)
  }
  
  # Calculate Silhouette Scores using all principal components
  # Exclude the 'label' column when computing distances
  pc_columns <- grep("^PC", names(PC), value = TRUE)
  silhouette_scores <- silhouette(as.numeric(as.factor(PC$label)), dist(PC[, pc_columns]))
  
  # Calculate the mean silhouette score
  mean_sil <- mean(silhouette_scores[, 3])
  
  # Print the mean silhouette score
  print(paste("Mean Silhouette Score using all PCs:", round(mean_sil, 3)))
  
  return(mean_sil)
}

# Get continents using rnaturalearth
get_continents <- function(data) {
  
  # Load the continent boundaries from rnaturalearth
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  
  # Convert the input data to an sf object using its longitude and latitude
  data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)
  
  # Perform a spatial join to attribute the continent to each point
  data_with_continents <- st_join(data_sf, world["continent"], join = st_within)
  
  # Extract the continent name from the join
  data$Continent <- data_with_continents$continent
  data$Continent <- make.names(data$Continent)
  return(data)
}

#Compare smoted daata vs original data
test_smote_classif_performance <- function(soil_otu_data, class_label) {
  
  # Initialize list to store confusion matrices
  all_confusion_matrices <- list()
  
  # Create train-test folds using the original data (soil_otu_data)
  folds <- createFolds(soil_otu_data[[class_label]], k = 5, returnTrain = TRUE)
  
  # Loop over balance types
  for (i in c("balance", "extreme")) {
    
    # Loop through the datasets: smoted_otu_data and soil_otu_data
    for (j_name in c("smoted_otu_data", "soil_otu_data")) {
      
      # Loop over the folds for cross-validation
      for (k in seq_along(folds)) {
        # Get the train and test indices for the current fold
        train_indices <- folds[[k]]
        test_indices <- setdiff(seq_len(nrow(soil_otu_data)), train_indices)
        
        # Subset the training and testing data from the original soil_otu_data
        train_data <- soil_otu_data[train_indices, ]
        test_data <- soil_otu_data[test_indices, ]
        
        # If we're working with smoted_otu_data, apply SMOTE to the training data
        if (j_name == "smoted_otu_data") {
          train_data <- smote_data(train_data, cluster_data = F, smote_type = i, class_label = class_label)
        }
        
        # Remove 'Continent' or 'Country' depending on class_label
        if (class_label == "Country") {
          train_data <- train_data[, !names(train_data) %in% c("Continent", "Longitude", "Latitude")]
          test_data <- test_data[, !names(test_data) %in% c("Continent", "Longitude", "Latitude")]
        } else {
          train_data <- train_data[, !names(train_data) %in% c("Country", "Longitude", "Latitude")]
          test_data <- test_data[, !names(test_data) %in% c("Country", "Longitude", "Latitude")]
        }
        
        # Apply Boruta feature selection on the training set
        formula <- as.formula(paste(class_label, "~ ."))
        boruta_result <- Boruta(formula, data = train_data, doTrace = 1, maxRuns = 100)
        features <- getSelectedAttributes(boruta_result, withTentative = TRUE)
        
        # Ensure features are numeric
        numeric_features <- names(train_data)[sapply(train_data, is.numeric) & names(train_data) %in% features]
        
        # Train a random forest model using the selected numeric features
        model <- randomForest(as.formula(paste(class_label, "~ .")), data = train_data[, c(numeric_features, class_label)])
        
        # Predict the test data using the trained model
        predictions <- predict(model, test_data[, numeric_features])
        
        # Align factor levels between predictions and actuals for confusion matrix
        test_data[[class_label]] <- factor(test_data[[class_label]], levels = levels(predictions))
        
        # Generate the confusion matrix for this fold
        cm <- confusionMatrix(predictions, test_data[[class_label]])
        
        # Store the confusion matrix in the list
        all_confusion_matrices[[paste(i, j_name, "Fold", k, sep = "_")]] <- cm
        
        # Print progress
        print(paste(i, j_name, "Fold", k, "completed."))
      }
    }
  }
  
  # Initialize matrices to store the results
  balance_smote <- matrix(0, nrow = length(levels(soil_otu_data[[class_label]])), ncol = 11)
  extreme_smote <- matrix(0, nrow = length(levels(soil_otu_data[[class_label]])), ncol = 11)
  balance_original <- matrix(0, nrow = length(levels(soil_otu_data[[class_label]])), ncol = 11)
  extreme_original <- matrix(0, nrow = length(levels(soil_otu_data[[class_label]])), ncol = 11)
  
  # Aggregating confusion matrix results
  for (fold in c("balance_smoted_otu_data_Fold_1", "balance_smoted_otu_data_Fold_2", "balance_smoted_otu_data_Fold_3", "balance_smoted_otu_data_Fold_4", "balance_smoted_otu_data_Fold_5")) {
    balance_smote <- balance_smote + as.matrix(all_confusion_matrices[[fold]]$byClass)
  }
  for (fold in c("balance_soil_otu_data_Fold_1", "balance_soil_otu_data_Fold_2", "balance_soil_otu_data_Fold_3", "balance_soil_otu_data_Fold_4", "balance_smoted_otu_data_Fold_5")) {
    balance_original <- balance_original + as.matrix(all_confusion_matrices[[fold]]$byClass)
  }
  for (fold in c("extreme_soil_otu_data_Fold_1", "extreme_soil_otu_data_Fold_2", "extreme_soil_otu_data_Fold_3", "extreme_soil_otu_data_Fold_4", "extreme_smoted_otu_data_Fold_5")) {
    extreme_original <- extreme_original + as.matrix(all_confusion_matrices[[fold]]$byClass)
  }
  for (fold in c("extreme_smoted_otu_data_Fold_1", "extreme_smoted_otu_data_Fold_2", "extreme_smoted_otu_data_Fold_3", "extreme_smoted_otu_data_Fold_4", "extreme_smoted_otu_data_Fold_5")) {
    extreme_smote <- extreme_smote + as.matrix(all_confusion_matrices[[fold]]$byClass)
  }
  
  
  balance_smote <- balance_smote / 5
  extreme_smote <- extreme_smote / 5
  balance_original <- balance_original / 5
  extreme_original <- extreme_original / 5
  
  # Remove "Class: " from the rownames using gsub
  cleaned_rownames <- gsub("Class: ", "", rownames(extreme_smote))
  extreme_smote <- data.frame(class_label = cleaned_rownames, extreme_smote, row.names = NULL)
  balance_smote <- data.frame(class_label = cleaned_rownames, balance_smote, row.names = NULL)
  extreme_original <- data.frame(class_label = cleaned_rownames, extreme_original, row.names = NULL)
  balance_original <- data.frame(class_label = cleaned_rownames, balance_original, row.names = NULL)
  
  extreme_smote["Source"] <- rep("extreme_smote", nrow(extreme_smote))
  balance_smote["Source"] <- rep("balance_smote", nrow(balance_smote))
  extreme_original["Source"] <- rep("extreme_original", nrow(extreme_original))
  balance_original["Source"] <- rep("balance_original", nrow(balance_original))
  
  total_cm_by_class <- rbind(extreme_smote, balance_smote, extreme_original, balance_original)
  
  long_cm <- pivot_longer(total_cm_by_class, cols = names(total_cm_by_class[names(total_cm_by_class) %in% c("Balanced.Accuracy")]), values_to = "Values", names_to = "Metric")
  
  p <- ggplot(long_cm, aes(x = class_label, y = Values)) +
    geom_point() +
    geom_line(group = 1) +
    facet_wrap(~Source)
  
  print(p)
}

# If first test doesn't work, try this one
test_smote_classif_performance1 <- function(soil_otu_data, class_label, test) {
  
  # Initialize list to store confusion matrices
  all_confusion_matrices <- list()
  
  # Create train-test folds using the original data (soil_otu_data)
  folds <- createFolds(soil_otu_data[[class_label]], k = 5, returnTrain = TRUE)
    
  # Loop through the datasets: smoted_otu_data and soil_otu_data
  for (j_name in c("smoted_otu_data", "soil_otu_data")) {
      
    # Loop over the folds for cross-validation
    for (k in seq_along(folds)) {
      # Get the train and test indices for the current fold
      train_indices <- folds[[k]]
      test_indices <- setdiff(seq_len(nrow(soil_otu_data)), train_indices)
        
      # Subset the training and testing data from the original soil_otu_data
      train_data <- soil_otu_data[train_indices, ]
      test_data <- soil_otu_data[test_indices, ]
        
      # If we're working with smoted_otu_data, apply SMOTE to the training data
      if (j_name == "smoted_otu_data") {
        train_data <- smote_data(train_data, cluster_data = F, smote_type = test, class_label = class_label)
      }
        
      # Remove 'Continent' or 'Country' depending on class_label
      if (class_label == "Country") {
        train_data <- train_data[, !names(train_data) %in% c("Continent", "Longitude", "Latitude")]
        test_data <- test_data[, !names(test_data) %in% c("Continent", "Longitude", "Latitude")]
      } else {
        train_data <- train_data[, !names(train_data) %in% c("Country", "Longitude", "Latitude")]
        test_data <- test_data[, !names(test_data) %in% c("Country", "Longitude", "Latitude")]
      }
        
      # Apply Boruta feature selection on the training set
      formula <- as.formula(paste(class_label, "~ ."))
      boruta_result <- Boruta(formula, data = train_data, doTrace = 1, maxRuns = 10000)
      features <- getSelectedAttributes(boruta_result, withTentative = TRUE)
        
      # Ensure features are numeric
      numeric_features <- names(train_data)[sapply(train_data, is.numeric) & names(train_data) %in% features]
        
      # Train a random forest model using the selected numeric features
      model <- randomForest(as.formula(paste(class_label, "~ .")), data = train_data[, c(numeric_features, class_label)])
        
      # Predict the test data using the trained model
      predictions <- predict(model, test_data[, numeric_features])
        
      # Align factor levels between predictions and actuals for confusion matrix
      test_data[[class_label]] <- factor(test_data[[class_label]], levels = levels(predictions))
        
      # Generate the confusion matrix for this fold
      cm <- confusionMatrix(predictions, test_data[[class_label]])
        
      # Store the confusion matrix in the list
      all_confusion_matrices[[paste(j_name, "Fold", k, sep = "_")]] <- cm
        
      # Print progress
      print(paste(j_name, "Fold", k, "completed."))
    }
  }
  
  # Initialize matrices to store the results
  test_smote <- matrix(0, nrow = length(levels(soil_otu_data[[class_label]])), ncol = 11)
  original <- matrix(0, nrow = length(levels(soil_otu_data[[class_label]])), ncol = 11)

  
  # Aggregating confusion matrix results for "smoted_otu_data"
  for (fold in 1:5) {
    test_smote <- test_smote + as.matrix(all_confusion_matrices[[paste("smoted_otu_data_Fold", fold, sep = "_")]]$byClass)
  }
  
  # Aggregating confusion matrix results for "soil_otu_data"
  for (fold in 1:5) {
    original <- original + as.matrix(all_confusion_matrices[[paste("soil_otu_data_Fold", fold, sep = "_")]]$byClass)
  }
  
  test_smote <- test_smote / 5
  original <- original / 5
  
  # Remove "Class: " from the rownames using gsub
  cleaned_rownames <- gsub("Class: ", "", rownames(test_smote))
  test_smote <- data.frame(class_label = cleaned_rownames, test_smote, row.names = NULL)
  original <- data.frame(class_label = cleaned_rownames, original, row.names = NULL)
  
  test_smote["Source"] <- rep("Oversampling all countries to 10 samples (Yali's)", nrow(test_smote))
  original["Source"] <- rep("original", nrow(original))
  
  total_cm_by_class <- rbind(test_smote, original)
  
  long_cm <- pivot_longer(total_cm_by_class, cols = names(total_cm_by_class[names(total_cm_by_class) %in% c("Balanced.Accuracy")]), values_to = "Values", names_to = "Metric")
  
  p <- ggplot(long_cm, aes(x = class_label, y = Values)) +
    geom_point() +
    geom_line(group = 1) +
    facet_wrap(~Source)
  
  print(p)
}

# load meta data and taxa
meta_path <- "~/Master/15hp_project/mGPS-master/mGPS-master/Data/Data/Soil/Dataset_01_22_2018_enviro.csv"
taxa_path <- "~/Master/15hp_project/mGPS-master/mGPS-master/Data/Data/Soil/Dataset_01_22_2018_taxa.csv"

  
soil_otu_data <- read.csv("~/Master/15hp_project/mGPS-master/mGPS-master/soil_scripts/serious_time/3_attribute_geo/iterations/2x_interpol/interpolated_dataset.csv")


# Define test settings for SMOTE
test <- list(Argentina = 85/(5*2), Australia = 1, 
                Chile = 85/(7*2), Mexico = 85/(8*2), 
                Morocco = 85/(9*2), Panama = 85/(5*2),
                People.s.Republic.of.China = 85/(6*2), Spain = 85/(19*2), 
                Tunisia = 85/(6*2), United.Kingdom = 85/(11*2),
                United.States.of.America = 1)


# Assess "balanced" vs "Extreme" SMOTEing on the dataset (5x CV)
test_smote_classif_performance1(soil_otu_data, class_label = "Country", test = test)

# Perform the SMOTE that performed the best
smoted_otu_data <- smote_data(soil_otu_data, cluster_data = F, smote_type = test, class_label = "Country")

# Plot the sample dispersion on class-label
plot_continental_fractions(smoted_otu_data, cluster_data = F, plot = T)
plot_continental_fractions(soil_otu_data, cluster_data = F, plot = T)

# Plot the coordinates of SMOTEd and original samples
plot_all_coordinates(soil_otu_data, smoted_otu_data, plot = T)

plot_pca(smoted_otu_data, label = "Country", plot = T)
plot_pca(soil_otu_data, label = "Country", plot = T)

# Save smoted data
write_csv(smoted_otu_data, "~/Master/15hp_project/mGPS-master/mGPS-master/soil_scripts/serious_time/1_SMOTE/iterations/2x_interpol/interpolated_smoted_dataset.csv")
