set.seed(1)
# Load Required Libraries
library(caret)
library(geosphere)
library(dplyr)
library(xgboost)
library(e1071)
library(doParallel)
library(leaflet)
library(tibble)
library(smotefamily)
library(imbalance)
library(Boruta)
library(rworldmap)    # For coastlines
library(maps)         # For map.where
library(sp)           # For spatial distance calculations
library(htmlwidgets)  # For saving Leaflet maps
library(ggplot2)
library(reshape2)
library(tidyr)
library(sf)
library(sp)
library(rnaturalearth)
library(rnaturalearthdata)


# Define Helper Functions (as in your original script)
data_normalise <- function(df) {
  return(df / rowSums(df))
}

data_preprocess_f <- function(train_f, target_in, hierarchy, remove_small) {
  set.seed(1)
  metasub_data <- droplevels(train_f)
  
  if(length(hierarchy) == 4){
    metasub_data[, hierarchy[1]] <- factor(metasub_data[, hierarchy[1]])
    metasub_data[, hierarchy[2]] <- factor(metasub_data[, hierarchy[2]])
  } else {
    metasub_data[, hierarchy[1]] <- factor(metasub_data[, hierarchy[1]])
  }
  
  remove_small <- as.numeric(remove_small)
  if (remove_small > 0){
    small_cities <- names(which(summary(metasub_data[, target_in]) < remove_small))
    remove_samples <- which(metasub_data[, target_in] %in% c("antarctica", small_cities))
    if (length(remove_samples) != 0){
      metasub_data <- droplevels(metasub_data[-c(remove_samples), ])
    }
  }
  
  for (i in 1:length(hierarchy)){
    empty <- which(metasub_data[, hierarchy[i]] == "" | is.na(metasub_data[, hierarchy[i]]))
    if (length(empty) != 0){
      metasub_data <- metasub_data[-c(empty),]
    }
  }
  
  metasub_data[, hierarchy[1]] <- droplevels(metasub_data[, hierarchy[1]])
  if(length(hierarchy) == 4){
    metasub_data[, hierarchy[2]] <- droplevels(metasub_data[, hierarchy[2]])
  }
  
  print("Metadata has been pre-processed ...")
  return(metasub_data)
}

species_select <- function(x, y, remove_correlated = TRUE, cores = 1, maxRuns = 100, confidence = 0.01) {
  set.seed(1)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  if (remove_correlated) {
    cor_matrix <- cor(x, use = "complete.obs")
    if(any(is.na(cor_matrix))) {
      warning("NA values found in correlation matrix. Check for constant variables.")
      cor_matrix[is.na(cor_matrix)] <- 0
    }
    
    correlated <- findCorrelation(cor_matrix, cutoff = 0.98, verbose = TRUE, names = FALSE)
    
    print("Correlated features indices:")
    print(correlated)
    
    if (length(correlated) > 0) {
      x <- x[, -correlated, drop = FALSE]
      print(paste("Correlated features removed:", length(correlated)))
    } else {
      print("No correlated features to remove.")
    }
  }
  
  set.seed(123)
  boruta_output <- Boruta(x, y, maxRuns = maxRuns, doTrace = 2, ntree = 50)
  
  print(boruta_output)
  plot(boruta_output, las = 2)
  
  final_features <- getSelectedAttributes(boruta_output, withTentative = FALSE)
  print(paste("Selected features:", toString(final_features)))
  
  stopCluster(cl)
  registerDoParallel(1)
  
  return(list(results = boruta_output, selectedFeatures = final_features))
}

# Define the mGPS Function (as in your original script)
mGPS <- function(training = training_data,
                 testing = testing_data,
                 classTarget = "city",
                 hierarchy = c('continent', 'city', 'Latitude', 'Longitude'),
                 variables,
                 nthread = 1,
                 coast = NULL) {
  if (is.null(training)){
    return(message("No training set given"))
  } else {
    training <- droplevels(training)
    message("Training mGPS...")
    
    set.seed(1234)
    folds <- createFolds(training[, classTarget], k = 5, returnTrain = TRUE)
    
    trControlClass <- trainControl(
      method = "cv",
      number = 5,
      verboseIter = FALSE,
      returnData = FALSE,
      search = "grid",
      savePredictions = "final",
      classProbs = TRUE,
      allowParallel = TRUE,
      index = folds
    )
    
    trControl <- trainControl(
      method = "cv",
      number = 5,
      verboseIter = FALSE,
      returnData = FALSE,
      search = "grid",
      savePredictions = "final",
      allowParallel = TRUE,
      index = folds
    )
    
    tune_grid <- expand.grid(
      nrounds = c(300, 600),
      eta = c(0.05, 0.1),
      max_depth = c(3, 6, 9),
      gamma = 0,
      colsample_bytree = c(0.6, 0.8),
      min_child_weight = c(1),
      subsample = 0.7
    )
    
    if (length(hierarchy) == 4) {
      Xgb_region <- train(x = training[, variables], y = training[, hierarchy[1]],
                          method = "xgbTree",
                          trControl = trControlClass,
                          tuneGrid = tune_grid,
                          nthread = nthread)
      
      l1_train <- data.frame(training[, variables], 
                             Xgb_region[["pred"]][order(Xgb_region$pred$rowIndex), levels(training[, hierarchy[1]])])
    } else {
      l1_train <- training[, variables]
    }
    
    Xgb_class <- train(x = l1_train, y = training[, classTarget],
                       method = "xgbTree",
                       trControl = trControlClass,
                       tuneGrid = tune_grid,
                       nthread = nthread)
    
    l2_train <- data.frame(l1_train, 
                           Xgb_class[["pred"]][order(Xgb_class$pred$rowIndex), levels(training[, classTarget])])
    
    Xgb_latitude <- train(x = l2_train, y = training[, "Latitude"],
                          method = "xgbTree",
                          trControl = trControl,
                          tuneGrid = tune_grid,
                          nthread = nthread)
    
    l3_train <- data.frame(l2_train, 
                           "latPred" = Xgb_latitude[["pred"]][order(Xgb_latitude$pred$rowIndex), "pred"])
    
    Xgb_longitude <- train(x = l3_train, y = training[, "Longitude"],
                           method = "xgbTree",
                           trControl = trControl,
                           tuneGrid = tune_grid,
                           nthread = nthread)
  }
  
  if (is.null(testing)){
    model <- function(test, variables) {
      regProbs <- predict(Xgb_region, newdata = test[, variables], type = "prob")
      l1_test <- data.frame(test[, variables], regProbs)
      classPred <- predict(Xgb_class, newdata = l1_test)
      classProbs <- predict(Xgb_class, newdata = l1_test, type = "prob")
      l2_test <- data.frame(l1_test, classProbs)
      latPred <- predict(Xgb_latitude, newdata = l2_test)
      l3_test <- data.frame(l2_test, latPred)
      longPred <- predict(Xgb_longitude, newdata = l3_test)
      return(list(classPred, latPred, longPred))
    }
    message("No test set...returning trained mGPS model function")
    return(list(Xgb_region, Xgb_class, Xgb_latitude, Xgb_longitude, "model" = model))
  } else {
    message("Generating predictions")
    if (length(hierarchy) == 4) {
      regProbs <- predict(Xgb_region, newdata = testing[, variables], type = "prob")
      l1_test <- data.frame(testing[, variables], regProbs)
    } else {
      l1_test <- testing[, variables]
    }
    classPred <- predict(Xgb_class, newdata = l1_test)
    classProbs <- predict(Xgb_class, newdata = l1_test, type = "prob")
    l2_test <- data.frame(l1_test, classProbs)
    latPred <- predict(Xgb_latitude, newdata = l2_test)
    l3_test <- data.frame(l2_test, latPred)
    longPred <- predict(Xgb_longitude, newdata = l3_test)
    
    longPred[longPred > 180] <- 180
    longPred[longPred < -180] <- -180
    latPred[latPred > 90] <- 90
    latPred[latPred < -90] <- -90
    
    if (!is.null(coast)) {
      toAdjust <- which(is.na(maps::map.where(database = "world", longPred, latPred)))
      adjusted <- mapply(find_coast, long = longPred[toAdjust], lat = latPred[toAdjust])
      longPred[toAdjust] <- adjusted[1,]
      latPred[toAdjust] <- adjusted[2,]
    }
    
    return(list(classPred, latPred, longPred))
  }
}

haversine_distance <- function(lat1, lon1, lat2, lon2) {
  # Earth radius in kilometers
  R <- 6371
  
  # Convert latitude and longitude from degrees to radians
  lat1 <- lat1 * pi / 180
  lon1 <- lon1 * pi / 180
  lat2 <- lat2 * pi / 180
  lon2 <- lon2 * pi / 180
  
  # Haversine formula
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  distance <- R * c
  
  return(distance)  # Distance in kilometers
}

calculate_distances <- function(latitudes1, longitudes1, latitudes2, longitudes2) {
  distances <- haversine_distance(latitudes1, longitudes1, latitudes2, longitudes2)
  return(distances)
}

pull_land <- function(land_preds, hierarchy) {
  coastlines <- cbind("x" = maps::SpatialLines2map(rworldmap::coastsCoarse)$x, 
                      "y" = maps::SpatialLines2map(rworldmap::coastsCoarse)$y)
  coastlines <- coastlines[complete.cases(coastlines), ]
  coastlines <- coastlines[coastlines[, 1] < 180, ]
  
  find_coast <- function(long, lat) {
    distances_from_coastline <- sp::spDistsN1(coastlines, c(long, lat), longlat = TRUE)
    closest_point <- which.min(distances_from_coastline)
    new_coords <- coastlines[closest_point, ]
    return(new_coords)
  }
  
  toAdjust <- which(is.na(maps::map.where(database = "world", land_preds$longPred, land_preds$latPred)))
  
  print("Points to adjust:")
  print(toAdjust)
  
  if (length(toAdjust) > 0) {
    adjusted <- mapply(find_coast, long = land_preds$longPred[toAdjust], lat = land_preds$latPred[toAdjust])
    
    print("Adjusted coordinates:")
    print(adjusted)
    
    if (length(adjusted) != 0) {
      land_preds$longPred[toAdjust] <- adjusted[1, ]
      land_preds$latPred[toAdjust] <- adjusted[2, ]
    }
  }
  
  print(hierarchy)
  
  if (class(hierarchy) == 'character' && hierarchy[3] %in% colnames(land_preds)) {
    for (i in 1:nrow(land_preds)) {
      if (length(hierarchy) == 3) {
        land_preds[i, "Distance_from_origin"] <- geosphere::distm(
          c(land_preds[i, "longPred"], land_preds[i, "latPred"]), 
          c(land_preds[i, hierarchy[3]], land_preds[i, hierarchy[2]]), 
          fun = geosphere::distHaversine
        ) / 1000
      } else if (length(hierarchy) == 4) {
        land_preds[i, "Distance_from_origin"] <- geosphere::distm(
          c(land_preds[i, "longPred"], land_preds[i, "latPred"]), 
          c(land_preds[i, hierarchy[4]], land_preds[i, hierarchy[3]]), 
          fun = geosphere::distHaversine
        ) / 1000
      }
    }
  }
  
  return(land_preds)
}

run_mGPS_cv <- function(dataset, dataset_name) {
  cat("========================================\n")
  cat("Processing Dataset:", dataset_name, "\n")
  cat("========================================\n")
  
  # Preprocess the dataset
  dataset <- na.omit(dataset)
  
  set.seed(1234)  # For reproducibility
  folds <- createFolds(dataset$Country, k = 5, returnTrain = TRUE)
  
  # Initialize Vectors to Store Information from the Cross-Validation Loop
  errors <- c()
  classification_accuracies <- c()
  country_tags <- c()
  cm_values <- list()
  all_preds <- list(latPred = c(), longPred = c(), classPred = character())  # Initialize as character
  all_actuals <- c()  # To store actual labels for overall confusion matrix
  
  # Initialize a Data Frame to Store Distance-Based Metrics Per Fold
  distance_metrics_per_fold <- data.frame(
    Fold = integer(),
    Within_100km = integer(),
    Within_250km = integer(),
    Within_500km = integer(),
    Within_1000km = integer(),
    Within_1500km = integer(),
    Within_3000km = integer(),
    stringsAsFactors = FALSE
  )
  
  # Define Distance Thresholds
  distance_thresholds <- c(100, 250, 500, 1000, 1500, 3000)
  
  # Cross-Validation Loop
  for (i in 1:5) {
    cat("Processing Fold", i, "...\n")
    
    fold <- folds[[i]]
    train <- dataset[fold, ]
    test <- dataset[-fold, ]
    
    if ("Continent" %in% colnames(test)) {
      test$Continent <- as.factor(test$Continent)
      train$Continent <- as.factor(train$Continent)
    }
    test$Country <- as.factor(test$Country)
    train$Country <- as.factor(train$Country)
    
    # Feature Selection
    features <- species_select(
      x = train[!names(train) %in% c("Country", "Continent", "Latitude", "Longitude")],
      y = train$Country,
      maxRuns = 10000,
      cores = 4
    )
    
    variables <- features$selectedFeatures
    
    # Model Training and Prediction
    preds <- mGPS(
      training = train,
      testing = test,
      classTarget = "Country",
      nthread = 8,
      hierarchy = c("Continent", "Country", "Latitude", "Longitude"),
      variables = variables
    )
    
    preds <- list(
      classPred = preds[[1]],
      longPred = preds[[3]],
      latPred = preds[[2]]
    )
    
    # Pull Land and Ensure Predictions are Valid
    preds <- pull_land(preds, c("Continent", "Country", "Latitude", "Longitude"))
    
    pred_countries <- get_geography(preds)
    
    pred_countries <- as.factor(pred_countries$Country)
    pred_countries <- make.names(pred_countries)
    
    pred_countries <- ifelse(pred_countries == "NA.", NA, pred_countries)
    
    # Ensure Factor Levels Match
    pred_countries <- factor(pred_countries, levels = levels(test$Country))
    
    pred_countries <- as.factor(ifelse(is.na(as.character(pred_countries)), 
                                       "Not_represented_in_data", 
                                       as.character(pred_countries)))    
    
    pred_countries <- factor(pred_countries, levels = levels(test$Country))
    
    # Compute Confusion Matrix for the Current Fold
    cm <- confusionMatrix(preds$classPred, test$Country)
    
    # Store the Entire Confusion Matrix
    cm_values[[i]] <- cm
    
    # Calculate the Error Between the Predicted Coordinates and the Nearest Border of the Actual Country
    predicted_lat <- preds$latPred
    predicted_lon <- preds$longPred
    actual_lat <- test$Latitude
    actual_lon <- test$Longitude
    
    # Initialize a vector to store distances to country borders
    distance_to_border <- numeric(nrow(test))
    
    # Loop over each test instance
    for (j in 1:nrow(test)) {
      pred_lon <- predicted_lon[j]
      pred_lat <- predicted_lat[j]
      actual_country <- as.character(test$Country[j])
      
      # Compute distance to the border of the actual country
      distance_to_border[j] <- get_distance_to_country_border(pred_lon, pred_lat, actual_country)
    }
    
    # Convert distances from meters to kilometers
    distance_to_border_km <- distance_to_border / 1000  # Convert meters to kilometers
    
    # Set distance to 0 where classification is correct
    correct_classification <- pred_countries == test$Country
    distance_to_border_km[correct_classification] <- 0
    
    # Use distance_to_border_km as your error metric
    errors <- c(errors, distance_to_border_km)
    
    # Tag Errors with Their Corresponding Country
    country_tags <- c(country_tags, as.character(pred_countries))
    
    # Calculate Classification Accuracy
    num_matches <- sum(correct_classification)
    classification_accuracy <- num_matches / length(pred_countries)
    classification_accuracies <- c(classification_accuracies, classification_accuracy)
    
    # Store Predictions and Actuals for Overall Confusion Matrix
    all_preds$latPred <- c(all_preds$latPred, preds$latPred)
    all_preds$longPred <- c(all_preds$longPred, preds$longPred)
    all_preds$classPred <- c(all_preds$classPred, as.character(preds$classPred))  # Store as character
    all_actuals <- c(all_actuals, as.character(test$Country))
    
    # Compute Distance-Based Metrics for the Current Fold
    counts_within_thresholds <- sapply(distance_thresholds, function(threshold) sum(distance_to_border_km <= threshold))
    
    # Add the Counts to the Distance Metrics Data Frame
    distance_metrics_per_fold <- rbind(
      distance_metrics_per_fold,
      data.frame(
        Fold = i,
        Within_100km = counts_within_thresholds[1],
        Within_250km = counts_within_thresholds[2],
        Within_500km = counts_within_thresholds[3],
        Within_1000km = counts_within_thresholds[4],
        Within_1500km = counts_within_thresholds[5],
        Within_3000km = counts_within_thresholds[6],
        stringsAsFactors = FALSE
      )
    )
    
    # Print Per-Fold Confusion Matrix
    cat("Confusion Matrix for Fold", i, ":\n")
    print(cm$table)
    cat("\n")
  }
  
  # Calculate Overall Statistics
  mean_error <- mean(errors, na.rm = TRUE)
  median_error <- median(errors, na.rm = TRUE)
  min_error <- min(errors, na.rm = TRUE)
  max_error <- max(errors, na.rm = TRUE)
  mean_classification_accuracy <- mean(classification_accuracies)
  
  # Calculate Overall Distance-Based Metrics
  distance_counts_overall <- sapply(distance_thresholds, function(threshold) sum(errors <= threshold, na.rm = TRUE))
  names(distance_counts_overall) <- paste0("Within_", distance_thresholds, "km")
  
  distance_percentages_overall <- sapply(distance_thresholds, function(threshold) {
    (sum(errors <= threshold, na.rm = TRUE) / length(errors)) * 100
  })
  names(distance_percentages_overall) <- paste0("Within_", distance_thresholds, "km_%")
  
  # Output Cross-Validation Results
  cat("=== Cross-Validation Results for", dataset_name, "===\n")
  cat("Mean Distance Error (km):", mean_error, "\n")
  cat("Median Distance Error (km):", median_error, "\n")
  cat("Min Distance Error (km):", min_error, "\n")
  cat("Max Distance Error (km):", max_error, "\n")
  cat("Mean Classification Accuracy:", mean_classification_accuracy, "\n\n")
  
  # Compute and Print Overall Confusion Matrix
  overall_cm <- confusionMatrix(
    factor(all_preds$classPred, levels = levels(dataset$Country)),
    factor(all_actuals, levels = levels(dataset$Country))
  )
  
  cat("=== Overall Confusion Matrix for", dataset_name, "===\n")
  print(overall_cm$table)
  cat("\n")
  cat("=== Overall Metrics for", dataset_name, "===\n")
  print(overall_cm$overall)
  cat("\n")
  cat("=== Per-Class Metrics for", dataset_name, "===\n")
  print(overall_cm$byClass)
  cat("\n")
  
  # Compute per-class accuracies
  per_class_accuracy <- diag(overall_cm$table) / rowSums(overall_cm$table)
  
  # Create a data frame with class labels and per-class accuracies
  per_class_accuracy_df <- data.frame(
    Class = rownames(overall_cm$table),
    Accuracy = per_class_accuracy
  )
  
  # Plot per-class accuracies using ggplot2
  p <- ggplot(per_class_accuracy_df, aes(x = reorder(Class, Accuracy), y = Accuracy)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = paste("Per-Class Accuracy for", dataset_name),
         x = "Class",
         y = "Accuracy") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  # Display the plot
  print(p)
  
  # Save the plot as a PNG file, adjust height based on number of classes
  num_classes <- nrow(per_class_accuracy_df)
  plot_height <- 6 + 0.1 * num_classes
  
  ggsave(filename = paste0("per_class_accuracy_", dataset_name, ".png"), plot = p, width = 10, height = plot_height)
  
  # Output Distance-Based Metrics
  cat("=== Overall Distance-Based Metrics for", dataset_name, "===\n")
  print(as.data.frame(t(distance_counts_overall)))
  cat("\n")
  cat("=== Overall Distance-Based Metrics (Percentages) for", dataset_name, "===\n")
  print(as.data.frame(t(distance_percentages_overall)))
  cat("\n")
  
  # Save Confusion Matrices to CSV Files
  for (i in 1:5) {
    fold_cm <- cm_values[[i]]
    write.csv(fold_cm$table, paste0("confusion_matrix_", dataset_name, "_fold_", i, ".csv"), row.names = TRUE)
  }
  
  # Save Overall Confusion Matrix
  write.csv(overall_cm$table, paste0("confusion_matrix_", dataset_name, "_overall.csv"), row.names = TRUE)
  
  # Save Distance-Based Metrics to CSV Files
  write.csv(distance_metrics_per_fold, paste0("distance_metrics_", dataset_name, "_per_fold.csv"), row.names = FALSE)
  write.csv(data.frame(t(distance_counts_overall)), paste0("distance_metrics_", dataset_name, "_overall_counts.csv"), row.names = FALSE)
  write.csv(data.frame(t(distance_percentages_overall)), paste0("distance_metrics_", dataset_name, "_overall_percentages.csv"), row.names = FALSE)
  
  # Step 1: Create a data frame that contains both the original and predicted coordinates
  line_data <- data.frame(
    lng_orig = dataset$Longitude,
    lat_orig = dataset$Latitude,
    lng_pred = all_preds$longPred,
    lat_pred = all_preds$latPred
  )
  
  # Step 2: Calculate the distance between each pair of points (in meters)
  line_data$distance <- distHaversine(
    matrix(c(line_data$lng_orig, line_data$lat_orig), ncol = 2),
    matrix(c(line_data$lng_pred, line_data$lat_pred), ncol = 2)
  )
  
  # Step 3: Create a color palette function based on the distances
  pal <- colorNumeric(
    palette = "viridis",  # Color palette
    domain = line_data$distance
  )
  
  map <- leaflet() %>%
    addTiles() %>%
    addCircleMarkers(
      lng = dataset$Longitude,
      lat = dataset$Latitude,
      radius = 2,
      popup = "Original Data",
      color = "blue",
      stroke = TRUE,
      fill = FALSE,       # Makes the circles transparent (no fill)
      opacity = 0.5
    ) %>%
    addCircleMarkers(
      lng = all_preds$longPred,
      lat = all_preds$latPred,
      radius = 2,
      popup = "Predictions",
      color = "red",
      stroke = TRUE,
      fill = FALSE,       # Makes the circles transparent (no fill)
      opacity = 0.5
    ) %>%
    addControl(
      html = paste0("<h3>", dataset_name, "</h3>"),
      position = "topright"
    )
  
  # Display the map
  print(map)
  
  # Return a list of metrics
  return(list(
    mean_error = mean_error,
    median_error = median_error,
    min_error = min_error,
    max_error = max_error,
    mean_classification_accuracy = mean_classification_accuracy,
    distance_counts_overall = distance_counts_overall,
    distance_percentages_overall = distance_percentages_overall,
    overall_confusion_matrix = overall_cm
  ))
}


get_geography <- function(preds) {
  
  lon <- preds$longPred
  lat <- preds$latPred
  
  df <- data.frame(lon = lon, 
                   lat = lat)
  
  # Load the country boundaries with continent and country information
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  
  # Check if necessary columns exist in the 'world' dataset
  required_world_columns <- c("continent", "name_en")
  missing_world_columns <- setdiff(required_world_columns, names(world))
  if (length(missing_world_columns) > 0) {
    stop(paste("The 'world' dataset is missing the following required columns:", 
               paste(missing_world_columns, collapse = ", ")))
  }
  
  # Convert the input data to an sf object using Longitude and Latitude
  data_sf <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  
  # Perform a spatial join to attribute continent and country to each point
  world_subset <- world[, c("continent", "name_en")]
  
  data_with_geography <- st_join(data_sf, world_subset, join = st_intersects)
  
  # Assign the continent and country to the original data
  df[["Continent"]] <- make.names(data_with_geography$continent)
  df[["Country"]] <- make.names(data_with_geography$name_en)
  
  # Ensure that the new columns are character vectors
  df[["Continent"]] <- as.character(df[["Continent"]])
  df[["Country"]] <- as.character(df[["Country"]])
  
  # Optionally, remove the spatial geometry to return a regular dataframe
  df <- st_drop_geometry(df)
  
  return(df)
}

get_distance_to_country_border <- function(longitude, latitude, country_name) {
  # Get the country's geometry
  countries <- ne_countries(scale = "medium", returnclass = "sf")
  country_geom <- countries[countries$name == country_name, ]
  
  # Check if country was found
  if (nrow(country_geom) == 0) {
    warning(paste("Country not found:", country_name))
    return(NA)
  }
  
  # Create an sf point from the predicted coordinates
  point <- st_sfc(st_point(c(longitude, latitude)), crs = st_crs(4326))
  
  # Transform geometries to a projected CRS (e.g., EPSG:3857) for accurate distance measurement
  country_geom_proj <- st_transform(country_geom, crs = 3857)
  point_proj <- st_transform(point, crs = 3857)
  
  # Compute the distance from the point to the country's border
  # First, check if point is inside the country
  inside_country <- st_within(point_proj, country_geom_proj, sparse = FALSE)
  
  if (inside_country) {
    # If point is inside the country, distance is 0
    return(0)
  } else {
    # If point is outside the country, compute the distance to the country's border
    distance <- st_distance(point_proj, country_geom_proj)
    return(as.numeric(distance))  # Convert units to numeric (meters)
  }
}

# Load Your Datasets
interpolation_smote <- read.csv("~/Master/15hp_project/mGPS-master/mGPS-master/soil_scripts/serious_time/1_SMOTE/full_dataset/SMOTEd_soil_samples.csv")
smote <- read.csv("~/Master/15hp_project/mGPS-master/mGPS-master/soil_scripts/serious_time/2_inteprolation/full_dataset/interpolated_SMOTEd_soil_samples.csv")
original <- read.csv("C:/Users/andre/OneDrive/Documents/Master/15hp_project/mGPS-master/mGPS-master/Data/Data/Soil/combined.csv")

interpolation_smote$Continent <- as.factor(interpolation_smote$Continent)
smote$Continent <- as.factor(smote$Continent)
original$Continent <- as.factor(original$Continent)
interpolation_smote$Country <- as.factor(interpolation_smote$Country)
smote$Country <- as.factor(smote$Country)
original$Country <- as.factor(original$Country)


# Define a list of datasets to process
datasets <- list(
  interpolation_smote = interpolation_smote,
  smote = smote,
  original = original
)

# Initialize a list to store results
results_list <- list()

# Process each dataset
for (dataset_name in names(datasets)) {
  dataset <- datasets[[dataset_name]]
  
  # Run cross-validation and collect metrics
  results <- run_mGPS_cv(dataset, dataset_name)
  
  # Store the results
  results_list[[dataset_name]] <- results
}

# Compare Performance Metrics between 'data' and 'control'

# Example: Compare Mean Classification Accuracy
cat("=== Comparison of Mean Classification Accuracy ===\n")
accuracy_comparison <- data.frame(
  Dataset = names(results_list),
  Mean_Classification_Accuracy = sapply(results_list, function(x) x$mean_classification_accuracy)
)
print(accuracy_comparison)
cat("\n")

# Example: Compare Mean Distance Error
cat("=== Comparison of Mean Distance Error (km) ===\n")
error_comparison <- data.frame(
  Dataset = names(results_list),
  Mean_Distance_Error_km = sapply(results_list, function(x) x$mean_error)
)
print(error_comparison)
cat("\n")

# Example: Compare Distance-Based Metrics Percentages
cat("=== Comparison of Distance-Based Metrics (Percentages) ===\n")
distance_percentage_comparison <- data.frame(
  Dataset = character(),
  Within_100km = numeric(),
  Within_250km = numeric(),
  Within_500km = numeric(),
  Within_1000km = numeric(),
  Within_1500km = numeric(),
  Within_3000km = numeric(),
  stringsAsFactors = FALSE
)

for (dataset_name in names(results_list)) {
  distance_percentage_comparison <- rbind(
    distance_percentage_comparison,
    c(Dataset = dataset_name, 
      t(results_list[[dataset_name]]$distance_percentages_overall))
  )
}

# Convert columns to numeric except the first
distance_percentage_comparison[, -1] <- sapply(distance_percentage_comparison[, -1], as.numeric)

print(distance_percentage_comparison)
cat("\n")

names(distance_percentage_comparison) <- c("data_type", "<100km", "<250km", "<500km", "<1000km", "<1500km", "<3000km")

long_dist <- pivot_longer(distance_percentage_comparison, cols = -data_type, values_to = "Value", names_to = "Metric")
long_dist$Metric <- factor(long_dist$Metric, levels = c("<100km", "<250km", "<500km", "<1000km", "<1500km", "<3000km"))

ggplot(long_dist, aes(x = Metric, y = Value, col = data_type)) +
  geom_point() +
  geom_line(aes(group = data_type)) +
  labs(y = "Percent", title = "SMOTE + Interpolation vs SMOTE")

ggplot(error_comparison, aes(x = Dataset, y = Mean_Distance_Error_km)) +
  geom_line(group = "1") +
  geom_point()




