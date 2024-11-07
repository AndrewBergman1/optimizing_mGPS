# ================================
# Spatial Interpolation Comparison
# IDW vs Multiple Spline and Other Methods with 5-fold Cross Validation and Grid Interpolation
# ================================

# ----------------------------
# 1. Load Necessary Libraries
# ----------------------------

set.seed(1)

required_packages <- c(
  "tidyverse", "gstat", "sp", "caret", "leaflet", 
  "sf", "gridExtra", "reshape2", "fields", "mgcv", "MBA", 
  "nngeo", "RANN", "geosphere"  # Added 'geosphere' for distance calculations
)
installed_packages <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Load libraries
library(tidyverse)   # Data manipulation
library(gstat)       # Geostatistical modeling and interpolation
library(sp)          # Spatial data classes
library(caret)       # Machine learning utilities
library(leaflet)     # Interactive maps
library(sf)          # Handling spatial data
library(gridExtra)   # Arranging multiple plots
library(reshape2)    # Data reshaping
library(fields)      # Thin Plate Splines and RBF
library(mgcv)        # Generalized Additive Models
library(MBA)         # Multilevel B-splines for Data Interpolation
library(nngeo)       # Nearest Neighbor Interpolation
library(RANN)        # Fast Nearest Neighbor Search
library(geosphere)   # Distance calculations

# ----------------------------
# 2. Function Definitions
# ----------------------------

# 2.1. Load and Preprocess Data
load_data <- function(path) {
  set.seed(1)
  
  if (!file.exists(path)) {
    stop("The specified file does not exist.")
  }
  
  data <- read_csv(path)
  
  # Identify OTU columns (assuming they start with "OTU")
  otu_indices <- grep("^OTU", names(data))
  if (length(otu_indices) == 0) {
    stop("No OTU columns found in the data.")
  }
  
  # Extract and convert OTU columns to numeric
  otu_data <- data[, otu_indices] %>%
    mutate(across(everything(), as.numeric))
  
  # Extract geographic data and convert types
  geo_data <- data %>%
    select(Longitude, Latitude, Country, Continent) %>%
    mutate(
      Latitude = as.numeric(Latitude),
      Longitude = as.numeric(Longitude),
      Continent = as.factor(Continent),
      Country = as.factor(Country)
    )
  
  # Combine geographic data with OTU data
  combined_data <- bind_cols(geo_data, otu_data)
  
  return(combined_data)
}

# 2.2. Transform Data to Spatial Object
transform_to_spatial_object <- function(data) {
  
  set.seed(1)
  
  # Remove rows with missing Longitude or Latitude
  complete_data <- data %>%
    drop_na(Longitude, Latitude)
  
  # Convert to SpatialPointsDataFrame
  sp_object <- SpatialPointsDataFrame(
    coords = complete_data %>% select(Longitude, Latitude),
    data = complete_data %>% select(-Longitude, -Latitude),
    proj4string = CRS("+proj=longlat +datum=WGS84")  # WGS 84 projection
  )
  
  return(sp_object)
}

# 2.3. IDW Interpolation Function
perform_idw <- function(train_data, test_points, target, idp = 2.0) {
  
  set.seed(1)
  
  idw_result <- idw(
    formula = as.formula(paste(target, "~ 1")),
    locations = train_data,
    newdata = test_points,
    idp = idp
  )
  return(idw_result$var1.pred)
}

# 2.4. Thin Plate Spline (TPS) Interpolation Function
perform_tps <- function(train_data, test_points, target) {
  
  set.seed(1)
  
  # Extract coordinates and target values
  coords <- coordinates(train_data)
  target_values <- train_data@data[[target]]
  
  # Fit Thin Plate Spline model
  tps_model <- Tps(coords, target_values)
  
  # Predict at test points
  prediction <- predict(tps_model, as.matrix(coordinates(test_points)))
  
  return(prediction)
}

# 2.5. GAM with Thin Plate Regression Splines Interpolation Function (Updated)
perform_gam <- function(train_data, test_points, target) {
  
  set.seed(1)
  
  # Extract training data including Longitude and Latitude
  df_train <- data.frame(
    Longitude = coordinates(train_data)[,1],
    Latitude = coordinates(train_data)[,2],
    train_data@data
  )
  
  # Fit GAM model
  gam_model <- gam(
    as.formula(paste(target, "~ s(Longitude, Latitude, bs = 'tp')")),
    data = df_train
  )
  
  # Prepare newdata with Longitude and Latitude for prediction
  newdata <- data.frame(
    Longitude = coordinates(test_points)[,1],
    Latitude = coordinates(test_points)[,2]
  )
  
  # Predict at test points
  prediction <- predict(gam_model, newdata = newdata)
  
  return(prediction)
}

# 2.6. Multilevel B-splines (MBA) Interpolation Function (Updated)
perform_mba <- function(train_data, test_points, target) {
  
  set.seed(1)
  
  # Extract coordinates and target values
  coords <- coordinates(train_data)
  target_values <- train_data@data[[target]]
  
  # Create a data frame with proper column names
  mba_data <- data.frame(x = coords[,1], y = coords[,2], z = target_values)
  
  # Remove duplicate coordinates to prevent interpolation errors
  mba_data <- mba_data[!duplicated(mba_data[, c("x", "y")]), ]
  
  # Check if there are enough unique points for interpolation
  if(nrow(mba_data) < 3){
    warning("Not enough unique data points for MBA interpolation.")
    return(rep(NA, nrow(test_points)))
  }
  
  # Determine grid size based on data extent and desired resolution
  # For example, aiming for a grid cell size of approximately 0.1 degrees
  grid_cell_size <- 0.1  # Adjust as needed
  x_range <- range(mba_data$x)
  y_range <- range(mba_data$y)
  no.X <- ceiling((x_range[2] - x_range[1]) / grid_cell_size)
  no.Y <- ceiling((y_range[2] - y_range[1]) / grid_cell_size)
  
  # Perform MBA interpolation
  mba_result <- tryCatch({
    mba.surf(
      data = mba_data,
      no.X = no.X,  # Number of grid points along X-axis
      no.Y = no.Y,  # Number of grid points along Y-axis
      extend = FALSE
    )
  }, error = function(e) {
    warning(paste("MBA interpolation failed:", e$message))
    return(NULL)
  })
  
  # If interpolation failed, return NA predictions
  if (is.null(mba_result)) {
    return(rep(NA, nrow(test_points)))
  }
  
  # Retrieve the interpolated grid
  x_grid <- mba_result$x
  y_grid <- mba_result$y
  z_grid <- mba_result$z
  
  # Initialize prediction vector
  prediction <- numeric(nrow(test_points))
  
  # Loop through each test point to assign the closest interpolated value
  for (i in 1:nrow(test_points)) {
    test_coord <- coordinates(test_points)[i, ]
    
    # Find the index of the closest grid point in X
    x_idx <- which.min(abs(x_grid - test_coord[1]))
    
    # Find the index of the closest grid point in Y
    y_idx <- which.min(abs(y_grid - test_coord[2]))
    
    # Retrieve the interpolated value from the grid
    prediction[i] <- z_grid[y_idx, x_idx]
  }
  
  return(prediction)
}

# 2.7. Nearest Neighbor Interpolation Function
perform_nn <- function(train_data, test_points, target) {
  
  set.seed(1)
  
  # Extract coordinates and target values
  train_coords <- coordinates(train_data)
  test_coords <- coordinates(test_points)
  target_values <- train_data@data[[target]]
  
  # Use RANN for fast nearest neighbor search
  nn_result <- nn2(data = train_coords, query = test_coords, k = 1)
  
  # Retrieve the nearest neighbor's value
  prediction <- target_values[nn_result$nn.idx]
  
  return(prediction)
}

# 2.8. Radial Basis Function (RBF) Interpolation Function
perform_rbf <- function(train_data, test_points, target) {
  
  set.seed(1)
  
  # Extract coordinates and target values
  coords <- coordinates(train_data)
  target_values <- train_data@data[[target]]
  
  # Fit RBF model using the 'fields' package
  rbf_model <- Tps(coords, target_values)
  
  # Predict at test points
  prediction <- predict(rbf_model, as.matrix(coordinates(test_points)))
  
  return(prediction)
}

# 2.9. Cross-Validation Function for Multiple Methods (Modified for 5-fold CV)
cross_validate_interpolation <- function(sp_data, target, methods = c("IDW", "TPS", "GAM", "MBA", "NN", "RBF"), idp = 2.0, k_folds = 5) {
  
  set.seed(1)
  
  if (!(target %in% names(sp_data))) {
    stop(paste("Target column", target, "not found in the data."))
  }
  
  n <- nrow(sp_data)
  
  # Create k-folds
  folds <- createFolds(sp_data@data[[target]], k = k_folds, list = TRUE, returnTrain = FALSE)
  
  # Initialize vectors to store predictions for each method
  pred_idw <- rep(NA, n)
  pred_tps <- rep(NA, n)
  pred_gam <- rep(NA, n)
  pred_mba <- rep(NA, n)
  pred_nn <- rep(NA, n)
  pred_rbf <- rep(NA, n)
  actual_values <- sp_data@data[[target]]
  
  # Initialize counters for failures
  idw_failures <- 0
  tps_failures <- 0
  gam_failures <- 0
  mba_failures <- 0
  nn_failures <- 0
  rbf_failures <- 0
  
  # Loop through each fold for 5-fold CV
  for (fold_idx in 1:length(folds)) {
    cat("Processing Fold", fold_idx, "of", length(folds), "...\n")
    
    test_indices <- folds[[fold_idx]]
    train_indices <- setdiff(1:n, test_indices)
    
    train <- sp_data[train_indices, ]
    test <- sp_data[test_indices, ]
    
    # Perform interpolations for each method
    # IDW
    if ("IDW" %in% methods) {
      tryCatch({
        pred_idw[test_indices] <- perform_idw(train, test, target, idp)
      }, error = function(e) {
        idw_failures <<- idw_failures + length(test_indices)
        warning(paste("IDW failed for Fold", fold_idx, ":", e$message))
      })
    }
    
    # TPS
    if ("TPS" %in% methods) {
      tryCatch({
        pred_tps[test_indices] <- perform_tps(train, test, target)
      }, error = function(e) {
        tps_failures <<- tps_failures + length(test_indices)
        warning(paste("TPS failed for Fold", fold_idx, ":", e$message))
      })
    }
    
    # GAM
    if ("GAM" %in% methods) {
      tryCatch({
        pred_gam[test_indices] <- perform_gam(train, test, target)
      }, error = function(e) {
        gam_failures <<- gam_failures + length(test_indices)
        warning(paste("GAM failed for Fold", fold_idx, ":", e$message))
      })
    }
    
    # MBA
    if ("MBA" %in% methods) {
      tryCatch({
        pred_mba[test_indices] <- perform_mba(train, test, target)
      }, error = function(e) {
        mba_failures <<- mba_failures + length(test_indices)
        warning(paste("MBA failed for Fold", fold_idx, ":", e$message))
      })
    }
    
    # Nearest Neighbor
    if ("NN" %in% methods) {
      tryCatch({
        pred_nn[test_indices] <- perform_nn(train, test, target)
      }, error = function(e) {
        nn_failures <<- nn_failures + length(test_indices)
        warning(paste("Nearest Neighbor failed for Fold", fold_idx, ":", e$message))
      })
    }
    
    # RBF
    if ("RBF" %in% methods) {
      tryCatch({
        pred_rbf[test_indices] <- perform_rbf(train, test, target)
      }, error = function(e) {
        rbf_failures <<- rbf_failures + length(test_indices)
        warning(paste("RBF failed for Fold", fold_idx, ":", e$message))
      })
    }
    
    # Progress indicator within fold
    cat("  Completed Fold", fold_idx, "\n")
  }
  
  # Compile results into a dataframe
  results <- data.frame(
    pred_idw = pred_idw,
    pred_tps = pred_tps,
    pred_gam = pred_gam,
    pred_mba = pred_mba,
    pred_nn = pred_nn,
    pred_rbf = pred_rbf,
    actual = actual_values
  )
  
  # Report failures
  if (idw_failures > 0) {
    cat("IDW failed for", idw_failures, "observations.\n")
  }
  if (tps_failures > 0) {
    cat("TPS failed for", tps_failures, "observations.\n")
  }
  if (gam_failures > 0) {
    cat("GAM failed for", gam_failures, "observations.\n")
  }
  if (mba_failures > 0) {
    cat("MBA failed for", mba_failures, "observations.\n")
  }
  if (nn_failures > 0) {
    cat("Nearest Neighbor failed for", nn_failures, "observations.\n")
  }
  if (rbf_failures > 0) {
    cat("RBF failed for", rbf_failures, "observations.\n")
  }
  
  return(results)
}

# 2.10. Calculate RMSE and R-squared
calculate_metrics <- function(results, method) {
  
  set.seed(1)
  
  actual <- results$actual
  predicted <- results[[paste0("pred_", tolower(method))]]
  
  # Remove NA values
  valid_indices <- which(!is.na(predicted) & !is.na(actual))
  actual <- actual[valid_indices]
  predicted <- predicted[valid_indices]
  
  if (length(predicted) < 2) {
    warning(paste("Not enough valid data to calculate metrics for", method, "."))
    return(list(RMSE = NA, R_squared = NA))
  }
  
  # Calculate RMSE
  rmse <- sqrt(mean((actual - predicted)^2))
  
  # Calculate R-squared
  correlation <- cor(actual, predicted)
  r_squared <- correlation^2
  
  return(list(RMSE = rmse, R_squared = r_squared))
}

# 2.11. Plot Predicted vs Actual Values
plot_predicted_vs_actual <- function(results, feature, method) {
  
  set.seed(1)
  
  actual <- results$actual
  predicted <- results[[paste0("pred_", tolower(method))]]
  
  # Remove NA values
  plot_data <- data.frame(actual = actual, predicted = predicted) %>%
    drop_na()
  
  if (nrow(plot_data) < 2) {
    warning(paste("Not enough valid data to plot for", feature, "-", method))
    return(NULL)
  }
  
  # Fit linear model
  lm_model <- lm(predicted ~ actual, data = plot_data)
  
  # Extract R-squared
  r_squared <- summary(lm_model)$r.squared
  
  # Create scatter plot with regression line
  p <- ggplot(plot_data, aes(x = actual, y = predicted)) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    labs(
      title = paste(feature, "-", method, ": Predicted vs Actual"),
      x = "Actual Value",
      y = "Predicted Value",
      subtitle = paste("R² =", round(r_squared, 3))
    ) +
    theme_minimal()
  
  return(p)
}

# ----------------------------
# 3. Main Script Execution
# ----------------------------

# 3.1. Define the path to your CSV data
data_path <- "~/Master/15hp_project/mGPS-master/mGPS-master/soil_scripts/serious_time/3_attribute_geo/round_2/interpolated_dataset_w_countries.csv" 

# 3.2. Load Data
data <- load_data(data_path)

# 3.3. Subset Data for North America
#subset_data <- data %>%
#  filter(Continent == "North.America")  # Adjust if continent names differ

subset_data <- data

# 3.4. Transform to Spatial Object
sp_data <- transform_to_spatial_object(subset_data)

# 3.5. Define Features to Interpolate (OTU Columns)
features <- grep("^OTU", names(subset_data), value = TRUE)

# 3.6. Initialize a dataframe to store evaluation metrics for each method and feature
evaluation_metrics <- data.frame(
  Feature = character(),
  Method = character(),
  RMSE = numeric(),
  R_squared = numeric(),
  stringsAsFactors = FALSE
)

# 3.7. Initialize a list to store plots
plot_list <- list()

# 3.8. Perform 5-fold CV for Each Feature and Each Method
for (feature in features) {
  cat("Performing 5-fold CV for", feature, "...\n")
  
  # Perform 5-fold CV for IDW, TPS, GAM, MBA, Nearest Neighbor, and RBF
  cv_results <- cross_validate_interpolation(
    sp_data, 
    target = feature, 
    methods = c("IDW", "TPS", "GAM", "MBA", "NN", "RBF"), 
    idp = 2.0,
    k_folds = 5  # Specify number of folds
  )
  
  # Calculate metrics for IDW
  metrics_idw <- calculate_metrics(cv_results, "IDW")
  
  # Calculate metrics for TPS
  metrics_tps <- calculate_metrics(cv_results, "TPS")
  
  # Calculate metrics for GAM
  metrics_gam <- calculate_metrics(cv_results, "GAM")
  
  # Calculate metrics for MBA
  metrics_mba <- calculate_metrics(cv_results, "MBA")
  
  # Calculate metrics for Nearest Neighbor
  metrics_nn <- calculate_metrics(cv_results, "NN")
  
  # Calculate metrics for RBF
  metrics_rbf <- calculate_metrics(cv_results, "RBF")
  
  # Store evaluation metrics
  evaluation_metrics <- evaluation_metrics %>%
    add_row(
      Feature = feature,
      Method = "IDW",
      RMSE = metrics_idw$RMSE,
      R_squared = metrics_idw$R_squared
    ) %>%
    add_row(
      Feature = feature,
      Method = "TPS",
      RMSE = metrics_tps$RMSE,
      R_squared = metrics_tps$R_squared
    ) %>%
    add_row(
      Feature = feature,
      Method = "GAM",
      RMSE = metrics_gam$RMSE,
      R_squared = metrics_gam$R_squared
    ) %>%
    add_row(
      Feature = feature,
      Method = "MBA",
      RMSE = metrics_mba$RMSE,
      R_squared = metrics_mba$R_squared
    ) %>%
    add_row(
      Feature = feature,
      Method = "Nearest_NN",
      RMSE = metrics_nn$RMSE,
      R_squared = metrics_nn$R_squared
    ) %>%
    add_row(
      Feature = feature,
      Method = "RBF",
      RMSE = metrics_rbf$RMSE,
      R_squared = metrics_rbf$R_squared
    )
  
  # Generate and store plots for each method
  plot_idw <- plot_predicted_vs_actual(cv_results, feature, "IDW")
  plot_tps <- plot_predicted_vs_actual(cv_results, feature, "TPS")
  plot_gam <- plot_predicted_vs_actual(cv_results, feature, "GAM")
  plot_mba <- plot_predicted_vs_actual(cv_results, feature, "MBA")
  plot_nn <- plot_predicted_vs_actual(cv_results, feature, "NN")
  plot_rbf <- plot_predicted_vs_actual(cv_results, feature, "RBF")
  
  # Combine the plots in a grid layout
  combined_plot <- grid.arrange(
    plot_idw, plot_tps, plot_gam, plot_mba, plot_nn, plot_rbf, 
    ncol = 3, 
    top = paste("Predicted vs Actual for", feature)
  )
  
  plot_list[[feature]] <- combined_plot
}

# 3.9. Display Evaluation Metrics
print(evaluation_metrics)

# 3.10. Visualize Performance Comparison Across Methods

# a. RMSE Comparison
rmse_plot <- ggplot(evaluation_metrics, aes(x = Feature, y = RMSE, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(
    title = "RMSE Comparison Across Interpolation Methods", 
    x = "Feature", 
    y = "RMSE"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# b. R-squared Comparison
r2_plot <- ggplot(evaluation_metrics, aes(x = Feature, y = R_squared, color = Method)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_line(aes(group = Method), position = position_dodge(width = 0.5)) +
  theme_minimal() +
  labs(
    title = "R-squared Comparison Across Interpolation Methods", 
    x = "Feature", 
    y = "R-squared"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# c. Arrange RMSE and R-squared plots side by side
grid.arrange(rmse_plot, r2_plot, ncol = 2)

# Step 1: Remove rows with NA in R_squared
clean_evaluation_metrics <- evaluation_metrics %>%
  filter(!is.na(R_squared))

# Step 2: For each Feature, select the Method with the highest R_squared
high_r2 <- clean_evaluation_metrics %>%
  group_by(Feature) %>%
  slice_max(order_by = R_squared, with_ties = TRUE) %>%  # with_ties = FALSE ensures only one row per group
  ungroup()

ggplot(high_r2, aes(x = Feature, y = R_squared, fill = Method)) +
  geom_bar(stat = "identity", width = 0.6) +
  theme_minimal() +
  labs(
    title = "Highest R² Value per Feature by Method",
    x = "Feature",
    y = expression(R^2),
    fill = "Method"
  ) +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )

# ----------------------------
# 4. Grid Interpolation with 1000 km Constraint
# ----------------------------

# 4.1. Convert SpatialPointsDataFrame to sf object for easier handling
sp_sf <- st_as_sf(sp_data)

# 4.2. Create a grid over the study area
# Define grid resolution (e.g., 0.5 degrees)
grid_resolution <- 2  # Adjust as needed

# Create a grid covering the extent of the data
grid <- st_make_grid(sp_sf, cellsize = grid_resolution, square = TRUE)

# Convert grid to sf object
grid_sf <- st_sf(geometry = grid)

# 4.3. Compute which grid points are within 1000 km of any data point
# Since 'st_is_within_distance' uses meters, convert 1000 km to meters
distance_threshold <- 50000  # 20 km in meters

# Determine which grid points are within 1000 km of any data point
within_1000km <- st_is_within_distance(grid_sf, sp_sf, dist = distance_threshold)

# Convert the list to a logical vector indicating if any data point is within 1000 km
grid_sf$within_1000km <- lengths(within_1000km) > 0

# 4.4. Filter grid points to keep only those within 1000 km
filtered_grid_sf <- grid_sf %>%
  filter(within_1000km) %>%
  select(-within_1000km)

# 4.5. Transform filtered grid back to SpatialPoints for compatibility
filtered_grid_sp <- as(filtered_grid_sf, "Spatial")

# 4.6. Perform Interpolation on Filtered Grid Points
# Define the interpolation methods to use on the grid
grid_methods <- c("IDW", "TPS", "GAM", "MBA", "NN", "RBF")

# Initialize a list to store interpolation results for each feature and method
grid_interpolation_results <- list()

for (feature in features) {
  cat("Performing Grid Interpolation for", feature, "...\n")
  
  # Initialize a list to store predictions for each method
  predictions <- list()
  
  # Prepare training data (all original data points)
  train <- sp_data
  
  # IDW
  if ("IDW" %in% grid_methods) {
    cat("  Performing IDW...\n")
    predictions$IDW <- perform_idw(train, filtered_grid_sp, feature, idp = 2.0)
  }
  
  # TPS
  if ("TPS" %in% grid_methods) {
    cat("  Performing TPS...\n")
    predictions$TPS <- perform_tps(train, filtered_grid_sp, feature)
  }
  
  # GAM
  if ("GAM" %in% grid_methods) {
    cat("  Performing GAM...\n")
    predictions$GAM <- perform_gam(train, filtered_grid_sp, feature)
  }
  
  # MBA
  if ("MBA" %in% grid_methods) {
    cat("  Performing MBA...\n")
    predictions$MBA <- perform_mba(train, filtered_grid_sp, feature)
  }
  
  # Nearest Neighbor
  if ("NN" %in% grid_methods) {
    cat("  Performing Nearest Neighbor...\n")
    predictions$NN <- perform_nn(train, filtered_grid_sp, feature)
  }
  
  # RBF
  if ("RBF" %in% grid_methods) {
    cat("  Performing RBF...\n")
    predictions$RBF <- perform_rbf(train, filtered_grid_sp, feature)
  }
  
  # Combine predictions into a dataframe
  grid_results <- data.frame(
    Longitude = coordinates(filtered_grid_sp)[,1],
    Latitude = coordinates(filtered_grid_sp)[,2],
    predictions
  )
  
  # Store in the main list
  grid_interpolation_results[[feature]] <- grid_results
  
  cat("  Completed Grid Interpolation for", feature, ".\n")
}



## Extract and save the interpolated data from the methods that performed the best

clean_evaluation_metrics <- evaluation_metrics %>%
  filter(!is.na(R_squared))

# Step 2: For each Feature, select the Method with the highest R_squared
high_r2 <- clean_evaluation_metrics %>%
  group_by(Feature) %>%
  slice_max(order_by = R_squared, n = 1, with_ties = FALSE) %>%  # with_ties = FALSE ensures only one row per group
  ungroup()

toSave <- data.frame(Longitude = grid_results$Longitude,
                     Latitude = grid_results$Latitude
)


for (feature in features) {
  # Check if the feature exists in grid_interpolation_results
  if (!feature %in% names(grid_interpolation_results)) {
    warning(paste("Feature", feature, "not found in grid_interpolation_results"))
    next  # Skip to the next iteration
  }
  
  # Extract f
  f <- grid_interpolation_results[[feature]]
  
  # Check if the feature exists in high_r2
  if (!any(high_r2$Feature == feature)) {
    warning(paste("Feature", feature, "not found in high_r2"))
    next  # Skip to the next iteration
  }
  
  # Subset high_r2 for the current feature
  m <- high_r2[high_r2$Feature == feature, ]
  
  if (m$Method == "Nearest_NN") {
    m$Method = "NN"
  }
  
  # Check if 'Method' column exists
  if (!"Method" %in% colnames(m)) {
    warning(paste("'Method' column not found for feature", feature, "in high_r2"))
    next
  }
  
  # Extract best_method
  best_method <- m$Method
  
  # Assign to toSave
  toSave[[feature]] <- f[, best_method]  # Using [[ ]] is safer for non-syntactic names
}

# Ensure that 'toSave' has all the columns present in 'data', filling missing ones with NA
missing_columns <- setdiff(names(data), names(toSave))

# Add missing columns to 'toSave' with NA values
for (col in missing_columns) {
  toSave[[col]] <- NA
}

# Reorder columns in 'toSave' to match 'data'
toSave <- toSave[, names(data)]

# Now, combine 'data' and 'toSave'
tot_data <- bind_rows(data, toSave)

# Optional: If you want to remove specific columns like 'Country' and 'Continent' from 'tot_data'
tot_data <- tot_data %>% select(-Country, -Continent)

# Create leaflet map to visualize the data points
leaflet(tot_data) %>%
  addTiles() %>%
  addCircleMarkers(lng = tot_data$Longitude, lat = tot_data$Latitude, radius = 1)

# Save the combined data to a CSV file
write_csv(tot_data, "~/Master/15hp_project/mGPS-master/mGPS-master/soil_scripts/serious_time/2_inteprolation/iterations/2x_interpol/interpolated_dataset.csv")
