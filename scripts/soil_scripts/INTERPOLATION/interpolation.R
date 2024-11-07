#################################################

# Author: Andrew Bergman
# Course: BINP37, Optimizing mGPS 

################################################

# Description:
#   Determine the best interpolation method for all features in soil dataset.

set.seed(1)

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
perform_idw <- function(train_data, test_point, target, idp = 2.0) {
  
  set.seed(1)
  
  idw_result <- idw(
    formula = as.formula(paste(target, "~ 1")),
    locations = train_data,
    newdata = test_point,
    idp = idp
  )
  return(idw_result$var1.pred[1])
}

# 2.4. Thin Plate Spline (TPS) Interpolation Function
perform_tps <- function(train_data, test_point, target) {
  
  set.seed(1)
  
  # Extract coordinates and target values
  coords <- coordinates(train_data)
  target_values <- train_data@data[[target]]
  
  # Fit Thin Plate Spline model
  tps_model <- Tps(coords, target_values)
  
  # Predict at test point
  prediction <- predict(tps_model, as.matrix(coordinates(test_point)))
  
  return(prediction)
}

# 2.5. GAM with Thin Plate Regression Splines Interpolation Function
perform_gam <- function(train_data, test_point, target) {
  
  set.seed(1)
  
  # Extract data
  df_train <- as.data.frame(train_data)
  
  # Fit GAM model
  gam_model <- gam(
    as.formula(paste(target, "~ s(Longitude, Latitude, bs = 'tp')")),
    data = df_train
  )
  
  # Predict at test point
  prediction <- predict(gam_model, newdata = as.data.frame(test_point))
  
  return(prediction)
}

# 2.6. Multilevel B-splines (MBA) Interpolation Function
perform_mba <- function(train_data, test_point, target) {
  
  set.seed(1)
  
  # Extract coordinates and target values
  coords <- coordinates(train_data)
  target_values <- train_data@data[[target]]
  
  # Perform MBA interpolation
  mba_result <- mba.surf(
    data = cbind(coords, target_values),
    no.X = 1,
    no.Y = 1,
    extend = FALSE
  )
  
  # Find the closest grid point in the MBA result to the test point
  x_grid <- mba_result$x
  y_grid <- mba_result$y
  z_grid <- mba_result$z
  
  # Calculate the closest grid indices
  x_idx <- which.min(abs(x_grid - coordinates(test_point)[,1]))
  y_idx <- which.min(abs(y_grid - coordinates(test_point)[,2]))
  
  # Retrieve the interpolated value
  prediction <- z_grid[y_idx, x_idx]
  
  return(prediction)
}

# 2.7. Nearest Neighbor Interpolation Function
perform_nn <- function(train_data, test_point, target) {
  
  set.seed(1)
  
  # Extract coordinates and target values
  train_coords <- coordinates(train_data)
  test_coord <- coordinates(test_point)
  target_values <- train_data@data[[target]]
  
  # Use RANN for fast nearest neighbor search
  nn_result <- nn2(data = train_coords, query = test_coord, k = 1)
  
  # Retrieve the nearest neighbor's value
  prediction <- target_values[nn_result$nn.idx]
  
  return(prediction)
}

# 2.8. Radial Basis Function (RBF) Interpolation Function
perform_rbf <- function(train_data, test_point, target) {
  
  set.seed(1)
  
  # Extract coordinates and target values
  coords <- coordinates(train_data)
  target_values <- train_data@data[[target]]
  
  # Fit RBF model using the 'fields' package
  rbf_model <- Tps(coords, target_values)
  
  # Predict at test point
  prediction <- predict(rbf_model, as.matrix(coordinates(test_point)))
  
  return(prediction)
}

# 2.9. Cross-Validation Function for Multiple Methods
cross_validate_interpolation <- function(sp_data, target, methods = c("IDW", "TPS", "GAM", "MBA", "NN", "RBF"), idp = 2.0) {
  
  set.seed(1)
  
  if (!(target %in% names(sp_data))) {
    stop(paste("Target column", target, "not found in the data."))
  }
  
  n <- nrow(sp_data)
  
  # Initialize vectors to store predictions for each method
  pred_idw <- numeric(n)
  pred_tps <- numeric(n)
  pred_gam <- numeric(n)
  pred_mba <- numeric(n)
  pred_nn <- numeric(n)
  pred_rbf <- numeric(n)
  actual_values <- numeric(n)
  
  # Initialize counters for failures
  idw_failures <- 0
  tps_failures <- 0
  gam_failures <- 0
  mba_failures <- 0
  nn_failures <- 0
  rbf_failures <- 0
  
  # Loop through each observation for LOOCV
  for (i in 1:n) {
    # Split data into training and test sets
    train <- sp_data[-i, ]
    test <- sp_data[i, ]
    
    # Perform IDW interpolation if selected
    if ("IDW" %in% methods) {
      pred_idw[i] <- tryCatch(
        perform_idw(train, test, target, idp),
        error = function(e) {
          idw_failures <<- idw_failures + 1
          warning(paste("IDW failed for observation", i, ":", e$message))
          return(NA)
        }
      )
    }
    
    # Perform Thin Plate Spline interpolation if selected
    if ("TPS" %in% methods) {
      pred_tps[i] <- tryCatch(
        perform_tps(train, test, target),
        error = function(e) {
          tps_failures <<- tps_failures + 1
          warning(paste("TPS failed for observation", i, ":", e$message))
          return(NA)
        }
      )
    }
    
    # Perform GAM interpolation if selected
    if ("GAM" %in% methods) {
      pred_gam[i] <- tryCatch(
        perform_gam(train, test, target),
        error = function(e) {
          gam_failures <<- gam_failures + 1
          warning(paste("GAM failed for observation", i, ":", e$message))
          return(NA)
        }
      )
    }
    
    # Perform MBA interpolation if selected
    if ("MBA" %in% methods) {
      pred_mba[i] <- tryCatch(
        perform_mba(train, test, target),
        error = function(e) {
          mba_failures <<- mba_failures + 1
          warning(paste("MBA failed for observation", i, ":", e$message))
          return(NA)
        }
      )
    }
    
    # Perform Nearest Neighbor interpolation if selected
    if ("NN" %in% methods) {
      pred_nn[i] <- tryCatch(
        perform_nn(train, test, target),
        error = function(e) {
          nn_failures <<- nn_failures + 1
          warning(paste("Nearest Neighbor failed for observation", i, ":", e$message))
          return(NA)
        }
      )
    }
    
    # Perform Radial Basis Function interpolation if selected
    if ("RBF" %in% methods) {
      pred_rbf[i] <- tryCatch(
        perform_rbf(train, test, target),
        error = function(e) {
          rbf_failures <<- rbf_failures + 1
          warning(paste("RBF failed for observation", i, ":", e$message))
          return(NA)
        }
      )
    }
    
    # Store actual value
    actual_values[i] <- as.numeric(test@data[[target]])
    
    # Progress indicator (optional)
    if (i %% 10 == 0 || i == n) {
      cat("Processed", i, "of", n, "observations.\n")
    }
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
data_path <- "~/Master/15hp_project/mGPS-master/mGPS-master/soil_scripts/serious_time/SMOTE/oversampling_to_42_per_country/SMOTEd_soil_samples.csv"

# 3.2. Load Data
data <- read_csv(data_path)

# 3.3. Subset Data for North America
subset_data <- data %>%
  filter(Continent == "North.America")  # Adjust if continent names differ

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

# 3.8. Perform LOOCV for Each Feature and Each Method
for (feature in features) {
  cat("Performing LOOCV for", feature, "...\n")
  
  # Perform LOOCV for IDW, TPS, GAM, MBA, Nearest Neighbor, and RBF
  cv_results <- cross_validate_interpolation(
    sp_data, 
    target = feature, 
    methods = c("IDW", "TPS", "GAM", "MBA", "NN", "RBF"), 
    idp = 2.0
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
r2_plot <- ggplot(evaluation_metrics, aes(x = Feature, y = R_squared, col = Method)) +
  geom_point(stat = "identity", position = position_dodge()) +
  geom_line(group = "1") +
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
  slice_max(order_by = R_squared, n = 1, with_ties = FALSE) %>%  # with_ties = FALSE ensures only one row per group
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
