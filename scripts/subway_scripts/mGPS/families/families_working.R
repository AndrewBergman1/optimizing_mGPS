library(caret)
library(geosphere)
library(dplyr)
library(xgboost)
library(e1071)
library(doParallel)
library(leaflet)

data_normalise <- function(df) {
  return(df/rowSums(df))
  
}

data_preprocess_f <- function(train_f,target_in,hierarchy,remove_small){
  
  #Remove control samples
  metasub_data <- droplevels(train_f)
  
  if(length(hierarchy) == 4){
    metasub_data[,hierarchy[1]] <- factor(metasub_data[,hierarchy[1]])
    metasub_data[,hierarchy[2]] <- factor(metasub_data[,hierarchy[2]])
  }else{
    metasub_data[,hierarchy[1]] <- factor(metasub_data[,hierarchy[1]])
  }
  #remove sparse samples locations and dubiously labelled samples. 
  
  remove_small <- as.numeric(remove_small)
  
  if (remove_small > 0){
    
    small_cities <-  names(which(summary(metasub_data[,target_in]) < remove_small))
    
    remove_samples <- which(metasub_data[,target_in] %in%  c("antarctica", small_cities))
    
    if (length(remove_samples) != 0 ){
      metasub_data <- droplevels(metasub_data[-c(remove_samples), ])
    }
  }
  
  #Correction of identified misslabelling of data 
  
  for (i in 1:length(hierarchy)){
    empty <- which(metasub_data[,hierarchy[i]] == "" | is.na(metasub_data[,hierarchy[i]]))
    
    if (length(empty) != 0  ){
      metasub_data <- metasub_data[-c(empty),]
    }
  }
  
  metasub_data[,hierarchy[1]] <- droplevels(metasub_data[,hierarchy[1]])
  
  if(length(hierarchy) == 4){
    metasub_data[,hierarchy[2]] <- droplevels(metasub_data[,hierarchy[2]])
  }
  
  print("metadata has been pre-processed ...")
  return(metasub_data)
}

##Feature selection algorithm
species_select <- function(x,
                           y,
                           remove_correlated = TRUE,
                           subsets = NULL,
                           cores = 1) {
  # Start parallel processing
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # Handle correlations
  if (remove_correlated) {
    # Safely compute correlation matrix
    cor_matrix <- cor(x, use = "complete.obs")
    if(any(is.na(cor_matrix))) {
      warning("NA values found in correlation matrix. Check for constant variables.")
      cor_matrix[is.na(cor_matrix)] <- 0  # Replace NAs with 0 to allow `findCorrelation` to proceed
    }
    
    correlated <- findCorrelation(cor_matrix, cutoff = 0.98, verbose = TRUE, names = FALSE)
    
    print("Correlated features indices:")
    print(correlated)
    
    if (length(correlated) > 0) {
      x <- x[, -correlated, drop = FALSE]  # Ensure matrix structure is maintained
      print(paste("Correlated features removed: ", length(correlated)))
    } else {
      print("No correlated features to remove.")
    }
  }
  
  # Define subsets if not provided
  if(is.null(subsets)){
    len <- ncol(x)
    subsets <- c(floor(len/2), floor(len / 4), floor(len / 8), floor(len / 16), floor(len / 32), floor(len / 64))
  }
  
  # Define RFE control
  rfe_ctrl <- rfeControl(
    functions = rfFuncs,
    method = "cv",
    number =  5,
    verbose = TRUE,
    allowParallel = TRUE
  )
  
  set.seed(123)
  featureElimination <- rfe(x = x, y = y, sizes = subsets, rfeControl = rfe_ctrl, tuneLength = 2)
  
  var_importance <- varImp(featureElimination$fit, scale = FALSE)
  
  
  
  # Stop parallel processing
  stopCluster(cl)
  registerDoParallel(1)  # Reset to default single-core processing
  
  # Return both the featureElimination object and the names of selected features
  selected_features <- colnames(x[, featureElimination$optVariables, drop = FALSE])
  return(list(results = featureElimination, selectedFeatures = selected_features, importance = var_importance))
}

##Main mGPS algorithm 
mGPS <-
  function(training = training_data,
           testing = testing_data,
           classTarget = "city",
           hierarchy = c('continent','city','latitude','longitude'),
           variables,
           nthread = 1,
           coast = NULL) {
    #Check for training set
    
    if(is.null(training)){
      return(message("No training set given"))
    } else{
      
      
      training <- droplevels(training)
      
      #Train mGPS with 5-fold cross validation of training set for hyperparameter tuning. 
      message("Training mGPS...")
      
      set.seed(1234)
      folds <- createFolds(training[,classTarget], k = 5, returnTrain = T)
      
      
      
      trControlClass <-  trainControl(
        method = "cv",
        number = 5,  
        verboseIter = FALSE,
        returnData = FALSE,
        search = "grid",
        savePredictions = "final",
        classProbs = T,
        allowParallel = T,
        index = folds )
      
      
      
      trControl <-  trainControl(
        method = "cv",
        number = 5,  
        verboseIter = FALSE,
        returnData = FALSE,
        search = "grid",
        savePredictions = "final",
        allowParallel = T,
        index = folds)
      
      
      
      tune_grid <- expand.grid(
        nrounds = c(300,600),
        eta = c( 0.05, 0.1),
        max_depth = c(3,6,9),
        gamma = 0,
        colsample_bytree = c(0.6,0.8),
        min_child_weight = c(1),
        subsample = (0.7)
      )
      
      if(length(hierarchy) == 4){
        
        Xgb_region <- train(x = training[,variables],y = training[,hierarchy[1]],
                            method = "xgbTree",
                            trControl = trControlClass,
                            tuneGrid = tune_grid,
                            nthread = nthread)
        
        l1_train <- data.frame(training[,c(variables)],Xgb_region[["pred"]][order(Xgb_region$pred$rowIndex),levels(training[,hierarchy[1]]) ])
        
      }else{
        
        l1_train <- training[,c(variables)]
        
      }
      Xgb_class <- train(x = l1_train,y = training[,classTarget],
                         method = "xgbTree",
                         trControl = trControlClass,
                         tuneGrid = tune_grid,
                         nthread = nthread)
      
      l2_train <- data.frame(l1_train,Xgb_class[["pred"]][order(Xgb_class$pred$rowIndex),levels(training[,classTarget]) ])
      
      
      Xgb_latitude <- train(x = l2_train ,y = training[,"latitude"],
                            method = "xgbTree",
                            trControl = trControl,
                            tuneGrid = tune_grid,
                            nthread = nthread)
      
      l3_train <- data.frame(l2_train, "latPred" = Xgb_latitude[["pred"]][order(Xgb_latitude$pred$rowIndex),"pred" ])
      
      Xgb_longitude <- train(x = l3_train ,y = training[,"longitude"],
                             method = "xgbTree",
                             trControl = trControl,
                             tuneGrid = tune_grid,
                             nthread = nthread)
      
    }
    #check for test set, return trained model if no test set. 
    if(is.null(testing)){
      
      model <- function(test,variables){
        regProbs <- predict(Xgb_region, newdata = test[,variables],type ="prob")
        
        l1_test <- data.frame(test[,variables], regProbs)
        
        classPred <- predict(Xgb_class, newdata = l1_test)
        classProbs <- predict(Xgb_class, newdata = l1_test,type ="prob")
        
        l2_test <-  data.frame(l1_test, classProbs) 
        latPred <- predict(Xgb_latitude, newdata = l2_test)
        
        l3_test <- data.frame(l2_test, latPred)
        longPred <- predict(Xgb_longitude, newdata = l3_test)
        return(list(classPred, latPred, longPred))
        
      }
      message("No test set...returning trained mGPS model function")
      return(list(Xgb_region,Xgb_class,Xgb_latitude,Xgb_longitude,"model" = model))
    }else{
      message("Generating predictions")
      #generate mGPS predictions for test set
      if(length(hierarchy) == 4){
        regProbs <- predict(Xgb_region, newdata = testing[,variables],type ="prob")
        
        l1_test <- data.frame(testing[,variables], regProbs)
      }else{
        l1_test <- testing[,variables]
      }
      classPred <- predict(Xgb_class, newdata = l1_test)
      classProbs <- predict(Xgb_class, newdata = l1_test,type ="prob")
      
      l2_test <-  data.frame(l1_test, classProbs) 
      latPred <- predict(Xgb_latitude, newdata = l2_test)
      
      l3_test <- data.frame(l2_test, latPred)
      longPred <- predict(Xgb_longitude, newdata = l3_test)
      
      #uuids <- testData$uuid
      
      #adjust out of bounds predictions
      longPred[longPred > 180] <- 180
      longPred[longPred < -180] <- -180
      latPred[latPred > 90] <- 90
      latPred[latPred < -90] <- -90
      #Pull to nearest coastline if provided
      find_coast <- function(long, lat) {
        distances_from_coastline <-
          sp::spDistsN1(coast, c(long, lat), longlat = TRUE)
        
        closest_point <-  which.min(distances_from_coastline)
        new_coords <- coast[closest_point,]
        
        return(new_coords)
        
      }
      if (!is.null(coast)) {
        toAdjust <-
          which(is.na(maps::map.where(database = "world", longPred, latPred)))
        
        adjusted <-
          mapply(find_coast, long = longPred[toAdjust], lat = latPred[toAdjust])
        
        
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
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  distance <- R * c
  
  return(distance)  # Distance in kilometers
}


# Function to calculate distance between corresponding pairs of coordinates
calculate_distances <- function(latitudes1, longitudes1, latitudes2, longitudes2) {
  distances <- haversine_distance(latitudes1, longitudes1, latitudes2, longitudes2)
  return(distances)
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
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  distance <- R * c
  
  return(distance)  # Distance in kilometers
}

# Function to calculate distance between corresponding pairs of coordinates
calculate_distances <- function(latitudes1, longitudes1, latitudes2, longitudes2) {
  distances <- haversine_distance(latitudes1, longitudes1, latitudes2, longitudes2)
  return(distances)
}

pull_land <- function(land_preds,hierarchy){
  
  coastlines <- cbind("x"  = maps::SpatialLines2map(rworldmap::coastsCoarse)$x ,"y" =maps::SpatialLines2map(rworldmap::coastsCoarse)$y)
  coastlines <- coastlines[complete.cases(coastlines),]
  coastlines <- coastlines[coastlines[,1] < 180 ,]
  
  find_coast <- function(long, lat) {
    
    distances_from_coastline <-
      sp::spDistsN1(coastlines, c(long, lat), longlat = TRUE)
    
    closest_point <-  which.min(distances_from_coastline)
    new_coords <- coastlines[closest_point,]
    
    return(new_coords)
  }
  
  toAdjust <-
    which(is.na(maps::map.where(database = "world", land_preds$longPred, land_preds$latPred)))
  
  adjusted <-
    mapply(find_coast, long = land_preds$longPred[toAdjust], lat = land_preds$latPred[toAdjust])
  
  if (length(adjusted) != 0){
    land_preds$longPred[toAdjust] <- adjusted[1,]
    land_preds$latPred[toAdjust] <- adjusted[2,]
  }
  
  print(hierarchy)
  
  if (class(hierarchy)=='character' && hierarchy[3] %in% colnames(land_preds)){
    for (i in 1:nrow(land_preds)){
      
      if(length(hierarchy) == 3){
        
        land_preds[i,"Distance_from_origin"] <- geosphere::distm(c(land_preds[i,"longPred"],land_preds[i,"latPred"]), c(land_preds[i,hierarchy[3]],land_preds[i,hierarchy[2]]), fun = geosphere::distHaversine)/1000
        
      }else  if(length(hierarchy) == 4){
        land_preds[i,"Distance_from_origin"] <- geosphere::distm(c(land_preds[i,"longPred"],land_preds[i,"latPred"]), c(land_preds[i,hierarchy[4]],land_preds[i,hierarchy[3]]), fun = geosphere::distHaversine)/1000
      }
    }
  }
  
  return(land_preds)
}


extract_all_prefixes <- function(names) {
  # Split column names by "."
  split_names <- lapply(strsplit(names, "\\."), function(x) x[1])  # Get the first part of the split
  # Find all unique prefixes
  all_prefixes <- unique(unlist(split_names))
  return(all_prefixes)
}

aggregate_data_by_prefix <- function(data, all_prefixes) {
  # Initialize an empty data frame to store aggregated results
  aggregated_data <- data.frame(row.names = rownames(data))
  
  for (prefix in all_prefixes) {
    # Find columns that start with the current prefix
    cols_with_prefix <- grep(paste0("^", prefix, "\\."), names(data), value = TRUE)
    
    # Ensure we do not drop dimensions when subsetting a single column
    if (length(cols_with_prefix) > 0) {
      subset_data <- data[, cols_with_prefix, drop = FALSE]
      
      # Sum these columns and create a new column in the aggregated data
      aggregated_data[[prefix]] <- rowSums(subset_data, na.rm = TRUE)
    }
  }
  
  return(aggregated_data)
}

meta_abundance <- read.csv("C:/Users/andre/OneDrive/Documents/Master/15hp_project/mGPS-master/mGPS-master/Data/Data/Metasub/meta_taxa_read.csv")
meta_abundance <- data_preprocess_f(meta_abundance, "city", c("continent", "city", "latitude", "longitude"), 8)


########## EXPERIMENTING WITH WIDER MICROORGANISM CLASSIFICATION!!!!

all_prefixes <- extract_all_prefixes(names(meta_abundance[44:ncol(meta_abundance)]))
aggregated_meta_abundance <- aggregate_data_by_prefix(meta_abundance[, 44:ncol(meta_abundance)], all_prefixes)

# Select the first 43 columns
initial_columns <- meta_abundance[, 1:43]
# Combine with the aggregated data
meta_abundance <- cbind(initial_columns, aggregated_meta_abundance)
meta_abundance[44:ncol(meta_abundance)] <- data_normalise(meta_abundance[44:ncol(meta_abundance)])

trainIndex <- createDataPartition(meta_abundance$city, p = 0.8, list = FALSE)
trainData <- meta_abundance[trainIndex, ]
testData <- meta_abundance[-trainIndex, ]


trainData[44:ncol(trainData)] <- data_normalise(trainData[44:ncol(trainData)])
trainData[, 6] <- as.factor(trainData[, 6])
trainData[44:ncol(trainData)] <- na.omit(trainData[44:ncol(trainData)])

testData[44:ncol(testData)] <- data_normalise(testData[44:ncol(testData)])
testData[, 6] <- as.factor(testData[, 6])
testData[44:ncol(testData)] <- na.omit(testData[44:ncol(testData)])


features <- species_select(trainData[44:ncol(trainData)], trainData$city, remove_correlated = T, c(50,100,150,200,250), cores = 4)

top_features <- unlist(features[2])

preds <- mGPS(trainData, testData, classTarget = "city", hierarchy = c("continent", "city", "latitude", "longitude"), variables = top_features, nthread=8)

preds <- pull_land(preds, c("continent", "city", "latitude", "longitude"))


map <- leaflet() %>%
  addTiles() %>%  
  addCircleMarkers(lng = preds[[3]], lat = preds[[2]], radius = 2, popup = "Predicted Coordinates") %>% 
  addCircleMarkers(lng = testData$longitude, lat = testData$latitude, radius = 0.5, color = colors, popup = "Test Data Coordinates") 

map

num_matches <- sum(preds[[1]] == testData$city)
distances <- calculate_distances(testData$latitude, testData$longitude, preds[[2]], preds[[3]])

# Calculate statistics
mean_distance <- mean(distances)
median_distance <- median(distances)
min_distance <- min(distances)
max_distance <- max(distances)
num_distances_less_than_100 <- length(which(distances < 100))
num_distances_less_than_250 <- length(which(distances < 250))
num_distances_less_than_500 <- length(which(distances < 500))
num_distances_less_than_1000 <- length(which(distances < 1000))
num_distances_less_than_1500 <- length(which(distances < 1500))
total_distances <- length(distances)

# Print statistics
cat("Statistic\t\tValue\n")
cat("Mean distance\t\t", mean_distance, "\n")
cat("Median distance\t", median_distance, "\n")
cat("Minimum distance\t", min_distance, "\n")
cat("Maximum distance\t", max_distance, "\n")
cat("Distances < 100\t\t", num_distances_less_than_100, "\n")
cat("Distances < 250\t\t", num_distances_less_than_250, "\n")
cat("Distances < 500\t\t", num_distances_less_than_500, "\n")
cat("Distances < 1000\t", num_distances_less_than_1000, "\n")
cat("Distances < 1500\t", num_distances_less_than_1500, "\n")
cat("Total distances\t\t", total_distances, "\n")
cat("classif. acc.\t\t", num_matches/length(preds[[1]]), "\n")

