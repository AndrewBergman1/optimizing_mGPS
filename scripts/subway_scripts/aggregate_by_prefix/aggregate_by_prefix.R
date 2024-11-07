#################################################

# Author: Andrew Bergman
# Course: BINP37, Optimizing mGPS 

################################################

# Description:
#   Aggregate subway data by the first string when splittingh taxa name with "." 


# Extract all prefixes of all columns passed
extract_all_prefixes <- function(names) {
  # Split column names by "."
  split_names <- lapply(strsplit(names, "\\."), function(x) x[1])  # Get the first part of the split
  # Find all unique prefixes
  all_prefixes <- unique(unlist(split_names))
  return(all_prefixes)
}

# Aggregate the data by the extracted prefixes
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

# Load data
meta_abundance <- read.csv("C:/Users/andre/OneDrive/Documents/Master/15hp_project/mGPS-master/mGPS-master/Data/Data/Metasub/meta_taxa_read.csv")

# Preprocess data
meta_abundance <- data_preprocess_f(meta_abundance, "city", c("continent", "city", "latitude", "longitude"), 8)

# Get all prefixes
all_prefixes <- extract_all_prefixes(names(meta_abundance[44:ncol(meta_abundance)]))
aggregated_meta_abundance <- aggregate_data_by_prefix(meta_abundance[, 44:ncol(meta_abundance)], all_prefixes)


# Select the first 43 columns
initial_columns <- meta_abundance[, 1:43]
# Combine with the aggregated data
aggregated_meta_abundance <- cbind(initial_columns, aggregated_meta_abundance)

# Save data
write.csv(aggregated_meta_abundance, "C:/Users/andre/OneDrive/Documents/Master/15hp_project/mGPS-master/mGPS-master/Data/Data/Metasub/Families_metasub.csv")

