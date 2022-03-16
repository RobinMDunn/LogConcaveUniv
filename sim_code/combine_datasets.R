# Merge all .csv files with a given file name prefix into a single .csv file.
# Change start_name (line 8) and file name (line 52) for different data.

# Read in library
library(LogConcaveUniv)

# Set common file name start for all files to merge
start_name <- "^fig01_fully_NP_randproj_"

# Get all file names with given start_name
filenames <- list.files(path = "sim_data", pattern = start_name, full.names = T)

# Read in first file.
file_1 <- data.table::fread(file = filenames[1])

# Check if dataset includes test_stat or avg_ts columns.
# If yes, read in datasets individually, setting these columns to numeric.
# (Some large values may have been stored as characters.)
# If no, read in all data with map_df.
if("test_stat" %in% colnames(file_1)) {

  file_list <- vector("list", length(filenames))

  for(i in 1:length(filenames)) {
    file_list[[i]] <- data.table::fread(file = filenames[i])
    file_list[[i]]$test_stat <- as.numeric(file_list[[i]]$test_stat)
  }

  all_data <- do.call("rbind", file_list)

} else if("avg_ts" %in% colnames(file_1)) {

  file_list <- vector("list", length(filenames))

  for(i in 1:length(filenames)) {
    file_list[[i]] <- data.table::fread(file = filenames[i])
    file_list[[i]]$avg_ts <- as.numeric(file_list[[i]]$avg_ts)
  }

  all_data <- do.call("rbind", file_list)

} else {

  all_data <- list.files(path = "sim_data",
                         pattern = start_name, full.names = T) %>%
    purrr::map_df(~data.table::fread(.))

}

# Save new data frame
data.table::fwrite(all_data, file = "sim_data/fig01_fully_NP_randproj.csv")
