cat("Loading packages and paths\n")
library(cyCombine)
library(magrittr)
# Directory containing .fcs files
data_dir <- "/path/to/dir/fcs/"
# Markers of interest
#marker_names_df = read.csv("202008_.._XHLT2/data/cytof/pipeline_processed/singlet_FCSs/renamed/marker_metadata.csv")
# Removeing _ and - from marker names
#markers = gsub("[_-]","", marker_names_df$new_name_2)
fcs_file <- "/path/to/one/fcs/file.fcs"

fcs_data <-
  read.FCS(fcs_file) # Read FCS file and store in a temporary variable
fcs_param_data <- fcs_data@parameters@data  # Extract parameter data
marker_metadata <- fcs_param_data %>%
  mutate(desc = ifelse(is.na(desc), name, desc)) %>%   # Convert NA values to character "NA"
  mutate(marker_name_short = gsub("^[^_]+_", "", desc) ) %>%
  select(name, desc, marker_name_short) %>% setNames(c("channel_name", "marker_name", "marker_name_short"))  # Select channel name and marker name and rename the column names
markers = marker_metadata$marker_name 


# Data frame of FCS file names and batch ids
metadata = read.csv("/path/to/file_metadata.csv")

cat("Reading in the FCS files\n")
# Compile fcs files, down-sample, and preprocess
uncorrected <- prepare_data(data_dir = data_dir,
                             markers = markers,
                             metadata = metadata, # Can also be .csv file or data.frame object
                             sample_ids = NULL,
                             batch_ids = "pool_id",
                             filename_col = "file_name",
                             #condition = "Tumor_line",
                             down_sample = FALSE,
                             #sample_size = 500000,
                             seed = 473,
                             transform = TRUE,
                             cofactor = 5,
                             derand = TRUE) 

cat("Saving the uncorrected data RDS\n")
saveRDS(uncorrected, file = "cycombine_uncorrected.RDS")

cat("Running batch correction\n")
# Run batch correction
corrected <- uncorrected %>%
  batch_correct(markers = markers,
                xdim = 8,
                ydim = 8,
                norm_method = "scale", # "rank" is recommended when combining data with heavy batch effects
                rlen = 10, # Consider a larger value, if results are not convincing (e.g. 100)
                covar = NULL,
                parametric = TRUE)
cat("Saving the corrected data RDS\n")
saveRDS(corrected, file = "cycombine_corrected.RDS")



