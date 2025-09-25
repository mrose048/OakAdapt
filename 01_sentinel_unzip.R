##%######################################################%##
#                                                          #
####               Sentinel-2 Processing                ####
#                                                          #
##%######################################################%##

# libraries
library(terra)
library(dplyr)
library(sf)
library(stringr)
library(fs)
library(utils)


# Function to unzip all Sentinel-2 .zip files
unzip_safe_folders <- function(zip_dir, output_dir) {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # List all .zip files in the directory
  zip_files <- list.files(zip_dir, pattern = "\\.zip$", full.names = TRUE, recursive = TRUE)
  
  # Loop through each zip file
  for (zip_file in zip_files) {
    # Define expected output folder name
    safe_folder <- file.path(output_dir, gsub("\\.zip$", "", basename(zip_file)))
    
    # Check if the folder already exists
    if (dir.exists(safe_folder)) {
      message(paste("Skipping:", basename(zip_file), "- Already unzipped"))
    } else {
      message(paste("Unzipping:", basename(zip_file)))
      unzip(zip_file, exdir = output_dir)
    }
  }
}

##%######################################################%##
#                                                          #
####                  Unzip - Randall                   ####
#                                                          #
##%######################################################%##

# Randall zipped files and directory for unzipped files
zip_directory <- "H:/My Drive/16_Oak_Woodland/1_Data/Sentinel2_folder/01_Randall_raw/All Sentinel Downloads 3 31/"       # Change this to the folder containing Sentinel-2 .zip files
output_directory <- "H:/My Drive/16_Oak_Woodland/1_Data/Sentinel2_folder/02_Randall_unzip/"  # Change this to where you want to extract the files

unzip_safe_folders(zip_directory, output_directory)

##%######################################################%##
#                                                          #
####                 Unzip - Dye Creek                  ####
#                                                          #
##%######################################################%##

# Dye Creek zipped files and directory for unzipped files
zip_directory <- "H:/My Drive/16_Oak_Woodland/1_Data/Sentinel2_folder/01_DyeCreek_raw/"       # Change this to the folder containing Sentinel-2 .zip files
output_directory <- "H:/My Drive/16_Oak_Woodland/1_Data/Sentinel2_folder/02_DyeCreek_unzip/"  # Change this to where you want to extract the files

unzip_safe_folders(zip_directory, output_directory)


