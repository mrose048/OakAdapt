##%######################################################%##
#                                                          #
####        Process Sentinel-2 imagery - Randall        ####
#                                                          #
##%######################################################%##

# libraries
library(terra)
library(dplyr)
library(sf)
library(stringr)
library(fs)
library(utils)


##%######################################################%##
#                                                          #
####                  Randall Preserve                  ####
#                                                          #
##%######################################################%##


# Define paths
randall_safe_dir <- "H:/My Drive/16_Oak_Woodland/1_Data/Sentinel2_folder/02_Randall_unzip/"  # Change this to your directory containing .SAFE folders
randall_output_dir <- "H:/My Drive/16_Oak_Woodland/1_Data/Sentinel2_folder/03_Randall_envi/"     # Where to store ENVI files

# read a sentinel image file for crs
og_crs <- terra::rast(
  'H:/My Drive/16_Oak_Woodland/1_Data/Sentinel2_folder/02_Randall_unzip/S2A_MSIL2A_20150922T184436_N0500_R027_T11SLU_20231017T001227.SAFE/GRANULE/L2A_T11SLU_A001311_20150922T184430/IMG_DATA/R20m/T11SLU_20150922T184436_AOT_20m.jp2'
) %>%
  crs()

# randall shapefile
randall_shape <- terra::vect(
  'H:/My Drive/16_Oak_Woodland/1_Data/00_Shapefiles/TNC_Data/Protected_by_TNC_April_2022_dis/Protected_by_TNC_April_2022_dis.shp'
) %>%
  project(og_crs)


# --- FIND img & IMG_DATA FOLDERS ---
img_paths <- list.dirs(randall_safe_dir, recursive = TRUE) %>%
  grep("IMG_DATA$", ., value = TRUE)



##%######################################################%##
#                                                          #
####           Troubleshooting Randall issue            ####
#                                                          #
##%######################################################%##

#### Troubleshooting issue with repeated vegetation fraction across multiple time steps

# Convert to Date format
fraction_to_date <- function(frac_year) {
  year <- floor(frac_year)
  remainder <- frac_year - year
  days_in_year <- ifelse(leap_year(year), 366, 365)
  date <- as.Date(paste0(year, "-01-01")) + round(remainder * days_in_year)
  return(date)
}

# dates with same vegetation fraction after unmixing
bad_dates_df <- data.table::fread('bad_dates.csv') %>% tibble()

# Full list of bad dates
bad_dates <- bad_dates_df$date

# Convert fractional years to Dates
bad_date_objs <- sapply(bad_dates, fraction_to_date)

# Convert Dates to character format matching Sentinel filenames
bad_date_strings <- format(as.Date(bad_date_objs, origin = "1970-01-01"), "%Y%m%d")

# Create one regular expression pattern that matches any of the bad dates
pattern <- paste(bad_date_strings, collapse = "|")


# Filter matched files based on the combined pattern
bad_folders <- img_paths[grepl(pattern, img_paths)]

# Print matched bad file paths
print(bad_folders)

# call in a few of the raster files from the bad folders
bad_rast1 <- rast('Sentinel2_folder/02_Randall_unzip/S2A_MSIL2A_20221015T183341_N0400_R027_T11SLU_20221016T001059.SAFE/GRANULE/L2A_T11SLU_A038205_20221015T183634/IMG_DATA/R10m/T11SLU_20221015T183341_B03_10m.jp2')

bad_rast2 <- rast(
  'Sentinel2_folder/02_Randall_unzip/S2B_MSIL2A_20170708T184429_N0500_R027_T11SLU_20230819T121422.SAFE/GRANULE/L2A_T11SLU_A001769_20170708T184427/IMG_DATA/R10m/T11SLU_20170708T184429_B03_10m.jp2'
)

##%######################################################%##
#                                                          #
####           Processing Randall ENVI files            ####
#                                                          #
##%######################################################%##


# Initialize vector to store skipped files due to errors
skipped_files <- vector()
missing_bands <- vector()

# --- PROCESS EACH img ---
for (img in img_paths) {
  
  # --- DEFINE OUTPUT FILE PATH ---
  pattern <- "S2[AB]_MSIL2A_(\\d{8}T\\d{6}).*?_T(\\d{2}[A-Z]{3})"
  matches <- str_match(img, pattern)
  
  if (any(is.na(matches))) {
    message("Skipping ", img, " (could not extract tile/date)")
    skipped_files <<- c(skipped_files, img)  # Log skipped img
    next
  }
  
  date_str <- matches[1,2]  # Extract YYYYMMDDTHHMMSS
  tile_id <- matches[1,3]   # Extract TXXYYY Tile ID
  
  envi_filename <- paste0("S2_", tile_id, "_", date_str, "_L2A_e")
  envi_out <- file.path(randall_output_dir, envi_filename)
  
  # --- CHECK IF OUTPUT ALREADY EXISTS ---
  if (file.exists(envi_out)) {
    message("Skipping ", img, " (output file already exists: ", envi_out, ")")
    skipped_files <<- c(skipped_files, img)  # Log skipped img
    next
  }
  
  # --- FIND IMG_DATA DIRECTORY ---
  img_data_path <- list.dirs(img, recursive = TRUE) %>%
    grep("IMG_DATA", ., value = TRUE)
  
  if (length(img_data_path) == 0) {
    message("Skipping ", img, " (IMG_DATA folder not found)")
    skipped_files <<- c(skipped_files, img)  # Log skipped img
    next
  }
  
  # --- FIND JP2 FILES IN BOTH R10m AND R20m ---
  band_files_10m <- list.files(img_data_path, pattern = "_(B02|B03|B04|B08)_10m.jp2$", full.names = TRUE)
  band_files_20m <- list.files(img_data_path, pattern = "_(B01|B05|B06|B07|B8A|B11|B12)_20m.jp2$", full.names = TRUE)
  
  if (length(band_files_10m) == 0 || length(band_files_20m) == 0) {
    message("Skipping ", img, " (missing 10m or 20m bands)")
    skipped_files <<- c(skipped_files, img)  # Log skipped img
    next
  }
  
  # --- RESAMPLE 20M BANDS ---
  resampled_20m_bands <- lapply(band_files_20m, function(band) {
    band_raster <- tryCatch({
      rast(band)
    }, error = function(e) {
      message("Error reading 20m band: ", band)
      missing_bands <<- c(missing_bands, band)
      NULL
    })
    
    if (!is.null(band_raster)) {
      resampled_band <- tryCatch({
        resample(band_raster, rast(band_files_10m[2]), method = "bilinear")
      }, error = function(e) {
        message("Error resampling 20m band: ", band)
        missing_bands <<- c(missing_bands, band)
        NULL
      })
      return(resampled_band)
    } else {
      return(NULL)
    }
  })
  
  # Remove NULL values
  resampled_20m_bands <- resampled_20m_bands[!sapply(resampled_20m_bands, is.null)]
  
  # --- COMBINE BANDS ---
  band_files <- c(
    resampled_20m_bands[grepl("B01", band_files_20m)],
    lapply(band_files_10m, function(band) tryCatch(rast(band), error = function(e) NULL)),
    resampled_20m_bands[grepl("B05", band_files_20m)],
    resampled_20m_bands[grepl("B06", band_files_20m)],
    resampled_20m_bands[grepl("B07", band_files_20m)],
    resampled_20m_bands[grepl("B8A", band_files_20m)],
    resampled_20m_bands[grepl("B11", band_files_20m)],
    resampled_20m_bands[grepl("B12", band_files_20m)]
  )
  
  # Remove NULL values
  band_files <- band_files[!sapply(band_files, is.null)]
  
  if (length(band_files) != 11) {
    message("Skipping ", img, " (missing bands, found ", length(band_files), ")")
    skipped_files <<- c(skipped_files, img)
    next
  }
  
  # --- STACK & EXPORT ---
  tryCatch({
    s2_stack <- rast(band_files)
    s2_crop <- crop(s2_stack, randall_shape)
    writeRaster(s2_crop, filename = envi_out, filetype = "ENVI", datatype = "FLT4S", overwrite = TRUE)
    
    # --- ADD METADATA ---
    # Define wavelength metadata
    wavelength_units <- "wavelength units = Micrometers"
    wavelength_values <- "wavelength = {0.443, 0.490, 0.560, 0.665, 0.705, 0.740, 0.783, 0.842, 0.865, 1.610, 2.190}"
    
    # Define HDR file path
    hdr_file <- file.path(randall_output_dir, paste0(basename(envi_out), ".hdr"))
    hdr_file <- normalizePath(hdr_file, winslash = "/", mustWork = FALSE)
    
    # Read existing .hdr file if it exists
    if (file.exists(hdr_file)) {
      hdr_content <- readLines(hdr_file)
    } else {
      hdr_content <- c()
    }
    
    # Check if wavelength metadata already exists
    has_wavelength_units <- any(grepl("^wavelength units", hdr_content))
    has_wavelength_values <- any(grepl("^wavelength =", hdr_content))
    
    # Append only if missing
    if (!has_wavelength_units) {
      hdr_content <- c(hdr_content, wavelength_units)
    }
    if (!has_wavelength_values) {
      hdr_content <- c(hdr_content, wavelength_values)
    }
    
    # Write back to the HDR file
    writeLines(hdr_content, hdr_file)
    
  }, error = function(e) {
    message("Error processing img ", img, ": ", e$message)
    skipped_files <<- c(skipped_files, img)
  })
  
  gc()
}

# Summary of skipped and missing bands
message("Skipped images due to errors: ", paste(skipped_files, collapse = ", "))
message("Missing bands: ", paste(missing_bands, collapse = ", "))

# select .jp2 files from skipped_files
skipped_files_jp2 <- lapply(skipped_files, function(x) {
  # Extract the IMG_DATA path
  img_data_path <- list.dirs(x, recursive = TRUE) %>%
    grep("IMG_DATA", ., value = TRUE)
  
  # Find all .jp2 files in the IMG_DATA directory
  jp2_files <- list.files(img_data_path, pattern = "\\.jp2$", full.names = TRUE)
  
  return(jp2_files)
})



##%######################################################%##
#                                                          #
####                     Dye Creek                      ####
#                                                          #
##%######################################################%##

# Define paths
dye_safe_dir <- "H:/My Drive/16_Oak_Woodland/1_Data/Sentinel2_folder/02_DyeCreek_unzip/"  # Change this to your directory containing .SAFE folders
dye_output_dir <- "H:/My Drive/16_Oak_Woodland/1_Data/Sentinel2_folder/03_DyeCreek_envi/"     # Where to store ENVI files

# read a sentinel image file for crs
og_crs <- terra::rast(
  'H:/My Drive/16_Oak_Woodland/1_Data/Sentinel2_folder/02_DyeCreek_unzip/S2A_MSIL2A_20160425T184922_N0500_R113_T10TEK_20231020T025906.SAFE/GRANULE/L2A_T10TEK_A004400_20160425T185405/IMG_DATA/R20m/T10TEK_20160425T184922_AOT_20m.jp2'
) %>%
  crs()

# dye creek shapefile
dye_shape <- terra::vect('H:/My Drive/16_Oak_Woodland/1_Data/00_Shapefiles/TNC_Data/Dye_Creek_Preserve_clean_proj.shp')  %>%
  project(og_crs)

# --- FIND img & IMG_DATA FOLDERS ---
img_paths <- list.dirs(dye_safe_dir, recursive = TRUE) %>%
  grep("IMG_DATA$", ., value = TRUE)

# Initialize vector to store skipped files due to errors
skipped_files <- vector()
missing_bands <- vector()

# --- PROCESS EACH img ---
for (img in img_paths) {
  
  # --- DEFINE OUTPUT FILE PATH ---
  pattern <- "S2[AB]_MSIL2A_(\\d{8}T\\d{6}).*?_T(\\d{2}[A-Z]{3})"
  matches <- str_match(img, pattern)
  
  if (any(is.na(matches))) {
    message("Skipping ", img, " (could not extract tile/date)")
    skipped_files <<- c(skipped_files, img)  # Log skipped img
    next
  }
  
  date_str <- matches[1,2]  # Extract YYYYMMDDTHHMMSS
  tile_id <- matches[1,3]   # Extract TXXYYY Tile ID
  
  envi_filename <- paste0("S2_", tile_id, "_", date_str, "_L2A_e")
  envi_out <- file.path(dye_output_dir, envi_filename)
  
  # --- CHECK IF OUTPUT ALREADY EXISTS ---
  if (file.exists(envi_out)) {
    message("Skipping ", img, " (output file already exists: ", envi_out, ")")
    skipped_files <<- c(skipped_files, img)  # Log skipped img
    next
  }
  
  # --- FIND IMG_DATA DIRECTORY ---
  img_data_path <- list.dirs(img, recursive = TRUE) %>%
    grep("IMG_DATA", ., value = TRUE)
  
  if (length(img_data_path) == 0) {
    message("Skipping ", img, " (IMG_DATA folder not found)")
    skipped_files <<- c(skipped_files, img)  # Log skipped img
    next
  }
  
  # --- FIND JP2 FILES IN BOTH R10m AND R20m ---
  band_files_10m <- list.files(img_data_path, pattern = "_(B02|B03|B04|B08)_10m.jp2$", full.names = TRUE)
  band_files_20m <- list.files(img_data_path, pattern = "_(B01|B05|B06|B07|B8A|B11|B12)_20m.jp2$", full.names = TRUE)
  
  if (length(band_files_10m) == 0 || length(band_files_20m) == 0) {
    message("Skipping ", img, " (missing 10m or 20m bands)")
    skipped_files <<- c(skipped_files, img)  # Log skipped img
    next
  }
  
  # --- RESAMPLE 20M BANDS ---
  resampled_20m_bands <- lapply(band_files_20m, function(band) {
    band_raster <- tryCatch({
      rast(band)
    }, error = function(e) {
      message("Error reading 20m band: ", band)
      missing_bands <<- c(missing_bands, band)
      NULL
    })
    
    if (!is.null(band_raster)) {
      resampled_band <- tryCatch({
        resample(band_raster, rast(band_files_10m[2]), method = "bilinear")
      }, error = function(e) {
        message("Error resampling 20m band: ", band)
        missing_bands <<- c(missing_bands, band)
        NULL
      })
      return(resampled_band)
    } else {
      return(NULL)
    }
  })
  
  # Remove NULL values
  resampled_20m_bands <- resampled_20m_bands[!sapply(resampled_20m_bands, is.null)]
  
  # --- COMBINE BANDS ---
  band_files <- c(
    resampled_20m_bands[grepl("B01", band_files_20m)],
    lapply(band_files_10m, function(band) tryCatch(rast(band), error = function(e) NULL)),
    resampled_20m_bands[grepl("B05", band_files_20m)],
    resampled_20m_bands[grepl("B06", band_files_20m)],
    resampled_20m_bands[grepl("B07", band_files_20m)],
    resampled_20m_bands[grepl("B8A", band_files_20m)],
    resampled_20m_bands[grepl("B11", band_files_20m)],
    resampled_20m_bands[grepl("B12", band_files_20m)]
  )
  
  # Remove NULL values
  band_files <- band_files[!sapply(band_files, is.null)]
  
  if (length(band_files) != 11) {
    message("Skipping ", img, " (missing bands, found ", length(band_files), ")")
    skipped_files <<- c(skipped_files, img)
    next
  }
  
  # --- STACK & EXPORT ---
  tryCatch({
    s2_stack <- rast(band_files)
    s2_crop <- crop(s2_stack, dye_shape)
    writeRaster(s2_crop, filename = envi_out, filetype = "ENVI", datatype = "FLT4S", overwrite = TRUE)
    
    # --- ADD METADATA ---
    # Define wavelength metadata
    wavelength_units <- "wavelength units = Micrometers"
    wavelength_values <- "wavelength = {0.443, 0.490, 0.560, 0.665, 0.705, 0.740, 0.783, 0.842, 0.865, 1.610, 2.190}"
    
    # Define HDR file path
    hdr_file <- file.path(dye_output_dir, paste0(basename(envi_out), ".hdr"))
    hdr_file <- normalizePath(hdr_file, winslash = "/", mustWork = FALSE)
    
    # Read existing .hdr file if it exists
    if (file.exists(hdr_file)) {
      hdr_content <- readLines(hdr_file)
    } else {
      hdr_content <- c()
    }
    
    # Check if wavelength metadata already exists
    has_wavelength_units <- any(grepl("^wavelength units", hdr_content))
    has_wavelength_values <- any(grepl("^wavelength =", hdr_content))
    
    # Append only if missing
    if (!has_wavelength_units) {
      hdr_content <- c(hdr_content, wavelength_units)
    }
    if (!has_wavelength_values) {
      hdr_content <- c(hdr_content, wavelength_values)
    }
    
    # Write back to the HDR file
    writeLines(hdr_content, hdr_file)
    
  }, error = function(e) {
    message("Error processing img ", img, ": ", e$message)
    skipped_files <<- c(skipped_files, img)
  })
  
  gc()
}

# Summary of skipped and missing bands
message("Skipped imgs due to errors: ", paste(skipped_files, collapse = ", "))
message("Missing bands: ", paste(missing_bands, collapse = ", "))

# Investigate missing files

files <- list.files('Sentinel2_folder/03_DyeCreek_envi', pattern = 'e$', full.names = F)
folders <- list.dirs('Sentinel2_folder/02_DyeCreek_unzip/', recursive = F, full.names = F)

# Extract tile + datetime from files
file_keys <- gsub("S2_(\\w+?)_(\\d{8}T\\d{6})_L2A.*", "\\1_\\2", files)

# Extract tile + datetime from folders
folder_keys <- gsub(".*_(\\d{8}T\\d{6})_.*_T(\\w+?)_.*", "\\2_\\1", folders)

# Folders that don't have a match in files
missing_folders <- folders[!(folder_keys %in% file_keys)]

# View missing folders
missing_folders
