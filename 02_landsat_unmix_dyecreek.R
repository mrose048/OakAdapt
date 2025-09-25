##%######################################################%##
#                                                          #
####       Mixture Residual Analysis - Dye Creek        ####
#                                                          #
##%######################################################%##

#Load libraries

library(caTools)
library(reticulate)
require(terra)
require(dplyr)
require(sf)
require(tidyverse)
require(ggspatial)

setwd('H:/My Drive/16_Oak_Woodland/')

# original crs
og_crs <- terra::rast('1_Data/Landsat/03_DyeCreek_mosaics/LC08_L2SP_044032_20130603_20200912_02_T1.tif') %>% crs()

# Dye Creek shapefile with correct projected coordinate system -- always check!
shape <- vect('1_Data/00_Shapefiles/TNC_Data/Dye_Creek_Preserve_clean_proj.shp') 


##%######################################################%##
#                                                          #
####         Landsat 8 & 9 - Spectral Unmixing          ####
#                                                          #
##%######################################################%##

#Set filename for data. endmembers, and bad bands
# data_file = '03_DyeCreek_mosaics/LC08_L2SP_041036_20130427_20200912_02_T1_Scale' # figure out what this data file needs to be, ENVI file for?
endmember_file = '00_Endmembers/EM_MR_L8_txtNH.txt'
badbands_file = '00_Badbands/L8bbl.csv'

#Set unit sum constraint weight
w = 1

#Multiply EMs by a scale factor?
rescale = 1/10000

files <- list.files('03_DyeCreek_mosaics/', pattern = 'T1_Scale$', full.names = TRUE)

# select only the LC files
files <- files[grep("LC", files)]

outfiles = paste('04_DyeCreek_unmixed/', basename(files),'_Unmixed',sep="")

#Read bad bands list
bands = read.csv(badbands_file, header = F, sep = ",", stringsAsFactors = F)

#Read endmember file
G = read.csv(endmember_file, header = F, sep = "")

for(i in 1:length(files)) {
  
  print(i)
  
  if(!file.exists(outfiles[i])) {
    
    # Read spatial data file
    D1 <- rast(files[i])  # Read as a raster
    crs_D1 <- crs(D1)  # Store the original CRS
    
    # Get wavelengths from the first column of the endmember file
    wavelengths <- G[, 1]
    G2 <- G[, 2:ncol(G)]
    
    # Rescale endmembers
    G3 <- G2 * rescale
    
    # Subset using bad bands
    D2 <- as.array(D1)[, , bands == 1]  # Extract valid bands
    G4 <- G3[bands == 1, ]
    
    # Reshape arrays
    D3 <- array_reshape(D2, c(dim(D2)[1] * dim(D2)[2], dim(D2)[3]), order = "C")
    
    # Add unit sum constraint
    D4 <- t(as.matrix(cbind(D3, rep(1, nrow(D3)))))  # Transpose to add unit sum constraint
    G5 <- as.matrix(rbind(G4, rep(1, ncol(G4))))
    
    # Compute residual
    r <- D4 - (G5 %*% solve(t(G5) %*% G5) %*% t(G5) %*% D4)
    
    # Compute mixture fractions
    u <- solve(t(G5) %*% G5) %*% t(G5) %*% D4
    
    # Transpose the matrix to match D3 format
    u2 <- t(u)
    
    # Reshape back to original dimensions
    u_reshape <- array_reshape(u2, c(dim(D1)[1], dim(D1)[2], dim(u2)[2]), order = "C")
    
    # Convert back to raster while preserving CRS
    u_raster <- rast(u_reshape, ext = ext(D1), crs = crs_D1)
    
    # crop and mask to the extent of DyeCreek preserve boundaries
    u_raster <- u_raster %>%
      crop(shape) %>%
      mask(shape)
    
    # write to file
    terra::writeRaster(
      u_raster,
      outfiles[i],
      overwrite = TRUE,
      datatype = "FLT4S",
      filetype = 'ENVI'
    )
    
    #Export result
    # write.ENVI(u_reshape, outfiles[i], interleave = c('bsq'))
    
    gc()
    
  }
}


##%######################################################%##
#                                                          #
####           Landsat 7 - Spectral Unmixing            ####
#                                                          #
##%######################################################%##

endmember_file = '00_Endmembers/EM_MR_L7_txtNH.txt'
badbands_file = '00_Badbands/L7bbl.csv'

#Set unit sum constraint weight
w = 1

#Multiply EMs by a scale factor?
rescale = 1/10000

files <- list.files('03_DyeCreek_mosaics/', pattern = 'T1_Scale$', full.names = TRUE)

# select only files with LE07 pattern
files <- files[grepl('LE07', files)]

outfiles = paste('04_DyeCreek_unmixed/', basename(files),'_Unmixed',sep="")

#Read bad bands list
bands = read.csv(badbands_file, header = F, sep = ",", stringsAsFactors = F)

#Read endmember file
G = read.csv(endmember_file, header = F, sep = "")

for(i in 1:length(files)) {
  
  print(i)
  
  if(!file.exists(outfiles[i])) {
    
    # Read spatial data file
    D1 <- rast(files[i])  # Read as a raster
    crs_D1 <- crs(D1)  # Store the original CRS
    
    # Get wavelengths from the first column of the endmember file
    wavelengths <- G[, 1]
    G2 <- G[, 2:ncol(G)]
    
    # Rescale endmembers
    G3 <- G2 * rescale
    
    # Subset using bad bands
    D2 <- as.array(D1)[, , bands == 1]  # Extract valid bands
    G4 <- G3[bands == 1, ]
    
    # Reshape arrays
    D3 <- array_reshape(D2, c(dim(D2)[1] * dim(D2)[2], dim(D2)[3]), order = "C")
    
    # Add unit sum constraint
    D4 <- t(as.matrix(cbind(D3, rep(1, nrow(D3)))))  # Transpose to add unit sum constraint
    G5 <- as.matrix(rbind(G4, rep(1, ncol(G4))))
    
    # Compute residual
    r <- D4 - (G5 %*% solve(t(G5) %*% G5) %*% t(G5) %*% D4)
    
    # Compute mixture fractions
    u <- solve(t(G5) %*% G5) %*% t(G5) %*% D4
    
    # Transpose the matrix to match D3 format
    u2 <- t(u)
    
    # Reshape back to original dimensions
    u_reshape <- array_reshape(u2, c(dim(D1)[1], dim(D1)[2], dim(u2)[2]), order = "C")
    
    # Convert back to raster while preserving CRS
    u_raster <- rast(u_reshape, ext = ext(D1), crs = crs_D1)
    
    # crop and mask to the extent of DyeCreek preserve boundaries
    u_raster <- u_raster %>%
      crop(shape) %>%
      mask(shape)
    
    # write to file
    terra::writeRaster(
      u_raster,
      outfiles[i],
      overwrite = TRUE,
      datatype = "FLT4S",
      filetype = 'ENVI'
    )
    
    #Export result
    # write.ENVI(u_reshape, outfiles[i], interleave = c('bsq'))
  }
}


##%######################################################%##
#                                                          #
####           Landsat 5 - Spectral Unmixing            ####
#                                                          #
##%######################################################%##

# problem files: 
endmember_file = '00_Endmembers/EM_MR_L5_txtNH.txt'
badbands_file = '00_Badbands/L5bbl.csv'

#Set unit sum constraint weight
w = 1

#Multiply EMs by a scale factor?
rescale = 1/10000

files <- list.files('03_DyeCreek_mosaics/', pattern = 'T1_Scale$', full.names = TRUE)

# select only files with LT pattern
files <- files[grepl('LT05', files)]

outfiles = paste('04_DyeCreek_unmixed/', basename(files),'_Unmixed',sep="")

#Read bad bands list
bands = read.csv(badbands_file, header = F, sep = ",", stringsAsFactors = F)

#Read endmember file
G = read.csv(endmember_file, header = F, sep = "")

for(i in 1:length(files)) {
  
  print(i)
  
  if(!file.exists(outfiles[i])) {
    
    # Read spatial data file
    D1 <- rast(files[i])  # Read as a raster
    crs_D1 <- crs(D1)  # Store the original CRS
    
    # Get wavelengths from the first column of the endmember file
    wavelengths <- G[, 1]
    G2 <- G[, 2:ncol(G)]
    
    # Rescale endmembers
    G3 <- G2 * rescale
    
    # Subset using bad bands
    D2 <- as.array(D1)[, , bands == 1]  # Extract valid bands
    G4 <- G3[bands == 1, ]
    
    # Reshape arrays
    D3 <- array_reshape(D2, c(dim(D2)[1] * dim(D2)[2], dim(D2)[3]), order = "C")
    
    # Add unit sum constraint
    D4 <- t(as.matrix(cbind(D3, rep(1, nrow(D3)))))  # Transpose to add unit sum constraint
    G5 <- as.matrix(rbind(G4, rep(1, ncol(G4))))
    
    # Compute residual
    r <- D4 - (G5 %*% solve(t(G5) %*% G5) %*% t(G5) %*% D4)
    
    # Compute mixture fractions
    u <- solve(t(G5) %*% G5) %*% t(G5) %*% D4
    
    # Transpose the matrix to match D3 format
    u2 <- t(u)
    
    # Reshape back to original dimensions
    u_reshape <- array_reshape(u2, c(dim(D1)[1], dim(D1)[2], dim(u2)[2]), order = "C")
    
    # Convert back to raster while preserving CRS
    u_raster <- rast(u_reshape, ext = ext(D1), crs = crs_D1)
    
    # crop and mask to the extent of DyeCreek preserve boundaries
    u_raster <- u_raster %>%
      crop(shape) %>%
      mask(shape)
    
    # write to file
    terra::writeRaster(
      u_raster,
      outfiles[i],
      overwrite = TRUE,
      datatype = "FLT4S",
      filetype = 'ENVI'
    )
    
    #Export result
    # write.ENVI(u_reshape, outfiles[i], interleave = c('bsq'))
    
  }
}

##%######################################################%##
#                                                          #
####           Landsat 4 - Spectral Unmixing            ####
#                                                          #
##%######################################################%##

endmember_file = '00_Endmembers/EM_MR_L4_txtNH.txt'
badbands_file = '00_Badbands/L5bbl.csv'

#Set unit sum constraint weight
w = 1

#Multiply EMs by a scale factor?
rescale = 1/10000

files <- list.files('03_DyeCreek_mosaics/', pattern = 'T1_Scale$', full.names = TRUE)

# select only files with LT pattern
files <- files[grepl('LT04', files)]

outfiles = paste('04_DyeCreek_unmixed/', basename(files),'_Unmixed',sep="")

#Read bad bands list
bands = read.csv(badbands_file, header = F, sep = ",", stringsAsFactors = F)

#Read endmember file
G = read.csv(endmember_file, header = F, sep = "")

for(i in 1:length(files)) {
  
  print(i)
  
  if(!file.exists(outfiles[i])) {
    
    # Read spatial data file
    D1 <- rast(files[i])  # Read as a raster
    crs_D1 <- crs(D1)  # Store the original CRS
    
    # Get wavelengths from the first column of the endmember file
    wavelengths <- G[, 1]
    G2 <- G[, 2:ncol(G)]
    
    # Rescale endmembers
    G3 <- G2 * rescale
    
    # Subset using bad bands
    D2 <- as.array(D1)[, , bands == 1]  # Extract valid bands
    G4 <- G3[bands == 1, ]
    
    # Reshape arrays
    D3 <- array_reshape(D2, c(dim(D2)[1] * dim(D2)[2], dim(D2)[3]), order = "C")
    
    # Add unit sum constraint
    D4 <- t(as.matrix(cbind(D3, rep(1, nrow(D3)))))  # Transpose to add unit sum constraint
    G5 <- as.matrix(rbind(G4, rep(1, ncol(G4))))
    
    # Compute residual
    r <- D4 - (G5 %*% solve(t(G5) %*% G5) %*% t(G5) %*% D4)
    
    # Compute mixture fractions
    u <- solve(t(G5) %*% G5) %*% t(G5) %*% D4
    
    # Transpose the matrix to match D3 format
    u2 <- t(u)
    
    # Reshape back to original dimensions
    u_reshape <- array_reshape(u2, c(dim(D1)[1], dim(D1)[2], dim(u2)[2]), order = "C")
    
    # Convert back to raster while preserving CRS
    u_raster <- rast(u_reshape, ext = ext(D1), crs = crs_D1)
    
    # crop and mask to the extent of DyeCreek preserve boundaries
    u_raster <- u_raster %>%
      crop(shape) %>%
      mask(shape)
    
    # write to file
    terra::writeRaster(
      u_raster,
      outfiles[i],
      overwrite = TRUE,
      datatype = "FLT4S",
      filetype = 'ENVI'
    )
    
    #Export result
    # write.ENVI(u_reshape, outfiles[i], interleave = c('bsq'))
    
  }
}


##%######################################################%##
#                                                          #
####                   File renaming                    ####
#                                                          #
##%######################################################%##


# rename files in '04_DyeCreek_unmixed/' so that the first string of characters is added to the end of the file

# TODO fix so that file extension does not get messed up!
rename_files <- function(directory) {
  # Get list of files in the directory
  files <- list.files(directory, full.names = TRUE)
  
  for (file in files) {
    # Extract the filename without the directory path
    filename <- basename(file)
    
    # Ensure it follows the expected Landsat naming format
    parts <- unlist(strsplit(filename, "_"))
    
    if (length(parts) < 8) {
      next  # Skip files that don't match the expected format
    }
    
    # Extract components
    sensor_mission <- parts[1]  # e.g., LC08
    product <- parts[2]          # e.g., L2SP
    path_row <- parts[3]         # e.g., 041036
    date1 <- parts[4]            # e.g., 20130630
    date2 <- parts[5]            # e.g., 20200912
    version <- parts[6]          # e.g., 02
    tier <- parts[7]             # e.g., T1
    scale_unmixed <- paste(parts[8:length(parts)], collapse = "_")  # e.g., Scale_Unmixed
    
    # New filename
    new_filename <- paste(product, date1, date2, version, tier, sensor_mission, path_row, scale_unmixed, sep = "_")
    
    # Construct full file paths
    new_filepath <- file.path(directory, new_filename)
    
    # Rename the file
    file.rename(file, new_filepath)
    
    # Print confirmation
    cat("Renamed:", filename, "->", new_filename, "\n")
  }
}



# Set directory path where the files are stored
directory_path <- "04_DyeCreek_unmixed/"  # Change this to your actual directory

# Run the function
rename_files(directory_path)

##%######################################################%##
#                                                          #
####                   Layer stacking                   ####
#                                                          #
##%######################################################%##

library(lubridate)
library(stringr)

# clear image list
a_files <- list.files(
  path = "H:/My Drive/16_Oak_Woodland/1_Data/Landsat/TNC_Landsat_Paths/Dye_Creek/A_new",
  pattern = ".jpg",
  full.names = FALSE
) %>%
  strsplit("_") %>%
  sapply(function(x) paste(x[1], x[4], sep = "_")) %>%
  unique()


# Stack all files in "04_DyeCreek_unmixed/" that were in the original clear folder and save as Envi raster
files <- list.files(
  path = "1_Data/Landsat/04_DyeCreek_unmixed/",
  pattern = "Unmixed$",
  full.names = TRUE
)

# Extract satellite and date information from the file names 
file_info <- basename(files) %>%
  strsplit("_") %>%
  sapply(function(x) paste(x[6], x[2], sep = "_"))

# Select files where the satellite and image date match those in a_files
matching_files <- files[file_info %in% a_files | grepl("2025", file_info)]

# write matching files to file
writeLines(matching_files, "1_Data/Landsat/DyeCreek_LS_dates.txt")

# derive fractional year from file names
calc_frac_year <- function(filename) {
  match <- regmatches(filename, regexpr('\\d{8}', filename))
  
  if (length(match) == 0) {
    warning(paste('No date found in filename:', filename))
    return(NA)
  }
  # date object
  date_obj <- ymd(match)
  
  #Extract year and day of year
  year <- year(date_obj)
  doy <- yday(date_obj)
  
  # Check if it's a leap year
  days_in_year <- ifelse(leap_year(year), 366, 365)
  
  # Compute fractional year
  frac_year <- year + (doy - 1) / days_in_year
  
  return(frac_year)
}

# Use calc_frac_years on files
frac_years <- calc_frac_year(files)

# Read the rasters
r <- lapply(files, rast)

# stack the 2nd layer in each element of r (vegetation fraction)
r2 <- lapply(r, function(x) x[[2]])

# stack the rasters
veg <- rast(r2)
names(veg) <- frac_years

# remove problem dates (286, 209, 208, 26, 633)
# "1999.99178082192", "1995.78630136986", "1995.74246575342", "1986.45479452055", "2016.3306010929"
veg2 <- veg %>% 
  subset(names(veg) != 1999.99178082192) %>% 
  subset(names(veg) != 1995.78630136986) %>% 
  subset(names(veg) != 1995.74246575342) %>% 
  subset(names(veg) != 1986.45479452055) %>% 
  subset(names(veg) != 2016.3306010929)

# mask of NA values in veg - this really truncates the study extent in the east, not sure what is happening
mask <- any(is.na(veg2))

writeRaster(
  mask,
  '05_DyeCreek_temporal_stack/DyeCreek_Landsat_mask',
  overwrite = TRUE,
  filetype = 'ENVI'
)

# # mask out NA values
masked_veg <- mask(veg2, mask, maskvalues = 1)


# # Save as ENVI file
# writeRaster(
#   masked_veg,
#   '05_DyeCreek_temporal_stack/DyeCreek_Landsat_1982_2024_masked',
#   datatype = "FLT4S",
#   overwrite = TRUE,
#   filetype = 'ENVI'
# )

# Remove agricultural areas from the stack
nlcd_dye <- rast('00_Shapefiles/NLCD_mask/dye_nlcd.tif')

# Create a mask for **non-agricultural and developed land**
non_ag_mask <- classify(nlcd_dye, cbind(c(81, 82, 21, 22, 23, 24), NA), others=1) 

# Crop non_ag_mask to match the sentinel imagery extent
non_ag_mask <- non_ag_mask %>% crop(masked_veg)

# Resample non_ag_mask to match s2 resolution (without reprojecting)
non_ag_mask <- resample(non_ag_mask, masked_veg, method="near")

# mask out non-agricultural areas
masked_veg_non_ag <- mask(masked_veg, non_ag_mask)

# elevation rasters
elev_rast <- list.files('00_DEM/', pattern = '.tif$', full.names = TRUE) %>%
  lapply(rast)

elev_rast <- do.call(terra::mosaic, elev_rast) %>%
  project(masked_veg[[1]]) %>%
  resample(masked_veg[[1]], method = 'bilinear')

high_elev <- elev_rast >100
high_elev[high_elev == 0] <- NA

# mask based on DEM
veg_final <- mask(masked_veg_non_ag, high_elev)


# write vegetation stack without ag to file
writeRaster(veg_final,
            '05_DyeCreek_temporal_stack/DyeCreek_LS_1982_2025_veg_stack',
            datatype = "FLT4S",
            overwrite = TRUE,
            filetype = 'ENVI')


# Change wavelength names in ENVI hdr file
hdr_file <- '05_DyeCreek_temporal_stack/DyeCreek_LS_1982_2025_veg_stack.hdr'
# 
# # Read the header file
hdr <- readLines(hdr_file)
# 
# # Add a line for wavelengths = {
hdr <- c(hdr[1:grep('data ignore value = nan', hdr)], 'wavelength = {')
# 
# # repeat lines 14 through 863
for (i in 15:863) {
  hdr <- c(hdr, hdr[i])
}
# 
# # write hdr file
writeLines(hdr,
           '05_DyeCreek_temporal_stack/DyeCreek_LS_1982_2025_veg_stack.hdr')


##%######################################################%##
#                                                          #
####              Temporal Mixture Models               ####
#                                                          #
##%######################################################%##

# Temporal endmembers
endmember_file = '1_Data/Landsat/05_DyeCreek_temporal_stack/DyeCreek_LS_1982_2025_TEMs.txt'

#Set unit sum constraint weight
w = 1

#Multiply EMs by a scale factor?
# rescale = 1/10000

# Create a dataframe with 837 columns with 1
bands = data.frame(matrix(1, nrow = 1, ncol = 832))

#Read endmember file
G = read.csv(endmember_file, header = F, sep = "", skip = 6)

# remove last three rows (May 2025) for consistency with Randall data
G <- G[1:(nrow(G) - 3), ]


ggplot(G) +
  geom_line(aes(x = V1, y = V3), color = 'darkgreen') +
  geom_line(aes(x = V1, y = V4), color = 'red') +
  geom_line(aes(x = V1, y = V5), color = 'blue') +
  theme_minimal()

# Figure for paper
library(ggplot2)

# Reshape the data for cleaner ggplot code
library(tidyr)
G_long <- pivot_longer(G, cols = V2:V5, names_to = "class", values_to = "veg_fraction") %>%
  filter(class != 'V2')
frac_year_to_date <- function(frac_year) {
  year <- floor(frac_year)
  remainder <- frac_year - year
  as.Date(paste0(year, "-01-01")) + round(remainder * 365.25)
}

# Apply to your column
G_long$date <- frac_year_to_date(G_long$V1)
G_long$year <- format(G_long$date, "%Y")
G_long$doy <- as.numeric(format(G_long$date, "%j"))

get_water_year <- function(date) {
  year <- as.numeric(format(date, "%Y"))
  month <- as.numeric(format(date, "%m"))
  ifelse(month >= 10, year + 1, year)
}

G_long$water_year <- get_water_year(G_long$date)

max_dates_by_class <- G_long %>%
  group_by(class, water_year) %>%
  filter(veg_fraction == max(veg_fraction, na.rm = TRUE)) %>%
  slice(1) %>%  # in case of ties
  ungroup()

# Get DOY relative to water year (Oct 1 = DOY 1)
G_long$wy_doy <- as.numeric(G_long$date - as.Date(paste0(G_long$water_year - 1, "-10-01"))) + 1

max_dates_by_class$wy_doy <- as.numeric(max_dates_by_class$date - as.Date(paste0(
  max_dates_by_class$water_year - 1, "-10-01"
))) + 1

class_means <- max_dates_by_class %>%
  group_by(class) %>%
  summarise(mean_peak = mean(wy_doy, na.rm = TRUE))

ggplot(max_dates_by_class, aes(x = water_year, y = wy_doy, color = class)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  # Add average peak timing per class
  geom_hline(data = class_means, aes(yintercept = mean_peak, color = class),
             linetype = "dashed", size = 1) +
  scale_color_manual(
    values = c(
      "V3" = "#39FF14",
      "V5" = "#00FFFF",
      "V4" = "red"
    )
  ) +
  labs(
    x = "Water Year",
    y = "Day of Water Year \n of Max Vegetation",
    color = "Vegetation Type",
    title = ""
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_minimal(base_size = 18) +
  theme(
    plot.background = element_rect(fill = "#1a1a1a", color = NA),
    # dark gray background
    panel.background = element_rect(fill = "#1a1a1a", color = NA),
    panel.grid.major = element_line(color = "#333333"),
    # subtle grid lines
    panel.grid.minor = element_line(color = "#2a2a2a"),
    axis.text = element_text(color = "white"),
    axis.title = element_text(
      color = "white",
      size = 18,
      face = "bold"
    ),
    legend.background = element_rect(fill = "#1a1a1a"),
    legend.key = element_rect(fill = "#1a1a1a"),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white", face = "bold"),
    legend.position = "none"
  )

ggplot(G_long, aes(x = V1, y = veg_fraction, color = class)) +
  geom_line(size = .65) +
  scale_color_manual(
    values = c(
      V3 = "#39FF14",
      V4 = "red",
      V5 = "cyan"
    ),
    labels = c("Evergreen", "Annual", "Deciduous") # customize labels
  ) +
  labs(
    x = "Year",
    # replace with your actual variable
    y = "Vegetation Fraction",
    # replace with your actual variable
    color = ""
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.background = element_rect(fill = "#1a1a1a", color = NA),
    # dark gray background
    panel.background = element_rect(fill = "#1a1a1a", color = NA),
    panel.grid.major = element_line(color = "#333333"),
    # subtle grid lines
    panel.grid.minor = element_line(color = "#2a2a2a"),
    axis.text = element_text(color = "white"),
    axis.title = element_text(
      color = "white",
      size = 14,
      face = "bold"
    ),
    legend.background = element_rect(fill = "#1a1a1a"),
    legend.key = element_rect(fill = "#1a1a1a"),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white", face = "bold"),
    legend.position = "none"
  ) + annotate(
    "rect",
    xmin = 2012,
    xmax = 2015.7,
    ymin = -Inf,
    ymax = Inf,
    fill = "red",
    alpha = 0.25
  ) +
  annotate(
    "rect",
    xmin = 2020,
    xmax = 2022,
    ymin = -Inf,
    ymax = Inf,
    fill = "red",
    alpha = 0.25
  )

# Read spatial data file
D1 <- rast('1_Data/Landsat/05_DyeCreek_temporal_stack/DyeCreek_LS_1982_2025_veg_stack')  # Read as a raster
crs_D1 <- crs(D1)  # Store the original CRS

# Remove May imags
D1 <- D1[[1:(nlyr(D1) - 3)]]  # Exclude the last three layers (May 2025)

# Get wavelengths from the first column of the endmember file
wavelengths <- G[, 1]
G2 <- G[, 2:ncol(G)]

# Rescale endmembers
# G3 <- G2 * rescale
G3 <- G2

# Subset using bad bands
D2 <- as.array(D1)[, , bands == 1]  # Extract valid bands
G4 <- G3[bands == 1, ]

# Reshape arrays
D3 <- array_reshape(D2, c(dim(D2)[1] * dim(D2)[2], dim(D2)[3]), order = "C")

# Add unit sum constraint
D4 <- t(as.matrix(cbind(D3, rep(1, nrow(
  D3
)))))  # Transpose to add unit sum constraint
G5 <- as.matrix(rbind(G4, rep(1, ncol(G4))))

# Compute residual of time series - you will transform this with PCA and evaluate the first three PC
r <- D4 - (G5 %*% solve(t(G5) %*% G5) %*% t(G5) %*% D4)

# Transpose to matrix to match D3
#r2 <- t(r)
r2 <- t(r[-nrow(r), ])  # Exclude the last row (unit sum constraint)

# reshape back to original dimensions
r_reshape <- array_reshape(r2, c(dim(D1)[1], dim(D1)[2], dim(r2)[2]), order = "C")

# Convert back to raster while preserving CRS
r_raster <- rast(r_reshape, ext = ext(D1), crs = crs_D1)

# Compute mixture fractions
u <- solve(t(G5) %*% G5) %*% t(G5) %*% D4

# Transpose the matrix to match D3 format
u2 <- t(u)

# Reshape back to original dimensions
u_reshape <- array_reshape(u2, c(dim(D1)[1], dim(D1)[2], dim(u2)[2]), order = "C")

# Convert back to raster while preserving CRS
u_raster <- rast(u_reshape, ext = ext(D1), crs = crs_D1)

# layer 1 = unvegetation, layer 2 = evergreen, layer 3 = annual, and layer 4 = deciduous

names(r_raster) <- names(D1)

# write to file
terra::writeRaster(
  r_raster,
  '06_DyeCreek_tem_unmixed/DyeCreek_LS_1982_2025_TMR_time_series',
  overwrite = TRUE,
  datatype = "FLT4S",
  filetype = 'ENVI'
)

# write temporally unmixed fractions to file
terra::writeRaster(
  u_raster,
  '1_Data/Landsat/06_DyeCreek_tem_unmixed/full/DyeCreek_LS_1982_2025_TM_results',
  overwrite = TRUE,
  datatype = "FLT4S",
  filetype = 'ENVI'
)

# RMSE
# Square residuals, mean across bands, take square root
rmse_map <- sqrt(apply(r_reshape^2, c(1, 2), mean, na.rm = TRUE))

# Convert to raster
rmse_raster <- rast(rmse_map, ext = ext(D1), crs = crs_D1)
plot(rmse_raster)

# Compute RMSE from residual matrix
rmse_all <- sqrt(mean(r2^2, na.rm = TRUE))
rmse_all

# unmix results
u_raster <- rast('1_Data/Landsat/06_DyeCreek_tem_unmixed/DyeCreek_LS_1982_2025_TM_results')

##%######################################################%##
#                                                          #
####                  TMR PCA results                   ####
#                                                          #
##%######################################################%##

tmr_pca <- rast('1_Data/Landsat/06_DyeCreek_tem_unmixed/full/DyeCreek_LS_1982_2025_TMR_PCA')
tmr_pca <- tmr_pca[[1:3]]
names(tmr_pca) <- c('TMR1', 'TMR2', 'TMR3')

# Band 3 is best aligned with trends of loss vs. growth, low values correspond with declines in vegetation while higher vaues correspond to increases
plot(tmr_pca[[3]], range = c(-.3, .3))

# plot eigenvectors:
eigs <- read.table('1_Data/Landsat/06_DyeCreek_tem_unmixed/full/DyeCreek_LS_1982_2025_TMR_PCA_eigs.txt', fill = T, header = F) %>% tibble()


##%######################################################%##
#                                                          #
####           Temporal trends of the TMR PCS           ####
#                                                          #
##%######################################################%##

# Load libraries
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(shiny)

# Set file paths (adjust as needed)
veg_path <- "1_Data/Landsat/05_DyeCreek_temporal_stack/DyeCreek_LS_1982_2025_veg_stack"

# Load raster stacks
veg_stack <- rast(veg_path)

# Create sample points within study area
set.seed(2024)
sample_points <- spatSample(veg_stack[[1]], size = 1000, method = "random", as.points = TRUE)

# Extract veg time series and TMR values at points
veg_vals <- terra::extract(veg_stack, sample_points)
tmr_vals <- terra::extract(tmr_pca, sample_points)

# Combine into data frame
veg_tmr_sample_df <- cbind(veg_vals, tmr_vals[, 2:4]) # remove ID column in tmr_vals

# Rename and clean
colnames(veg_tmr_sample_df)[(ncol(veg_tmr_sample_df)-2):ncol(veg_tmr_sample_df)] <- c("TMR1", "TMR2", "TMR3")
veg_tmr_sample_df$ID <- NULL

# Add pixel ID
veg_tmr_sample_df$pixel_id <- seq_len(nrow(veg_tmr_sample_df))

# Pivot to long format
veg_tmr_sample_df_long <- veg_tmr_sample_df |>
  pivot_longer(
    cols = matches("^[0-9]{4}"),  # columns that start with a 4-digit year
    names_to = "date",
    values_to = "veg_score"
  ) |>
  mutate(
    date = as.numeric(date),
    TMR1_cat = case_when(
      TMR1 >= quantile(TMR1, 0.9, na.rm = TRUE) ~ "TMR1_high",
      TMR1 <= quantile(TMR1, 0.1, na.rm = TRUE) ~ "TMR1_low",
      TRUE ~ "TMR1_mid"
    ),
    TMR2_cat = case_when(
      TMR2 >= quantile(TMR2, 0.9, na.rm = TRUE) ~ "TMR2_high",
      TMR2 <= quantile(TMR2, 0.1, na.rm = TRUE) ~ "TMR2_low",
      TRUE ~ "TMR2_mid"
    ),
    TMR3_cat = case_when(
      TMR3 >= quantile(TMR3, 0.9, na.rm = TRUE) ~ "TMR3_high",
      TMR3 <= quantile(TMR3, 0.1, na.rm = TRUE) ~ "TMR3_low",
      TRUE ~ "TMR3_mid"
    )
  )


# Compute mean and sd of veg_score by year and TMR category
avg_trends <- bind_rows(
  veg_tmr_sample_df_long |> group_by(date, TMR_category = TMR1_cat) |> summarise(mean_v_score = mean(veg_score, na.rm = TRUE)),
  veg_tmr_sample_df_long |> group_by(date, TMR_category = TMR2_cat) |> summarise(mean_v_score = mean(veg_score, na.rm = TRUE)),
  veg_tmr_sample_df_long |> group_by(date, TMR_category = TMR3_cat) |> summarise(mean_v_score = mean(veg_score, na.rm = TRUE))
)


ui <- fluidPage(
  titlePanel("Vegetation Trends by TMR Category"),
  sidebarLayout(
    sidebarPanel(
      selectInput("tmr_component", "Choose TMR Component:",
                  choices = c("TMR1", "TMR2", "TMR3"),
                  selected = "TMR3"),
      checkboxGroupInput("tmr_level", "TMR Level:",
                         choices = c("High" = "high", "Low" = "low"),
                         selected = c("high", "low")),
      sliderInput("y_range", "Y-axis Range:", min = -0.2, max = 1, value = c(-0.1, 0.6), step = 0.05),
      checkboxInput("add_vlines", "Show Event Years", value = TRUE)
    ),
    mainPanel(
      plotOutput("trendPlot", height = "600px")
    )
  )
)

server <- function(input, output) {
  output$trendPlot <- renderPlot({
    selected_categories <- paste0(input$tmr_component, "_", input$tmr_level)

    plot_data <- avg_trends |>
      filter(TMR_category %in% selected_categories)

    p <- ggplot(plot_data, aes(x = date, y = mean_v_score, color = TMR_category)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c(
        "TMR1_high" = "red", "TMR1_low" = "blue",
        "TMR2_high" = "darkgreen", "TMR2_low" = "lightgreen",
        "TMR3_high" = "purple", "TMR3_low" = "orange"
      )) +
      theme_minimal(base_size = 16) +
      labs(
        title = paste("Vegetation Trends by", input$tmr_component),
        x = "Year", y = "Mean Vegetation Fraction",
        color = "Category"
      ) +
      ylim(input$y_range)

    if (input$add_vlines) {
      event_years <- c(2015, 2017, 2019, 2021, 2023, 2025)
      p <- p + geom_vline(xintercept = event_years, linetype = "dashed", color = "gray70")
    }

    p
  })
}

shinyApp(ui = ui, server = server)

# TMR 1 trend
tmr1_trend <- veg_tmr_sample_df_long |>
  mutate(
    TMR1_group = case_when(
      TMR1 <= quantile(TMR1, 0.05, na.rm = TRUE) ~ "Q1 (lowest)",
      TMR1 > quantile(TMR1, 0.05, na.rm = TRUE) & TMR1 <= quantile(TMR1, 0.35, na.rm = TRUE) ~ "Q2",
      TMR1 > quantile(TMR1, 0.35, na.rm = TRUE) & TMR1 <= quantile(TMR1, 0.65, na.rm = TRUE) ~ "Q3",
      TMR1 > quantile(TMR1, 0.65, na.rm = TRUE) & TMR1 <= quantile(TMR1, 0.95, na.rm = TRUE) ~ "Q4",
      TMR1 > quantile(TMR3, 0.95, na.rm = TRUE) ~ "Q5 (highest)",
      TRUE ~ NA_character_
    ),
    TMR1_group = factor(TMR1_group, levels = c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)"))
  )

quintile_trends <- tmr1_trend |>
  group_by(date, TMR1_group) |>
  summarise(
    mean_v_score = mean(veg_score, na.rm = TRUE),
    .groups = "drop"
  )

p_quintiles <- ggplot(
  quintile_trends %>% na.omit() %>% filter(TMR1_group %in% c('Q1 (lowest)', 'Q5 (highest)')),
  aes(x = date, y = mean_v_score, color = TMR1_group)
) +
  geom_line(size = 1.5) +
  scale_color_manual(values = c(
    "Q1 (lowest)" = "#0D0887FF",
    # "Q2" = "#7301A8FF",
    # "Q3" = "#BD3786FF",
    # "Q4" = "#ED7953FF",
    "Q5 (highest)" = "#FDC926FF"
  )) +
  labs(
    title = "Vegetation Trends by TMR1 Quintile",
    x = "Year",
    y = "Mean Vegetation Fraction",
    color = "TMR3 Quintile"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    legend.position = 'none',
    legend.background = element_rect(fill = "white", color = "black"),
    legend.direction = 'horizontal'
  ) +
  scale_x_continuous(breaks = c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2025)) 

print(p_quintiles)

# Prepare raster
r <- tmr_pca[[1]]
raster_df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
colnames(raster_df) <- c("x", "y", "value")

# Define uneven breaks once, based on TMR3 values (raster or sample points)
uneven_breaks <- quantile(tmr_pca[[1]][], probs = c(0, 0.05, 0.35, 0.65, 0.95, 1), na.rm = TRUE)
uneven_labels <- c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)")

# For raster data
raster_df <- as.data.frame(tmr_pca[[1]], xy = TRUE, na.rm = TRUE)
colnames(raster_df) <- c("x", "y", "value")

raster_df <- raster_df %>%
  mutate(
    quintile = cut(value,
                   breaks = uneven_breaks,
                   include.lowest = TRUE,
                   labels = uneven_labels)
  )


# Assign reversed plasma colors
plasma_colors <- c("#0D0887FF","#7301A8FF", "#BD3786FF",  "#ED7953FF","#FDC926FF")
names(plasma_colors) <- uneven_labels

legend_title <- "PC1 Quintile"
bbox <- st_bbox(shape)

# Plot
p <- ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = quintile)) +
  scale_fill_manual(values = plasma_colors, name = legend_title) +
  geom_sf(data = st_as_sf(shape), fill = NA, color = "black", lwd = 1) +
  annotation_scale(location = "bl", width_hint = 0.2, height = unit(0.6, "cm"), text_cex = 1.5) +
  annotation_north_arrow(location = "br", which_north = "true", height = unit(2.5, "cm"), width = unit(1.75, "cm")) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE) +
  theme_minimal() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.position = c(.8, .6),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24)
  )

p


library(patchwork)

# Combine the plots side-by-side or top-bottom
combined_plot <- (p_quintiles + p) +
  plot_layout(widths = c(0.4, 0.6)) 

print(combined_plot)

# save
ggsave(
  "H:/My Drive/16_Oak_Woodland/3_Figures/Manuscript_figures/Randall_S2_TMR3_residuals_quintiles.png",
  plot = combined_plot,
  width = 30,
  height = 12,
  dpi = 600
)


