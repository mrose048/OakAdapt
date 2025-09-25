#Load libraries
library(caTools)
library(reticulate)
require(terra)
require(dplyr)
require(sf)

setwd('H:/My Drive/16_Oak_Woodland/1_Data/')

# original crs
og_crs <- terra::rast('Landsat/03_Randall_mosaics/LC08_L2SP_041036_20130427_20200912_02_T1_Scale') %>% crs()

# Randall preserve shapefile
shape <- terra::vect('00_Shapefiles/TNC_Data/Protected_by_TNC_April_2022_dis/Protected_by_TNC_April_2022_dis.shp') %>%
  project(og_crs)



##%######################################################%##
#                                                          #
####         Landsat 8 & 9 - Spectral Unmixing          ####
#                                                          #
##%######################################################%##

#Set filename for data. endmembers, and bad bands

endmember_file = '00_Endmembers/EM_MR_L8_txtNH.txt'
badbands_file = '00_Badbands/L8bbl.csv'

#Set unit sum constraint weight
w = 1

#Multiply EMs by a scale factor?
rescale = 1/10000

files <- list.files('03_Randall_mosaics/', pattern = 'T1_Scale$', full.names = TRUE)

# select only the LC files
files <- files[grep("LC", files)]

outfiles = paste('04_Randall_unmixed/', basename(files),'_Unmixed',sep="")

#Read bad bands list
bands = read.csv(badbands_file, header = F, sep = ",", stringsAsFactors = F)

#Read endmember file
G = read.csv(endmember_file, header = F, sep = "")

for(i in 1:length(files)) {
  
  print(i)
  
 # if(!file.exists(outfiles[i])) {
    
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
    
    # crop and mask to the extent of Randall preserve boundaries
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
    
#  }
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

files <- list.files('03_Randall_mosaics/', pattern = 'T1_Scale$', full.names = TRUE)

# select only files with LE07 pattern
files <- files[grepl('LE07', files)]

outfiles = paste('04_Randall_unmixed/', basename(files),'_Unmixed',sep="")

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
    
    # crop and mask to the extent of Randall preserve boundaries
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

endmember_file = '00_Endmembers/EM_MR_L5_txtNH.txt'
badbands_file = '00_Badbands/L5bbl.csv'

#Set unit sum constraint weight
w = 1

#Multiply EMs by a scale factor?
rescale = 1/10000

files <- list.files('03_Randall_mosaics/', pattern = 'T1_Scale$', full.names = TRUE)

# select only files with LT pattern
files <- files[grepl('LT05', files)]

outfiles = paste('04_Randall_unmixed/', basename(files),'_Unmixed',sep="")

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
    
    # crop and mask to the extent of Randall preserve boundaries
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

files <- list.files('03_Randall_mosaics/', pattern = 'T1_Scale$', full.names = TRUE)

# select only files with LT pattern
files <- files[grepl('LT04', files)]

outfiles = paste('04_Randall_unmixed/', basename(files),'_Unmixed',sep="")

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
    
    # crop and mask to the extent of Randall preserve boundaries
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


# rename files in '04_Randall_unmixed/' so that the first string of characters is added to the end of the file

rename_files <- function(directory) {
  # Get list of files in the directory
  files <- list.files(directory, full.names = TRUE)
  
  for (file in files) {
    # Extract filename and extension
    filename <- basename(file)
    extension <- tools::file_ext(filename)
    filename_no_ext <- tools::file_path_sans_ext(filename)
    
    # Split by underscore
    parts <- unlist(strsplit(filename_no_ext, "_"))
    
    # Check for expected format (at least 8 parts)
    if (length(parts) < 8) {
      cat("Skipping malformed file:", filename, "\n")
      next
    }
    
    # Extract components
    sensor_mission <- parts[1]  # e.g., LC08
    product <- parts[2]         # e.g., L2SP
    path_row <- parts[3]        # e.g., 041036
    date1 <- parts[4]           # e.g., 20130630
    date2 <- parts[5]           # e.g., 20200912
    version <- parts[6]         # e.g., 02
    tier <- parts[7]            # e.g., T1
    scale_unmixed <- paste(parts[8:length(parts)], collapse = "_")  # remainder
    
    # Create new filename with extension
    new_filename <- paste(product, date1, date2, version, tier, sensor_mission, path_row, scale_unmixed, sep = "_")
    new_filename <- paste0(new_filename, ".", extension)
    new_filepath <- file.path(directory, new_filename)
    
    # If file with new name already exists, overwrite it
    if (file.exists(new_filepath)) {
      file.remove(new_filepath)
    }
    
    # Rename file
    success <- file.rename(file, new_filepath)
    if (success) {
      cat("Renamed:", filename, "->", new_filename, "\n")
    } else {
      cat("Failed to rename:", filename, "\n")
    }
  }
}



# Set directory path where the files are stored
directory_path <- "04_Randall_unmixed/"  # Change this to your actual directory

# Run the function
rename_files(directory_path)

##%######################################################%##
#                                                          #
####                   Layer stacking                   ####
#                                                          #
##%######################################################%##

library(lubridate)

# Clean and parse .jpg filenames
ls_a <- list.files(
  "H:/My Drive/16_Oak_Woodland/1_Data/Landsat/TNC_Landsat_Paths/Randall/A_new/",
  pattern = "\\.jpg$", full.names = FALSE
) %>%
  str_remove("\\.jpg$") %>%
  str_replace_all("L1TP|L1GS", "L2SP") %>%
  discard(~ grepl("T2", .) | grepl(" \\(1\\)$", .))

# Split by underscore and extract sensor (1st) and date (4th)
parts <- str_split_fixed(ls_a, "_", n = 7)
clear_ids <- paste0(parts[, 1], "_", parts[, 4])

files <- list.files("1_Data/Landsat/04_Randall_unmixed/", pattern = "Unmixed$", full.names = TRUE)

file_parts <- str_split_fixed(basename(files), "_", n = 9)
file_ids <- paste0(file_parts[, 6], "_", file_parts[, 2])  # LC08_20130427

selected_files <- files[file_ids %in% clear_ids]

# Get the base filenames
selected_fnames <- basename(selected_files)

# Split and extract: date = 2nd, sensor = 6th element
file_parts <- str_split_fixed(selected_fnames, "_", n = 9)
sensor_dates <- paste0(file_parts[, 6], "_", file_parts[, 2])

# write to file
writeLines(
  sensor_dates,
  "1_Data/Landsat/Randall_Landsat_dates.txt"
)

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

frac_years <- calc_frac_year(selected_files)


r <- lapply(selected_files, rast)

# stack the 2nd layer in each element of r (vegetation fraction)
r2 <- lapply(r, function(x) x[[2]])

# stack the rasters
veg <- rast(r2)
names(veg) <- frac_years

# mask of NA values in veg
mask <- any(is.na(veg))

writeRaster(
  mask,
  '05_Randall_temporal_stack/Randall_Landsat_mask',
  overwrite = TRUE,
  filetype = 'ENVI'
)

# mask out NA values
masked_veg <- mask(veg, mask, maskvalues = 1)


# Save as ENVI file
writeRaster(
  masked_veg,
  '05_Randall_temporal_stack/Randall_Landsat_1982_2025',
  datatype = "FLT4S",
  overwrite = TRUE,
  filetype = 'ENVI'
)


# Change wavelength names in ENVI hdr file
hdr_file <- '05_Randall_temporal_stack/Randall_Landsat_1982_2025.hdr'

# Read the header file
hdr <- readLines(hdr_file)

# Add a line for wavelengths = {
hdr <- c(hdr[1:grep('data ignore value = nan', hdr)], 'wavelength = {')

# repeat lines 14 through 482
for (i in 14:482) {
  hdr <- c(hdr, hdr[i])
}

# write hdr file
writeLines(hdr, '05_Randall_temporal_stack/Randall_Landsat_1982_2025.hdr')
# have to make a copy of the ENVI file in folder to match the new header


##%######################################################%##
#                                                          #
####      Calendar plot for Landsat data - Randall      ####
#                                                          #
##%######################################################%##

# Extract relevant info from filenames
df <- data.frame(file = selected_files) %>%
  mutate(
    base_name = basename(file),
    parts = str_split(base_name, "_", simplify = TRUE),
    acq_date_str = parts[, 2],
    sensor = parts[, 6],
    acq_date = ymd(acq_date_str),
    year = year(acq_date),
    doy = yday(acq_date)
  )

# Order sensors for plotting
df$sensor <- factor(df$sensor, levels = c("LT04", "LT05", "LE07", "LC08", "LC09"))

# Preview
print(df)

# Plot

ggplot(df, aes(x = doy, y = as.numeric(year), color = sensor)) +
  geom_point(size = 5) +  # Larger points with white outline
  theme_minimal(base_size = 16) +  # Base size for better readability
  labs(
    x = 'Day of Year',
    y = 'Year',
    color = 'Satellite'
  ) +
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#D95F02", "#E7298A", "#66A61E")) +
  theme(
    # --- Background and grid improvements ---
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_line(color = "gray50", linetype = "dashed"),
    panel.grid.minor = element_blank(),  # Removing minor grids for clarity
    
    # --- Axis and title adjustments ---
    axis.text = element_text(size = 28, color = "white"),
    axis.title = element_text(size = 30, face = "bold", color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # plot.title = element_text(hjust = 0.5, size = 22, face = "bold", color = "white"),
    
    # --- Legend modifications ---
    legend.position = 'bottom',
    legend.title = element_text(size = 30, face = "bold", color = "white"),
    legend.text = element_text(size = 28, color = "white"),
    legend.key = element_rect(fill = 'black', color = NA),
    legend.background = element_rect(fill = 'black', color = NA)
  ) +
  scale_y_continuous(breaks=seq(1980, 2025, 5))

# save image
ggsave(
  "H:/My Drive/16_Oak_Woodland/3_Figures/Manuscript_figures/Randall_Landsat_1982_2025_calendar_plot.png",
  width = 10,
  height = 10,
  dpi = 600
)


##%######################################################%##
#                                                          #
####              Temporal Mixture Models               ####
#                                                          #
##%######################################################%##

# Temporal endmembers
endmember_file  = 'Landsat/05_Randall_temporal_stack/Randall_Landsat_PCA_tems.txt'

#Set unit sum constraint weight
w = 1

#Multiply EMs by a scale factor?
# rescale = 1/10000

# Create a dataframe with 456 columns with 1
bands = data.frame(matrix(1, nrow = 1, ncol = 456))

#Read endmember file
G = read.csv(endmember_file, header = F, sep = "")

# Read spatial data file
D1 <- rast('Landsat/05_Randall_temporal_stack/Randall_Landsat_1982_2025')  # Read as a raster
crs_D1 <- crs(D1)  # Store the original CRS

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

# RMSE
# Square residuals, mean across bands, take square root
rmse_map <- sqrt(apply(r_reshape^2, c(1, 2), mean, na.rm = TRUE))

# Convert to raster
rmse_raster <- rast(rmse_map, ext = ext(D1), crs = crs_D1)
# plot(rmse_raster)

# Compute RMSE from residual matrix
rmse_all <- sqrt(mean(r2^2, na.rm = TRUE))
rmse_all

# Compute mixture fractions
u <- solve(t(G5) %*% G5) %*% t(G5) %*% D4

# Transpose the matrix to match D3 format
u2 <- t(u)

# Reshape back to original dimensions
u_reshape <- array_reshape(u2, c(dim(D1)[1], dim(D1)[2], dim(u2)[2]), order = "C")

# Convert back to raster while preserving CRS
u_raster <- rast(u_reshape, ext = ext(D1), crs = crs_D1)

names(r_raster) <- names(D1)

# write to file
terra::writeRaster(
  r_raster,
  '06_Randall_tem_unmixed/Randall_Landsat_1982_2025_TMR',
  overwrite = TRUE,
  datatype = "FLT4S",
  filetype = 'ENVI'
)

# write temporally unmixed fractions to file
terra::writeRaster(
  u_raster,
  '06_Randall_tem_unmixed/Randall_Landsat_1982_2025_TEM',
  overwrite = TRUE,
  datatype = "FLT4S",
  filetype = 'ENVI'
)

# Residual mixture fraction time series
tmr <- rast('06_Randall_tem_unmixed/Randall_Landsat_1982_2025_TMR')

# Temporal endmember fractions
tem <- rast('06_Randall_tem_unmixed/Randall_Landsat_1982_2025_TEM')

# Mask tmr based on fraction of non-annual veg (lyr2)
n_annual <- tem$lyr.2 + tem$lyr.4
mask <- n_annual > .15

# Mask tmr based on >15% non annual vegetation cover
tmr <- tmr %>%
  mask(mask, maskvalues = 0, updatevalue = NA)

# write raster
terra::writeRaster(
  tmr,
  '06_Randall_tem_unmixed/Randall_Landsat_1982_2025_TMR_oak_mask',
  overwrite = TRUE,
  datatype = "FLT4S",
  filetype = 'ENVI'
)



# Principal component of the residual mixture fractions
tmr_pca <- rast('06_Randall_tem_unmixed/Randall_Landsat_1982_2024_TMR_PCA') 


# Create ggplot of tmr_pca[[1]]
tmr_pca_df <- as.data.frame(tmr_pca[[1:3]], xy = TRUE)

# rename
names(tmr_pca_df) <- c('x', 'y', 'TMR1', 'TMR2', 'TMR3')

ggplot() +
  geom_raster(data = tmr_pca_df, aes(x = x, y = y, fill = TMR1)) +
  tidyterra::geom_spatvector(
    data = shape,
    fill = 'transparent',
    color = 'black',
    lwd = 0.5
  ) +
  scale_fill_distiller(
    palette = "BrBG",
    direction = -1,
    limits = c(-.5, .5),  # Set fill limits
    oob = scales::squish  # Prevents errors by squishing out-of-range values
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 1.5,
    unit_category = 'metric'
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )


# save
ggsave("H:/My Drive/16_Oak_Woodland/3_Figures/Randall_Landsat_1982_2024_TMR1.png", width = 10, height = 10, dpi = 300)

# TMR2
ggplot() +
  geom_raster(data = tmr_pca_df, aes(x = x, y = y, fill = TMR2)) +
  tidyterra::geom_spatvector(
    data = shape,
    fill = 'transparent',
    color = 'black',
    lwd = 0.5
  ) +
  scale_fill_distiller(
    palette = "BrBG",
    direction = -1,
    limits = c(-.5, .5),  # Set fill limits
    oob = scales::squish  # Prevents errors by squishing out-of-range values
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 1.5,
    unit_category = 'metric'
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

# save
ggsave("H:/My Drive/16_Oak_Woodland/3_Figures/Randall_Landsat_1982_2024_TMR2.png", width = 10, height = 10, dpi = 300)


# TMR 3
ggplot() +
  geom_raster(data = tmr_pca_df, aes(x = x, y = y, fill = TMR3)) +
  tidyterra::geom_spatvector(
    data = shape,
    fill = 'transparent',
    color = 'black',
    lwd = 0.5
  ) +
  scale_fill_distiller(
    palette = "BrBG",
    direction = -1,
    limits = c(-.5, .5),  # Set fill limits
    oob = scales::squish  # Prevents errors by squishing out-of-range values
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 1.5,
    unit_category = 'metric'
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

# save
ggsave("H:/My Drive/16_Oak_Woodland/3_Figures/Randall_Landsat_1982_2024_TMR3.png", width = 10, height = 10, dpi = 300)


# Extract dates from raster (if available, otherwise define them manually)
time_steps <- seq(from = as.Date("2013-01-01"), to = as.Date("2020-12-31"), length.out = nlyr(tmr_time_series))

# Compute mean temporal mixture residuals per time step
tmr_means <- as.data.frame(t(apply(as.matrix(tmr_time_series), 2, mean, na.rm = TRUE)))
colnames(tmr_means) <- "TMR"
tmr_means$Date <- time_steps  # Assign corresponding time labels

# Convert to long format for ggplot
tmr_means_long <- tmr_means %>%
  tidyr::gather(key = "Component", value = "Value", -Date)


# project to lat long
tmr_pca <- tmr_pca %>% project('+proj=longlat +datum=NAD83')

# oak plots
oak_plots <- readxl::read_xlsx(
  '00_Shapefiles/TNC_Data/Oak Plot Tree Master List (2).xlsx'
) %>% tibble()
blue_oak <- oak_plots %>%
  filter(Species == 'Blue Oak')

# omit observations where seedling count is NA
blue_oak <- blue_oak %>%
  filter(!is.na(`Seedling Count`)) %>%
  mutate(`Seedling and Sapling Count` = `Seedling Count` + `Sapling Count`) %>%
  # mutate a column to classify seedling and sapling count
  mutate(`Seedling and Sapling Class` = case_when(
    `Seedling and Sapling Count` == 0 ~ 'None',
    `Seedling and Sapling Count` >= 1 & `Seedling and Sapling Count` < 2 ~ 'Low',
    `Seedling and Sapling Count` >= 2 & `Seedling and Sapling Count` < 4 ~ 'Medium',
    `Seedling and Sapling Count` >= 4 ~ 'High'
  ))

# plots as vector
# set levels of Seedling and Sapling Class to None, Low, Medium, High
blue_oak$`Seedling and Sapling Class` <- factor(
  blue_oak$`Seedling and Sapling Class`,
  levels = c('High', 'Medium', 'Low', 'None')
)

# Write as vector
blue_oak_sv <- terra::vect(
  blue_oak %>% dplyr::select(
    'Ranch',
    'Tree ID',
    'longitude',
    'latitude',
    'Species',
    'Seedling and Sapling Count',
    'Seedling and Sapling Class'
  ) %>%
    rename(
      'ID' = 'Tree ID',
      'count' = 'Seedling and Sapling Count',
      'class' = 'Seedling and Sapling Class'
    ) %>% na.omit(),
  geom = c('longitude', 'latitude'),
  crs = '+proj=longlat +datum=NAD83'
)

blue_oak_sv$class <- as.character(blue_oak_sv$class)

# extract blue_oak_sp for the first three principal components of tmr_pca
tmr_pca <- tmr_pca[[1:3]]
names(tmr_pca) <- c('PC1', 'PC2', 'PC3')
tmr_pca_blue_oak <- terra::extract(tmr_pca, blue_oak_sv, bind = T) 

# as data frame
tmr_pca_blue_oak <- as.data.frame(tmr_pca_blue_oak) %>% tibble()

# GGplot showing oak plot class and first three principal components of tmr_pca
ggplot(tmr_pca_blue_oak, aes(x = class, y  = PC1)) +
  geom_boxplot() +
  theme_minimal()

##%######################################################%##
#                                                          #
####              NEW ANALYSIS - MAY 2025               ####
#                                                          #
##%######################################################%##


# Load libraries
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(shiny)

# Set file paths (adjust as needed)
veg_path <- "05_Randall_temporal_stack/Randall_Landsat_1982_2025"
tmr_pca <- "06_Randall_tem_unmixed/Randall_LS_TMR_PCA"

# Load raster stacks
veg_stack <- rast(veg_path)
tmr_pca <- rast(tmr_pca) 
tmr_pca[[1:3]]

# Create sample points within study area
set.seed(2024)

# Extract the third principal component layer
tmr3 <- tmr_pca[[3]]

# Create a mask of cells that are not NA
non_na_cells <- which(!is.na(values(tmr3)))

# Use those indices to sample from only valid cells
# Number of samples you want:
n_samples <- 2000

# Randomly sample indices from valid cells
set.seed(123)  # for reproducibility
sample_indices <- sample(non_na_cells, size = n_samples)

# Convert cell indices to xy coordinates
sample_points <- xyFromCell(tmr3, sample_indices)

# Convert to SpatVector
sample_points_vect <- vect(sample_points, type = "points", crs = crs(tmr3))


# Extract veg time series and TMR values at points
veg_vals <- terra::extract(veg_stack, sample_points_vect)
tmr_vals <- terra::extract(tmr_pca, sample_points_vect)

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
      TMR1 >= quantile(TMR1, 0.95, na.rm = TRUE) ~ "TMR1_high",
      TMR1 <= quantile(TMR1, 0.05, na.rm = TRUE) ~ "TMR1_low",
      TRUE ~ "TMR1_mid"
    ),
    TMR2_cat = case_when(
      TMR2 >= quantile(TMR2, 0.95, na.rm = TRUE) ~ "TMR2_high",
      TMR2 <= quantile(TMR2, 0.05, na.rm = TRUE) ~ "TMR2_low",
      TRUE ~ "TMR2_mid"
    ),
    TMR3_cat = case_when(
      TMR3 >= quantile(TMR3, 0.95, na.rm = TRUE) ~ "TMR3_high",
      TMR3 <= quantile(TMR3, 0.05, na.rm = TRUE) ~ "TMR3_low",
      TRUE ~ "TMR3_mid"
    )
  )


# Compute mean and sd of veg_score by year and TMR category
avg_trends <- bind_rows(
  veg_tmr_sample_df_long |> group_by(date, TMR_category = TMR1_cat) |> summarise(mean_v_score = mean(veg_score, na.rm = TRUE)),
  veg_tmr_sample_df_long |> group_by(date, TMR_category = TMR2_cat) |> summarise(mean_v_score = mean(veg_score, na.rm = TRUE)),
  veg_tmr_sample_df_long |> group_by(date, TMR_category = TMR3_cat) |> summarise(mean_v_score = mean(veg_score, na.rm = TRUE))
)


# ui <- fluidPage(
#   titlePanel("Vegetation Trends by TMR Category"),
#   sidebarLayout(
#     sidebarPanel(
#       selectInput("tmr_component", "Choose TMR Component:",
#                   choices = c("TMR1", "TMR2", "TMR3"),
#                   selected = "TMR3"),
#       checkboxGroupInput("tmr_level", "TMR Level:",
#                          choices = c("High" = "high", "Low" = "low"),
#                          selected = c("high", "low")),
#       sliderInput("y_range", "Y-axis Range:", min = -0.2, max = 1, value = c(-0.1, 0.6), step = 0.05),
#       checkboxInput("add_vlines", "Show Event Years", value = TRUE)
#     ),
#     mainPanel(
#       plotOutput("trendPlot", height = "600px")
#     )
#   )
# )
# 
# server <- function(input, output) {
#   output$trendPlot <- renderPlot({
#     selected_categories <- paste0(input$tmr_component, "_", input$tmr_level)
# 
#     plot_data <- avg_trends |>
#       filter(TMR_category %in% selected_categories)
# 
#     p <- ggplot(plot_data, aes(x = date, y = mean_v_score, color = TMR_category)) +
#       geom_line(linewidth = 1) +
#       scale_color_manual(values = c(
#         "TMR1_high" = "red", "TMR1_low" = "blue",
#         "TMR2_high" = "darkgreen", "TMR2_low" = "lightgreen",
#         "TMR3_high" = "purple", "TMR3_low" = "orange"
#       )) +
#       theme_minimal(base_size = 16) +
#       labs(
#         title = paste("Vegetation Trends by", input$tmr_component),
#         x = "Year", y = "Mean Vegetation Fraction",
#         color = "Category"
#       ) +
#       ylim(input$y_range)
# 
#     if (input$add_vlines) {
#       event_years <- c(1982, 1986, 1990, 1994, 1998, 2002, 2006, 2010, 2014, 2018, 2022, 2025)
#       p <- p + geom_vline(xintercept = event_years, linetype = "dashed", color = "gray70")
#     }
# 
#     p
#   })
# }
# 
# shinyApp(ui = ui, server = server)

# TMR 3 trend
tmr3_trend <- veg_tmr_sample_df_long |>
  mutate(
    TMR3_group = case_when(
      TMR3 <= quantile(TMR3, 0.05, na.rm = TRUE) ~ "Q1 (lowest)",
      TMR3 > quantile(TMR3, 0.05, na.rm = TRUE) & TMR3 <= quantile(TMR3, 0.35, na.rm = TRUE) ~ "Q2",
      TMR3 > quantile(TMR3, 0.35, na.rm = TRUE) & TMR3 <= quantile(TMR3, 0.65, na.rm = TRUE) ~ "Q3",
      TMR3 > quantile(TMR3, 0.65, na.rm = TRUE) & TMR3 <= quantile(TMR3, 0.95, na.rm = TRUE) ~ "Q4",
      TMR3 > quantile(TMR3, 0.95, na.rm = TRUE) ~ "Q5 (highest)",
      TRUE ~ NA_character_
    ),
    TMR3_group = factor(TMR3_group, levels = c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)"))
  )

quintile_trends <- tmr3_trend |>
  group_by(date, TMR3_group) |>
  summarise(
    mean_v_score = mean(veg_score, na.rm = TRUE),
    .groups = "drop"
  )

p_quintiles <- ggplot(
  quintile_trends %>% na.omit() %>% filter(TMR3_group %in% c('Q1 (lowest)', 'Q5 (highest)')),
  aes(x = date, y = mean_v_score, color = TMR3_group)
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
    title = "Vegetation Trends by TMR3 Quintile",
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
  scale_x_continuous(breaks = c(seq(1982, 2025, by = 6)))

print(p_quintiles)

# Prepare raster
r <- tmr_pca[[3]]
raster_df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
colnames(raster_df) <- c("x", "y", "value")

# Define uneven breaks once, based on TMR3 values (raster or sample points)
uneven_breaks <- quantile(tmr_pca[[3]][], probs = c(0, 0.05, 0.35, 0.65, 0.95, 1), na.rm = TRUE)
uneven_labels <- c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)")

# For raster data
raster_df <- as.data.frame(tmr_pca[[3]], xy = TRUE, na.rm = TRUE)
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
names(plasma_colors) <- quintile_labels

legend_title <- "PC3 Quintile"
bbox <- st_bbox(randall)

# Plot
p <- ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = quintile)) +
  scale_fill_manual(values = plasma_colors, name = legend_title) +
  geom_sf(data = st_as_sf(randall), fill = NA, color = "black", lwd = 1) +
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
  "H:/My Drive/16_Oak_Woodland/3_Figures/Manuscript_figures/RandallLS_TMR3_residuals_quintiles.png",
  plot = combined_plot,
  width = 30,
  height = 12,
  dpi = 600
)


library(terra)

# Select the raster layer
r <- tmr_pca[[3]]

# Define custom quantile breaks
uneven_breaks <- quantile(r[], probs = c(0, 0.05, 0.35, 0.65, 0.95, 1), na.rm = TRUE)

# Reclassify raster based on breaks into categories 1–5
r_class <- classify(r, rcl = matrix(c(
  uneven_breaks[1], uneven_breaks[2], 1,
  uneven_breaks[2], uneven_breaks[3], 2,
  uneven_breaks[3], uneven_breaks[4], 3,
  uneven_breaks[4], uneven_breaks[5], 4,
  uneven_breaks[5], uneven_breaks[6], 5
), ncol = 3, byrow = TRUE))

# Assign factor labels
levels(r_class) <- data.frame(
  ID = 1:5,
  Category = c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)")
)

writeRaster(r_class,
            filename = "08_Outputs/Randall_Landsat_change_quintiles.tif",
            filetype = "GTiff",  # GeoTIFF format
            overwrite = TRUE)

writeRaster(r,
            filename = "08_Outputs/Randall_Landsat_change_raw.tif",
            filetype = "GTiff",  # GeoTIFF format
            overwrite = TRUE)
