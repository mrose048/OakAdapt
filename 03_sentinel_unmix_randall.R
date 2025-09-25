##%######################################################%##
#                                                          #
####                 Sentinel Unmixing                  ####
#                                                          #
##%######################################################%##

#Load libraries
library(caTools)
library(reticulate)
require(terra)
require(dplyr)
require(sf)
library(stringr)
require(lubridate)
require(tidyverse)

setwd('H:/My Drive/16_Oak_Woodland/1_Data/')

#### Troubleshooting issue with repeated vegetation fraction across multiple time steps

files <- list.files('Sentinel2_folder/03_Randall_envi', pattern = 'e$', full.names = TRUE)

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

# Check output
print(bad_date_strings)

# Create one regular expression pattern that matches any of the bad dates
pattern <- paste(bad_date_strings, collapse = "|")

# Filter matched files based on the combined pattern
bad_files <- files[grepl(pattern, files)]

# Print matched bad file paths
print(bad_files)

# call in files as rasters
bad_files_r <- lapply(bad_files[6:10], rast)



##%######################################################%##
#                                                          #
####                  Randall Preserve                  ####
#                                                          #
##%######################################################%##


endmember_file = '00_Endmembers/EM_MR_Sentinel.txt'
badbands_file = '00_Badbands/S2bbl.csv'


#Set unit sum constraint weight
w = 1

#Multiply EMs & images by a scale factor?
rescale = 1/10000

files <- list.files('Sentinel2_folder/03_Randall_envi', pattern = 'e$', full.names = TRUE)

outfiles = paste('Sentinel2_folder/04_Randall_unmixed/', basename(files),'_Unmixed',sep="")

#Read bad bands list
bands = read.csv(badbands_file, header = F, sep = ",", stringsAsFactors = F)

#Read endmember file
G = read.csv(endmember_file, header = F, sep = "")

# check on S2_11SLU_20180208T183521_L2A_e, S2_11SLV_20171011T183311_L2A_e, S2_11SLV_20171230T183751_L2A_e, S2_11SLV_20191021T183411_L2A_e

for(i in 71:length(files)) {
  
  print(i)
  
  if(!file.exists(outfiles[i])) {
    
    # Read spatial data file
    D1 <- rast(files[i])  # Read as a raster
    crs_D1 <- crs(D1)  # Store the original CRS
    
    # scale for surface reflectance
    D1 <- D1/10000
    
    # Get wavelengths from the first column of the endmember file
    wavelengths <- G[, 1]
    #  G2 <- G[, 2:ncol(G)]
    
    G2 <- G[, 2:4]
    
    # Ensure all values are numeric
    G2 <- apply(G2, 2, function(x) as.numeric(as.character(x)))
    
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
    # u_raster <- u_raster
    
    # write to file
    terra::writeRaster(
      u_raster,
      outfiles[i],
      overwrite = TRUE,
      datatype = "FLT4S",
      filetype = 'ENVI'
    )
    
    gc()
    
  }
}

##%######################################################%##
#                                                          #
####                   Layer Stacking                   ####
#                                                          #
##%######################################################%##

# read a sentinel image file for crs
og_crs <- terra::rast(
  'H:/My Drive/16_Oak_Woodland/1_Data/Sentinel2_folder/02_Randall_unzip/S2A_MSIL2A_20150922T184436_N0500_R027_T11SLU_20231017T001227.SAFE/GRANULE/L2A_T11SLU_A001311_20150922T184430/IMG_DATA/R20m/T11SLU_20150922T184436_AOT_20m.jp2'
) %>%
  crs()

# randall shapefile
randall <- terra::vect(
  'H:/My Drive/16_Oak_Woodland/1_Data/00_Shapefiles/TNC_Data/Protected_by_TNC_April_2022_dis/Protected_by_TNC_April_2022_dis.shp'
) %>%
  project(og_crs)

# stack vegetation fraction images & rename with fractional year
library(lubridate)

# clear image list
ls_a <- list.files(
  'H:/My Drive/16_Oak_Woodland/1_Data/Sentinel2_folder/calendar_plots/Randall/A',
  pattern = '.jpg',
  full.names = F
)
# remove .jpg from items in ls_a
ls_a <- gsub('.jpg', '', ls_a)
ls_a <- gsub('preview', '', ls_a)
# Remove items containing ' (1)' and (2)
ls_a <- ls_a[!grepl("\\(1\\)|\\(2\\)", ls_a)]
# Just the date (YYYYMMDD)
ls_a <- gsub("\\-", "", ls_a)

# Stack all files in "04_Randall_unmixed/" that were in the original clear folder and save as Envi raster
files <- list.files('1_Data/Sentinel_2/04_Randall_unmixed/', pattern = 'Unmixed$', full.names = TRUE)

# Filter files that contain any of the clear Landsat dates
selected_files <- files[sapply(files, function(file) any(grepl(paste(ls_a, collapse = "|"), file)))]

# scene 11SLV
files_11SLV <- selected_files[grepl("11SLV", selected_files)]
files_11SLU <- selected_files[grepl("11SLU", selected_files)]

# Assuming files_11SLV and files_11SLU are vectors of file paths
# Extract dates from filenames
dates_11SLV <- str_extract(files_11SLV, "\\d{8}")  # Extracts date in format YYYYMMDD
dates_11SLU <- str_extract(files_11SLU, "\\d{8}")

# Match files based on the date
matched_files <- intersect(dates_11SLV, dates_11SLU)  # Get common dates between the two scenes

# write to file:
writeLines(matched_files, '1_Data/Sentinel_2/RandallS2_dates.txt')

#### Troubleshooting issue with repeated vegetation fraction across multiple time steps


# dates with same vegetation fraction after unmixing
# bad_dates_df <- data.table::fread('bad_dates.csv') %>% tibble()

# Full list of bad dates
# bad_dates <- bad_dates_df$date

# Convert fractional years to Dates
# bad_date_objs <- sapply(bad_dates, fraction_to_date)

# Convert Dates to character format matching Sentinel filenames
# bad_date_strings <- format(as.Date(bad_date_objs, origin = "1970-01-01"), "%Y%m%d")

# Check output
# print(bad_date_strings)

# Create one regular expression pattern that matches any of the bad dates
# pattern <- paste(bad_date_strings, collapse = "|")

# Filter matched files based on the combined pattern
#bad_files <- matched_files[grepl(pattern, matched_files)]

# Print matched bad file paths
#print(bad_files)

# call in files from files_11SLV that match
#bad_files_11SLV <- files_11SLV[grepl(paste(bad_date_strings, collapse = "|"), files_11SLV)]
#bad_files_11SLV_r <- lapply(bad_files_11SLV, rast) # issue is happening before before unmixing

#bad_files_11SLU <- files_11SLU[grepl(paste(bad_date_strings, collapse = "|"), files_11SLU)]
#bad_files_11SLU_r <- lapply(bad_files_11SLU, rast) # issue is happening before before unmixing


##%######################################################%##
#                                                          #
####       Mosaic the 2 scenes for matching days        ####
#                                                          #
##%######################################################%##

source('H:/My Drive/16_Oak_Woodland/2_Scripts/00_functions.R')


# Loop over each matched date
for (date in matched_files) {
  print(date)
  clipped_11SLV <- clip_sentinel_image(files_11SLV, randall, date)
  clipped_11SLU <- clip_sentinel_image(files_11SLU, randall, date)
  
  # Mosaic the images
  mosaic_stack <- merge(clipped_11SLV, clipped_11SLU)
  
  # Export each mosaicked image (for each date)
  writeRaster(
    mosaic_stack,
    paste0("Sentinel2_folder/04_Randall_unmixed/S2_unmix", date, "_mosaicked"),
    datatype = "FLT4S",
    overwrite = TRUE,
    filetype = 'ENVI'
  )
}

r <- list.files(
  'Sentinel2_folder/04_Randall_unmixed/',
  pattern = 'mosaicked$',
  full.names = TRUE
)

# Function to derive fractional year from file names
library(lubridate)
calc_frac_year <- function(filenames) {
  # Extract all occurrences of an 8-digit date (YYYYMMDD)
  matches <- regmatches(filenames, regexpr("\\d{8}", filenames))
  
  # Convert matches to date objects
  date_objs <- ymd(matches)
  
  # Identify invalid dates and return NA for them
  invalid_dates <- is.na(date_objs)
  if (any(invalid_dates)) {
    warning(paste("Invalid or missing dates in", sum(invalid_dates), "filenames"))
  }
  
  # Extract year and day of year
  years <- year(date_objs)
  doys <- yday(date_objs)
  
  # Check leap years and compute days in year
  days_in_years <- ifelse(leap_year(years), 366, 365)
  
  # Compute fractional year
  frac_years <- years + (doys - 1) / days_in_years
  
  # Assign NA to invalid entries
  frac_years[invalid_dates] <- NA
  
  return(frac_years)
}

frac_years <- calc_frac_year(r)

# Create a list of raster stacks
r <- lapply(r, rast)

# stack the 2nd layer in each element of r (vegetation fraction)
r2 <- lapply(r, function(x) x[[2]])

# common extent 
common_extent <- ext(r2[[101]])

# raster stack
r2 <- lapply(r2, crop, randall)
r2 <- lapply(r2, mask, randall)
r2 <- lapply(r2, function(r) crop(r, common_extent))
veg <- rast(r2)
names(veg) <- frac_years

# mask of NA values in veg
mask <- any(is.na(veg))

# mask out NA values
masked_veg <- mask(veg, mask, maskvalues = 1)

# Save as ENVI file
writeRaster(
  masked_veg,
  'Sentinel2_folder/05_Randall_temporal_stack/Randall_S2_stack',
  datatype = "FLT4S",
  overwrite = TRUE,
  filetype = 'ENVI'
)


# Change wavelength names in ENVI hdr file
hdr_file <- 'Sentinel2_folder/05_Randall_temporal_stack/Randall_S2_stack.hdr'

# Read the header file
hdr <- readLines(hdr_file)

# Add a line for wavelengths = {
hdr <- c(hdr[1:grep('data ignore value = nan', hdr)], 'wavelength = {')

# repeat lines 14 through 295
for (i in 15:296) {
  hdr <- c(hdr, hdr[i])
}

# write hdr file
writeLines(hdr, 'Sentinel2_folder/05_Randall_temporal_stack/Randall_S2_stack.hdr')


##%######################################################%##
#                                                          #
####         Troubleshooting Randall vegetation         ####
####                fraction time series                ####
#                                                          #
##%######################################################%##
require(tidyverse)

masked_veg <- rast('Sentinel2_folder/05_Randall_temporal_stack/Randall_S2_stack')

randall_blue_oaks <- vect('00_Shapefiles/TNC_Data/oak_plots.shp') %>%
  terra::project(masked_veg)

# extract vegetation time series
oak_ex <- terra::extract(masked_veg, randall_blue_oaks)
oak_ex2 <- tibble(oak_ex) %>%
  pivot_longer(
    cols = -c(ID),
    names_to = "date",
    values_to = "fraction"
  )

ggplot(oak_ex2 %>% filter(ID %in% c(1,100, 200, 300)), aes(x = date, y = fraction, group = ID)) +
  geom_line() +
  geom_point() +
  labs(
    title = "",
    x = "Date",
    y = "Vegetation Fractional Cover"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank()) + 
  facet_wrap(~ID)

# list of dates where vegetation fraction = 0.159176588
bad_dates <- oak_ex2 %>%
  filter(fraction >= 0.159176587 & fraction <= 0.159176589) %>%
  pull(date)

# bad_dates as data frame
bad_dates_df <- data.frame(date = bad_dates)

# save
write.csv(bad_dates_df, 'bad_dates.csv', row.names = FALSE)

##%######################################################%##
#                                                          #
####                   Calendar plot                    ####
#                                                          #
##%######################################################%##

# Create a data frame from dates
df <- tibble(
  acq_date_str = matched_files,
  acq_date = ymd(acq_date_str),
  year = year(acq_date),
  doy = yday(acq_date),
  sensor = "Sentinel-2"  # add static sensor label if needed
)

# Plot
ggplot(df, aes(x = doy, y = year, color = sensor)) +
  geom_point(size = 5) +
  theme_minimal(base_size = 16) +
  labs(
    x = 'Day of Year',
    y = 'Year',
    color = 'Satellite'
  ) +
  scale_color_manual(values = c("Sentinel-2" = "#1B9E77")) +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_line(color = "gray50", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 28, color = "white"),
    axis.title = element_text(size = 30, face = "bold", color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = 'bottom',
    legend.title = element_text(size = 30, face = "bold", color = "white"),
    legend.text = element_text(size = 28, color = "white"),
    legend.key = element_rect(fill = 'black', color = NA),
    legend.background = element_rect(fill = 'black', color = NA)
  ) +
  scale_y_continuous(breaks = seq(min(df$year), max(df$year)))

# save image
ggsave(
  "H:/My Drive/16_Oak_Woodland/3_Figures/Manuscript_figures/Randall_Sentinel2_calendar_plot.png",
  width = 10,
  height = 10,
  dpi = 600
)


##%######################################################%##
#                                                          #
####                 Temporal Unmixing                  ####
#                                                          #
##%######################################################%##

# RAN IN PYTHON BECAUSE R DOES NOT HAVE ENOUGH MEMORY FOR SUCH A LARGE STUDY AREA, SO MANY TIME STEPS AT 10M RESOLUTION

# Temporal endmembers
# endmember_file  = 'Sentinel2_folder/05_Randall_temporal_stack/Randall_S2_stack_PCA_TEMs.txt' # the file you created from ENVI PCA
# 
# #Set unit sum constraint weight
# w = 1
# 
# #Multiply EMs by a scale factor?
# # rescale = 1/10000
# 
# # Create a dataframe with 456 columns with 1
# bands = data.frame(matrix(1, nrow = 1, ncol = 281))
# 
# #Read endmember file
# G = read.csv(endmember_file, header = F, sep = "")
# 
# # Read spatial data file
# D1 <- rast('Sentinel2_folder/05_Randall_temporal_stack/Randall_S2_stack')  # Read vegetation fraction stack as a raster
# crs_D1 <- crs(D1)  # Store the original CRS
# 
# # Get wavelengths from the first column of the endmember file
# wavelengths <- G[, 1]
# G2 <- G[, 2:ncol(G)]
# 
# # Subset using bad bands
# D2 <- as.array(D1)[, , bands == 1]  # Extract valid bands
# gc()
# G3 <- G2[bands == 1, ]
# 
# # Reshape arrays
# D3 <- array_reshape(D2, c(dim(D2)[1] * dim(D2)[2], dim(D2)[3]), order = "C")
# rm(D2)
# gc()
# 
# # Add unit sum constraint
# D4 <- t(as.matrix(cbind(D3, rep(1, nrow(
#   D3
# )))))  # Transpose to add unit sum constraint
# G4 <- as.matrix(rbind(G3, rep(1, ncol(G3))))
# 
# # Compute residual of time series - you will transform this with PCA and evaluate the first three PC
# r <- D4 - (G4 %*% solve(t(G4) %*% G4) %*% t(G4) %*% D4)
# 
# # Transpose to matrix to match D3
# #r2 <- t(r)
# r2 <- t(r[-nrow(r), ])  # Exclude the last row (unit sum constraint)
# 
# # reshape back to original dimensions
# r_reshape <- array_reshape(r2, c(dim(D1)[1], dim(D1)[2], dim(r2)[2]), order = "C")
# 
# # Convert back to raster while preserving CRS
# r_raster <- rast(r_reshape, ext = ext(D1), crs = crs_D1)
# 
# # Compute mixture fractions
# u <- solve(t(G5) %*% G5) %*% t(G5) %*% D4
# 
# # Transpose the matrix to match D3 format
# u2 <- t(u)
# 
# # Reshape back to original dimensions
# u_reshape <- array_reshape(u2, c(dim(D1)[1], dim(D1)[2], dim(u2)[2]), order = "C")
# 
# # Convert back to raster while preserving CRS
# u_raster <- rast(u_reshape, ext = ext(D1), crs = crs_D1)
# 
# names(r_raster) <- names(D1)
# 
# # write to file
# terra::writeRaster(
#   r_raster,
#   'Sentinel2_folder/06_Randall_tem_unmixed/Randall_S2_TMR_time_series', # feel free to rename - this raster contains the temporal mixture model residuals
#   overwrite = TRUE,
#   datatype = "FLT4S",
#   filetype = 'ENVI'
# )
# 
# # write temporally unmixed fractions to file
# terra::writeRaster(
#   u_raster,
#   'Sentinel2_folder/06_Randall_tem_unmixed/Randall_S2_TEMs', # temporally unmixed fractions (should be 4 rasters)
#   overwrite = TRUE,
#   datatype = "FLT4S",
#   filetype = 'ENVI'
# )

library(terra)

# Load residual raster
res_rast <- rast("Sentinel_2/05_Randall_temporal_stack/Randall_S2_stack_MixtureResidual.bsq")

# Pixel-wise RMSE across layers (bands)
rmse_rast <- sqrt(app(res_rast, fun = function(x) mean(x^2, na.rm=TRUE)))

# Overall RMSE
rmse_all <- sqrt(mean(values(res_rast)^2, na.rm=TRUE))

print(rmse_all)
plot(rmse_rast)

tem <- rast('Sentinel_2/05_Randall_temporal_stack/Randall_S2_stack_Unmixed.bsq')

# Function to compute % of pixels in 0-1 for a single layer
percent_fun <- function(x) {
  n_in_range <- sum(x >= 0 & x <= 1, na.rm = TRUE)
  n_total <- sum(!is.na(x))
  (n_in_range / n_total) * 100
}

# Apply to all layers
percent_0_1 <- global(tem, fun = percent_fun)

# Add layer names
percent_0_1_df <- data.frame(layer = names(r), percent_in_0_1 = percent_0_1[,1])
percent_0_1_df


# Explore TMM results
tems <- data.table::fread('Sentinel_2/05_Randall_temporal_stack/Randall_S2_stack_PCA_TEMs.txt') %>% tibble()

# rename columns
colnames(tems) <- c('Year', 'Unvegetated', 'Evergreen', 'Deciduous', 'Annual')

# plot
plot_df <- tems %>%
  dplyr::select(-Unvegetated) %>%
  pivot_longer(-Year, names_to = 'Endmember', values_to = 'Score')

plot <- 
  ggplot(plot_df, aes(x = Year, y = Score, color = Endmember)) +
  geom_line(lwd = .9) +
  theme_minimal() +
  theme(
    text = element_text(size = 35, color = "white"),
    # White text
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "white"
    ),
    axis.text.y = element_text(color = "white"),
    plot.title = element_text(color = "white"),
    panel.background = element_rect(fill = "gray25", color = NA),
    # Black panel background
    plot.background = element_rect(fill = "gray25", color = NA),
    # Black plot background
    panel.grid.major = element_line(color = "gray50"),
    # Subtle grid lines
    panel.grid.minor = element_line(color = "gray40"),
    legend.position = 'bottom'
  ) +
  labs(title = '', x = 'Year', y = 'Fv') +
  scale_color_manual(values = c('red', 'cyan', 'green'))   # Adjust 'black' to 'white' for visibility

print(plot)

# unmixed rasters

unmixed_r <- rast('Sentinel2_folder/05_Randall_temporal_stack/Randall_S2_stack_Unmixed.bsq')
names(unmixed_r) <- c('Unvegetated', 'Deciduous', 'Evergreen', 'Annual')
# plot
plot(unmixed_r)


# what proportion of each vegetation type do we see for the oak plots?
randall_blue_oaks <- vect('00_Shapefiles/TNC_Data/oak_plots.shp') %>%
  terra::project(unmixed_r)

# extract vegetation time series
oak_ex <- terra::extract(unmixed_r, randall_blue_oaks)
oak_ex2 <- tibble(oak_ex) %>%
  pivot_longer(
    cols = -c(ID),
    names_to = "veg_type",
    values_to = "fraction"
  )

# plot 
ggplot(oak_ex2, aes(x = veg_type, y = fraction)) +
  geom_boxplot()

# there are some oak plots with very low deciduous fraction (negative)
# proposed mask = Deciduous > 0 --> produces a map that looks a lot like the oak woodland mask from LEMMA


# Randall oak woodland mask
oak_mask <- unmixed_r$Deciduous > 0

# Prepare oak_df with meaningful factor labels
oak_df <- as.data.frame(oak_mask, xy = TRUE) |>
  rename(oak = 3) |>
  mutate(oak = ifelse(oak == 1, "Oak woodland", NA),
         oak = factor(oak, levels = "Oak woodland"))

# Convert boundary to sf and ensure CRS match
randall_sf <- st_as_sf(randall)
randall_sf <- st_transform(randall_sf, crs = st_crs(oak_mask))  # make sure CRS match

# Create the plot
p_oak <- ggplot() +
  geom_raster(data = oak_df %>% na.omit(), aes(x = x, y = y, fill = oak)) +
  geom_sf(
    data = randall_sf,
    fill = NA,
    color = "black",
    lwd = 1
  ) +
  scale_fill_manual(
    values = c("Oak woodland" = "#228B22"),
    name = NULL,
    drop = FALSE
  ) +
  labs(title = "", x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.2,
    height = unit(0.6, "cm"),
    text_cex = 1.5
  ) +
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    height = unit(2.5, "cm"),
    width = unit(1.75, "cm")
  ) +
  coord_sf(
    xlim = c(bbox["xmin"], bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"]),
    expand = FALSE
  ) 

p_oak

# save
ggsave(
  "H:/My Drive/16_Oak_Woodland/3_Figures/Manuscript_figures/Randall_Sentinel2_oak_mask.png",
  plot = p_oak,
  width = 10,
  height = 10,
  dpi = 600
)


# mask Sentinel-2 and Landsat residuals to the oak mask
# Sentinel-2
s2_residuals <- rast('Sentinel2_folder/05_Randall_temporal_stack/Randall_S2_stack_MixtureResidual.bsq')

# mask
s2_residuals <- mask(s2_residuals, oak_mask, maskvalues = 0)

# save
writeRaster(
  s2_residuals,
  'Sentinel2_folder/06_Randall_tem_unmixed/Randall_S2_stack_MixtureResidual_oak_mask',
  overwrite = TRUE,
  datatype = "FLT4S",
  filetype = 'ENVI'
)

# Landsat
ls_residuals <- rast('06_Randall_tem_unmixed/Randall_Landsat_1982_2025_TMR')

# resample the S2 mask
mask_resampled <- resample(oak_mask, ls_residuals, method = "near")

ls_residuals <- ls_residuals %>% 
  mask(mask_resampled, maskvalues = 0)

# save
writeRaster(
  ls_residuals,
  '06_Randall_tem_unmixed/Randall_Landsat_1982_2025_TMR_oak_mask',
  overwrite = TRUE,
  datatype = "FLT4S",
  filetype = 'ENVI'
)

##%######################################################%##
#                                                          #
####                 Blue oak range map                 ####
#                                                          #
##%######################################################%##

# Load the blue oak range shapefile
blue_oak_range <- vect('00_Shapefiles/blue_oak_range/BlueOakFeature.shp') %>%
  terra::project(oak_mask)

##%######################################################%##
#                                                          #
####                  TMR PCA results                   ####
#                                                          #
##%######################################################%##

tmr_pca <- rast('Sentinel2_folder/06_Randall_tem_unmixed/Randall_S2_TMR_PCA')
tmr_pca <- tmr_pca[[1:3]]
names(tmr_pca) <- c('TMR1', 'TMR2', 'TMR3')

# Band 3 is best aligned with trends of loss vs. growth, low values correspond with declines in vegetation while higher vaues correspond to increases
plot(tmr_pca[[3]], range = c(-.3, .3))


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
veg_path <- "Sentinel2_folder/05_Randall_temporal_stack/Randall_S2_stack"
tmr_pca <- "Sentinel2_folder/06_Randall_tem_unmixed/Randall_S2_TMR_PCA"

# Load raster stacks
veg_stack <- rast(veg_path)
tmr_pca <- rast(tmr_pca) 
tmr_pca[[1:3]]

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
#       event_years <- c(2015, 2017, 2019, 2021, 2023, 2025)
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
  scale_x_continuous(breaks = c(seq(2015, 2025, by = 2)))

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
  "H:/My Drive/16_Oak_Woodland/3_Figures/Manuscript_figures/Randall_S2_TMR3_residuals_quintiles.png",
  plot = combined_plot,
  width = 30,
  height = 12,
  dpi = 600
)

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
            filename = "08_Outputs/Randall_Sentinel_change_quintiles.tif",
            filetype = "GTiff",  # GeoTIFF format
            overwrite = TRUE)

writeRaster(r,
            filename = "08_Outputs/Randall_Sentinel_change_raw.tif",
            filetype = "GTiff",  # GeoTIFF format
            overwrite = TRUE)

