##%######################################################%##
#                                                          #
####                 Sentinel Unmixing                  ####
#                                                          #
##%######################################################%##

# Load libraries

library(caTools)
library(reticulate)
require(terra)
require(dplyr)
require(sf)
require(elevatr)
require(tidyverse)
require(ggspatial)


setwd('H:/My Drive/16_Oak_Woodland/1_Data/')

##%######################################################%##
#                                                          #
####            Adjust Dye Creek study area             ####
####          based on agricultural land cover          ####
#                                                          #
##%######################################################%##

s2 <- list.files('Sentinel2_folder/03_DyeCreek_envi', pattern = 'e$', full.names = TRUE)
s2 <- s2[1] %>% rast()

# Dye Creek shapefile with correct projected coordinate system -- always check!
dyecreek <- vect('00_Shapefiles/TNC_Data/Dye_Creek_Preserve_clean_proj.shp') 

nlcd <- rast('H:/My Drive/UCR Drive/Franklin_grant/project/data/landcover/NLCD2016/NLCD_2016_Land_Cover_L48_20190424.img')
dyecreek_t <- dyecreek %>% project(crs(nlcd))

nlcd_dye <- nlcd %>% crop(dyecreek_t)
nlcd_dye <- nlcd_dye %>% project(crs(dyecreek))

writeRaster(
  nlcd_dye,
  '00_Shapefiles/NLCD_mask/dye_nlcd.tif',
  overwrite = TRUE
)

nlcd_dye <- rast('00_Shapefiles/NLCD_mask/dye_nlcd.tif')

# Create a mask for **non-agricultural and developed land**
non_ag_mask <- classify(nlcd_dye, cbind(c(81, 82, 21, 22, 23, 24), NA), others=1) 

# Crop non_ag_mask to match the sentinel imagery extent
non_ag_mask <- non_ag_mask %>% crop(s2)

# Resample non_ag_mask to match s2 resolution (without reprojecting)
non_ag_mask <- resample(non_ag_mask, s2, method="near")

# irrigated areas
north <- st_read(unzip("00_Shapefiles/TNC_Data/DyeCreek_irrigated/Irrigated area North Dye Ck.kmz", exdir = tempdir()), quiet = TRUE)
south <- st_read(unzip("00_Shapefiles/TNC_Data/DyeCreek_irrigated/Irrigated area South Dye Creek.kmz", exdir = tempdir()), quiet = TRUE)

north <- vect(north) %>%
  project(dyecreek)
south <- vect(south) %>%
  project(dyecreek)


##%######################################################%##
#                                                          #
####                 Spectral Unmixing                  ####
#                                                          #
##%######################################################%##

endmember_file = '00_Endmembers/EM_MR_Sentinel.txt'
badbands_file = '00_Badbands/S2bbl.csv'


#Set unit sum constraint weight
w = 1

#Multiply EMs & images by a scale factor?
rescale = 1/10000

files <- list.files('Sentinel2_folder/03_DyeCreek_envi', pattern = 'e$', full.names = TRUE)

outfiles = paste('Sentinel2_folder/04_DyeCreek_unmixed/', basename(files),'_Unmixed',sep="")

#Read bad bands list
bands = read.csv(badbands_file, header = F, sep = ",", stringsAsFactors = F)

#Read endmember file
G = read.csv(endmember_file, header = F, sep = "")

# files with errors
# "Sentinel2_folder/03_DyeCreek_envi/S2_10TEK_20190505T184929_L2A_e"
#  "Sentinel2_folder/03_DyeCreek_envi/S2_10TEK_20240913T190011_L2A_e"

for(i in 1:length(files)) {
  
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
    
    # crop and mask to the extent of DyeCreek preserve boundaries
    u_raster <- u_raster %>%
      mask(non_ag_mask)
    
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
####                   Layer Stacking                   ####
#                                                          #
##%######################################################%##

# stack vegetation fraction images & rename with fractional year
library(lubridate)

# clear image list
ls_a <- list.files(
  'H:/My Drive/16_Oak_Woodland/1_Data/Sentinel_2/calendar_plots/DyeCreek/A',
  pattern = '.jpg',
  full.names = F
)
# remove .jpg from items in ls_a
ls_a <- gsub('.jpg', '', ls_a)
ls_a <- gsub('preview', '', ls_a)
ls_a <- gsub('L1GS', 'L2SP', ls_a)
# Remove items containing ' (1)' and (2)
ls_a <- ls_a[!grepl("\\(1\\)|\\(2\\)", ls_a)]
# Just the date (YYYYMMDD)
ls_a <- gsub("\\.", "", ls_a)

# Stack all files in "04_DyeCreek_unmixed/" that were in the original clear folder and save as Envi raster
files <- list.files('1_Data/Sentinel_2/04_DyeCreek_unmixed/', pattern = 'Unmixed$', full.names = TRUE)

# Filter files that contain any of the clear Landsat dates
# selected_files <- files[sapply(files, function(file) any(grepl(paste(ls_a, collapse = "|"), file)))]
selected_files <- files

# write to selected_files to file
writeLines(selected_files, '1_Data/Sentinel_2/DyeCreek_S2_files.txt')


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

frac_years <- calc_frac_year(selected_files)
frac_years <- round(frac_years,3)

r <- lapply(selected_files, rast)

# stack the 2nd layer in each element of r (vegetation fraction)
r2 <- lapply(r, function(x) x[[2]])

# stack the rasters
veg <- rast(r2)
names(veg) <- frac_years

# mask of NA values in veg
mask <- any(is.na(veg))

# mask out NA values
masked_veg <- mask(veg, mask, maskvalues = 1)

# crop to dye creek extent
masked_veg <- masked_veg %>%
  crop(dyecreek) %>%
  mask(dyecreek)

# elevation rasters
elev_rast <- list.files('1_Data/00_DEM/', pattern = '.tif$', full.names = TRUE) %>%
  lapply(rast)

elev_rast <- do.call(terra::mosaic, elev_rast) %>%
  project(masked_veg[[1]]) %>%
  resample(masked_veg[[1]], method = 'bilinear')
  
high_elev <- elev_rast >100
high_elev[high_elev == 0] <- NA

# mask based on DEM
masked_veg <- mask(masked_veg, high_elev)
names(masked_veg) <- frac_years


# Save as ENVI file
writeRaster(
  masked_veg,
  'Sentinel2_folder/05_DyeCreek_temporal_stack/DyeCreek_S2_stack',
  datatype = "FLT4S",
  overwrite = TRUE,
  filetype = 'ENVI'
)


# Change wavelength names in ENVI hdr file
hdr_file <- 'Sentinel2_folder/05_DyeCreek_temporal_stack/DyeCreek_S2_stack.hdr'

# Read the header file
hdr <- readLines(hdr_file)

# Add a line for wavelengths = {
hdr <- c(hdr[1:grep('data ignore value = nan', hdr)], 'wavelength = {')

# repeat lines 15-380 (based on how many images you have, always check hdr file)
for (i in 15:391) {
  hdr <- c(hdr, hdr[i])
}

# write hdr file
writeLines(hdr, 'Sentinel2_folder/05_DyeCreek_temporal_stack/DyeCreek_S2_stack.hdr')



##%######################################################%##
#                                                          #
####                 Temporal Unmixing                  ####
#                                                          #
##%######################################################%##

source('H:/My Drive/16_Oak_Woodland/2_Scripts/00_custom_functions.R')

unmix_r <- unmix(
  endmember_file = 'Sentinel_2/05_DyeCreek_temporal_stack/DyeCreek_S2_stack_TEMs.txt', # endmember file with temporal endmembers
  D1_file = 'Sentinel_2/05_DyeCreek_temporal_stack/DyeCreek_S2_stack', # raster stack with vegetation fraction time series
  n_layers = 376, # number of layers in the raster stack (should be 376 for Dye Creek)
  w = 1, # unit sum constraint weight
  rescale = 1 # rescale endmembers and raster stack to surface reflectance
)

# write temporally unmixed fractions to file
terra::writeRaster(
  unmix_r$unmixing_results,
  'Sentinel2_folder/06_DyeCreek_tem_unmixed/DyeCreek_S2_TEMs', # temporally unmixed fractions (should be 4 rasters)
  overwrite = TRUE,
  datatype = "FLT4S",
  filetype = 'ENVI'
)

# write to file (full stack)
terra::writeRaster(
  unmix_r$residuals,
  'Sentinel2_folder/06_DyeCreek_tem_unmixed/DyeCreek_S2_TMR_time_series_full', # feel free to rename - this raster contains the temporal mixture model residuals
  overwrite = TRUE,
  datatype = "FLT4S",
  filetype = 'ENVI'
)

# write pre-fire to file
terra::writeRaster(
  unmix_r$residuals[[1:349]],
  'Sentinel2_folder/06_DyeCreek_tem_unmixed/DyeCreek_S2_TMR_time_series_pre_fire',
  # feel free to rename - this raster contains the temporal mixture model residuals
  overwrite = TRUE,
  datatype = "FLT4S",
  filetype = 'ENVI'
)

# write post-fire to file
terra::writeRaster(
  unmix_r$residuals[[350:376]],
  'Sentinel2_folder/06_DyeCreek_tem_unmixed/DyeCreek_S2_TMR_time_series_post_fire',
  # feel free to rename - this raster contains the temporal mixture model residuals
  overwrite = TRUE,
  datatype = "FLT4S",
  filetype = 'ENVI'
)

tem <- rast('Sentinel_2/06_DyeCreek_tem_unmixed/DyeCreek_S2_TEMs')

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


##%######################################################%##
#                                                          #
####             Temporal endmember trends              ####
#                                                          #
##%######################################################%##


tem <- read.csv(endmember_file, header = F, sep = "", skip = 6) %>% tibble()

# rename columns
colnames(tem) <- c('Year', 'Unvegetated', 'Evergreen',  'Deciduous', 'Annual')


# plot
plot_df <- tem %>%
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
  scale_color_manual(values = c('red', 'cyan', 'green')) +  # Adjust 'black' to 'white' for visibility
  ylim(-0.05, .55) +
  scale_x_continuous(breaks = seq(2016, 2025, 1)) # change for your time period

print(plot)

# save figure
ggsave(
  "H:/My Drive/16_Oak_Woodland/3_Figures/DyeCreek_S2_pre_fire_tEMs.png",
  plot,
  width = 14,
  height = 10,
  dpi = 300
)


##%######################################################%##
#                                                          #
####                Oak Mask - Dye Creek                ####
#                                                          #
##%######################################################%##

# tem_r[[1]] = unvegetated, tem_r[[2]] = evergreen, tem_r[[3]] = deciduous, tem_r[[4]] = annual


tem_r <- rast('Sentinel2_folder/06_DyeCreek_tem_unmixed/DyeCreek_S2_TEMs')
# mask based on deciduousness
decid <- tem_r[[3]]
annual <- tem_r[[4]]
evergreen <- tem_r[[2]] # aligns best with oak woodland

# Randall oak woodland mask
oak_mask <- evergreen > 0

# Prepare oak_df with meaningful factor labels
oak_df <- as.data.frame(oak_mask, xy = TRUE) |>
  rename(oak = 3) |>
  mutate(oak = ifelse(oak == 1, "Oak woodland", NA),
         oak = factor(oak, levels = "Oak woodland"))

# Convert boundary to sf and ensure CRS match
dye_sf <- st_as_sf(dyecreek)

bbox <- st_bbox(dye_sf) # Get bounding box of the dye creek preserve

# Create the plot
p_oak <- ggplot() +
  geom_sf(
    data = dye_sf,
    fill = 'white',
    color = "black",
    lwd = 1
  ) +
  geom_raster(data = oak_df %>% na.omit(), aes(x = x, y = y, fill = oak)) +
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
    width_hint = 0.17,
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
  "H:/My Drive/16_Oak_Woodland/3_Figures/Manuscript_figures/DyeCreek_Sentinel2_oak_mask.png",
  plot = p_oak,
  width = 10,
  height = 12,
  dpi = 600
)

##%######################################################%##
#                                                          #
####       Mask residual time series to oak mask        ####
#                                                          #
##%######################################################%##

# Define input rasters
s2_paths <- list(
  full      = "Sentinel2_folder/06_DyeCreek_tem_unmixed/DyeCreek_S2_TMR_time_series_full",
  pre_fire  = "Sentinel2_folder/06_DyeCreek_tem_unmixed/DyeCreek_S2_TMR_time_series_pre_fire",
  post_fire = "Sentinel2_folder/06_DyeCreek_tem_unmixed/DyeCreek_S2_TMR_time_series_post_fire"
)

# Output directory (adjust if needed)
out_dir <- "Sentinel2_folder/06_DyeCreek_tem_unmixed"

# Load and mask each raster
s2_residuals <- lapply(s2_paths, function(path) {
  rast(path) %>%
    mask(oak_mask, maskvalues = 0)
})

# Write masked rasters to file
names(s2_residuals) <- names(s2_paths)  # label list elements
for (name in names(s2_residuals)) {
  writeRaster(
    s2_residuals[[name]],
    filename = file.path(out_dir, paste0("DyeCreek_S2_", name, "_MixtureResidual_oak_mask")),
    overwrite = TRUE,
    datatype = "FLT4S",
    filetype = "ENVI"
  )
}

# Landsat
ls_residuals <- rast('06_DyeCreek_tem_unmixed/DyeCreek_LS_1982_2025_TMR_time_series')

# resample the S2 mask
mask_resampled <- resample(oak_mask, ls_residuals, method = "near")

ls_residuals <- ls_residuals %>% 
  mask(mask_resampled, maskvalues = 0)

# save
writeRaster(
  ls_residuals,
  '06_DyeCreek_tem_unmixed/DyeCreek_LS_1982_2025_TMR_oak_mask',
  overwrite = TRUE,
  datatype = "FLT4S",
  filetype = 'ENVI'
)

# pre fire
writeRaster(
  ls_residuals[[1:816]],
  '06_DyeCreek_tem_unmixed/DyeCreek_LS_1982_2025_TMR_pre_fire_oak_mask',
  overwrite = TRUE,
  datatype = "FLT4S",
  filetype = 'ENVI'
)

# post fire
writeRaster(
  ls_residuals[[817:845]],
  '06_DyeCreek_tem_unmixed/DyeCreek_LS_1982_2025_TMR_post_fire_oak_mask',
  overwrite = TRUE,
  datatype = "FLT4S",
  filetype = 'ENVI'
)

