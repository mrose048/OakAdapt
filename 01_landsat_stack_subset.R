# Landsat Collection 2: stacking and subsetting

setwd('H:/My Drive/16_Oak_Woodland/1_Data/')
# setwd('/Users/brookerose/Library/CloudStorage/GoogleDrive-mbrose@sdsu.edu/My Drive/16_Oak_Woodland/1_Data/')

# Libraries
library(terra)
library(tidyverse)
library(mapview)

##%######################################################%##
#                                                          #
####             Randall Preserve Stacking,             ####
####              Subsetting and rescaling              ####
#                                                          #
##%######################################################%##

# Randall data
randall_files <- list.files('01_Randall_downloads/', full.names = F, recursive = T) # edit path  based on where your downloads R

# Create folders for each landsat date and place files in those folders
df_rf <- data.frame(file = randall_files, file_base = basename(randall_files), label = 'Randall')

# Identifying file name for each scene/date
folders <- str_split_fixed(df_rf$file_base, 'T1_', 2)[,1]
# add T1 to folder names
folders <- paste0(folders, 'T1')
# remove anything with pattern T2
folders <- folders[!grepl('T2', folders)]
folders <- folders[!grepl('desktop', folders)]
folders <- unique(folders)

# create folder for each folder
for (i in 1:length(folders)) {
  dir.create(
    paste0(
      'H:/My Drive/16_Oak_Woodland/1_Data/02_Randall_Landsat_raw/', # edit path based on your directory
      folders[i], '/'
    )
  )
}

# copy files to folders

# Copy files to folders only if they do not already exist
for (i in 1:nrow(df_rf)) {
  print(i)
  
  from_path <- paste0(
    'H:/My Drive/16_Oak_Woodland/1_Data/01_Randall_downloads/',
    df_rf$file[i]
  )
  
  to_path <- paste0(
    'H:/My Drive/16_Oak_Woodland/1_Data/02_Randall_Landsat_raw/',
    basename(str_split_fixed(df_rf$file[i], 'T1_', 2)[, 1]),
    'T1/',
    df_rf$file_base[i]
  )
  
  # Check if file already exists before copying
  if (!file.exists(to_path)) {
    file.copy(from = from_path, to = to_path)
  } else {
    print(paste("File already exists:", to_path))
  }
}


# clear image list
ls_a <- list.files('H:/My Drive/16_Oak_Woodland/1_Data/TNC_Landsat_Paths/Randall/A_new/', pattern = '.jpg', full.names = F)
# remove .jpg from items in ls_a
ls_a <- gsub('.jpg', '', ls_a)
ls_a <- gsub('L1TP', 'L2SP', ls_a)
ls_a <- gsub('L1GS', 'L2SP', ls_a)

# remove T2 from ls_a
ls_a <- ls_a[!grepl('T2', ls_a)]

# Remove items containing ' (1)'
ls_a <- ls_a[!grepl(" \\(1\\)$", ls_a)]

# split for just the first part
ls_a <- apply(str_split_fixed(ls_a, "_", 5)[, 1:4], 1, paste, collapse = "_")

folders <- apply(str_split_fixed(folders, "_", 5)[, 1:4], 1, paste, collapse = "_")

# which elements in ls_a are missing from folders
missing <- ls_a[!ls_a %in% folders]

# which elements in folders are missing from A_new
missing2 <- folders[!folders %in% ls_a]
# Ensure missing_2 contains only unique folder names
missing2 <- unique(missing2) # 17 dates that were downloaded 

# write to excel file
write.csv(missing, 'TNC_Landsat_Paths/Randall/missing_files.csv')

##%######################################################%##
#                                                          #
####               Stack and crop images                ####
#                                                          #
##%######################################################%##

# folder directories
folders <- list.dirs('02_Randall_Landsat_raw', full.names = F)

# Filter folders that contain any of the clear Landsat image names
selected_folders <- folders[sapply(folders, function(folder) any(grepl(paste(ls_a, collapse = "|"), folder)))]

# Print selected folders
print(selected_folders)

# Check length of selected folders
length(selected_folders)

# Check length of clear list
length(ls_a) #ok!

tnc_shape_buffer <- vect('00_Shapefiles/TNC_Data/Randall_buffer.shp')

# for clear data, crop to the extent of Randall
for (i in 1:length(selected_folders)) {
  print(i)
  if (!file.exists(paste0('03_Randall_mosaics/',selected_folders[i], '.tif'))) {
    files <- list.files(
      paste0(
        '02_Randall_Landsat_raw/',
        selected_folders[i],
        '/'
      ),
      pattern = '_SR_B',
      full.names = T
    )
    # v <- vrt(files) # converts to a virtual raster dataset using gdal
    # vc <- v %>% crop(tnc_shape_buffer)
    r <- files %>% rast()
    # sr <- c(r)
    sc <- r %>% crop(tnc_shape_buffer) 
    writeRaster(sc,
                paste0('03_Randall_mosaics/', selected_folders[i], '.tif'),
                overwrite = T)
  }
}


# for(i in 1:length(problem_files)) {
#   print(i)
#   if (!file.exists(paste0('03_Randall_mosaics/', problem_files[i], '.tif'))) {
#     files <- list.files(
#       paste0(
#         '02_Randall_Landsat_raw/',
#         problem_files[i],
#         '/'
#       ),
#       pattern = '_SR_B',
#       full.names = T
#     )
#     r <- files %>% rast()
#     sc <- r %>% crop(tnc_shape_buffer) 
#     writeRaster(sc,
#                 paste0('03_Randall_mosaics/', problem_files[i], '.tif'),
#                 overwrite = T)
#   }
# }

#########################################

# At the time of analysis (March 2025), the band multiplier and additive factor were the same across all bands and satellites for Landsat Level 2 Collection 2 surface reflectance products

# Scale reflectance
files <- list.files(
  'H:/My Drive/16_Oak_Woodland/1_Data/03_Randall_mosaics/',
  full.names = T,
  pattern = 'T1.tif$'
)
og_dir <- list.dirs('H:/My Drive/16_Oak_Woodland/1_Data/02_Randall_Landsat_raw/',
                    full.names = T)

b_mult = 0.0000275 ## CHECK
b_add = -0.2

# Start timing
start_time <- Sys.time()

for (i in 1:length(files)) {
  # meta = list.files(og_dir[i+1], pattern = 'MTL.txt', full.names = T)
  print(i)
  # read txt file
  # meta = readLines(meta)
  
  if (!file.exists(gsub('.tif', '_Scale', files[i]))) {
    # Paste scale between file name and extension
    outfile <- gsub('.tif', '_Scale', files[i])
    
    r <- files[i] %>% rast()
    rout <- r # creating a raster object for saving the scaled rasters
    
    rout <- (r * b_mult) + b_add
    
    writeRaster(rout, outfile, filetype = 'ENVI', overwrite = T)
    gc()
  }
}

# End timing
end_time <- Sys.time()

# Calculate and print elapsed time
elapsed_time <- end_time - start_time
print(paste("Total time elapsed:", elapsed_time))

##%######################################################%##
#                                                          #
####    Dye Creek Stacking, Subsetting and rescaling    ####
#                                                          #
##%######################################################%##

# Define base directories
dye_dir <- '01_DyeCreek_downloads/'
google_drive_dir <- 'H:/My Drive/16_Oak_Woodland/1_Data/TNC_Landsat_Paths/Dye_Creek/'
amelie_hd_dir <- 'F:/Dye Creek TNC/Bulk Downloads/'

# Remove duplicate files with "(1)" or "(2)" in the name
dye_files <- list.files(dye_dir, full.names = TRUE, recursive = T)
dup_patterns <- c(" \\(1\\)\\.", " \\(2\\)\\.")

for (pattern in dup_patterns) {
  file.remove(dye_files[grepl(pattern, dye_files)])
}

# Get unique filenames from Google Drive
gd_files <- list.files(google_drive_dir,
                       recursive = TRUE,
                       full.names = FALSE) %>%
  basename() %>%
  strsplit("_") %>%
  sapply(function(x)
    paste(x[1], x[4], sep = "_")) %>%
  unique()

# Get unique filenames from Amélie's hard drive
ah_files <- list.files(amelie_hd_dir, full.names = FALSE) %>%
  strsplit("_") %>%
  sapply(function(x)
    paste(x[1], x[4], sep = "_")) %>%
  unique()

# Identify missing files
missing1 <- setdiff(ah_files, gd_files)  # Files on Amelie's hard drive but not in Google Drive

# Get unique filenames from A folder
a_files <- list.files(file.path(google_drive_dir, "A_new"),
                      pattern = ".jpg",
                      full.names = FALSE) %>%
  strsplit("_") %>%
  sapply(function(x)
    paste(x[1], x[4], sep = "_")) %>%
  unique()

# Files in A folder missing from Google Drive
# missing2 <- setdiff(a_files, gd_files) # none

# Save missing files list
# write.csv(missing2, file.path(google_drive_dir, "files_in_A_missing_from_Google_Drive.csv"), row.names = FALSE)

# Files in Google Drive missing from A folder
missing3 <- setdiff(gd_files, a_files) # fine, can remove

# Function to check missing files in B and C folders
check_missing_in_folder <- function(folder_name, missing_list) {
  folder_files <- list.files(
    file.path(google_drive_dir, folder_name, folder_name),
    pattern = ".jpg",
    full.names = FALSE
  ) %>%
    strsplit("_") %>%
    sapply(function(x)
      paste(x[1], x[4], sep = "_")) %>%
    unique()
  intersect(missing_list, folder_files)
}

# Files missing from A but present in B and C
missing4 <- check_missing_in_folder("B", missing3)  # Files in B folder
missing5 <- check_missing_in_folder("C", missing3)  # Files in C folder

# Save results
write.csv(
  missing4,
  file.path(google_drive_dir, "files_in_Google_Drive_that_are_in_B.csv"),
  row.names = FALSE
)
write.csv(
  missing5,
  file.path(google_drive_dir, "files_in_Google_Drive_that_are_in_C.csv"),
  row.names = FALSE
)

# Files missing from A, B, and C
missing6 <- setdiff(missing3, c(missing4, missing5))

# Save final list of missing files
write.csv(
  missing6,
  file.path(
    google_drive_dir,
    "files_in_Google_Drive_that_are_missing_from_A_B_C.csv"
  ),
  row.names = FALSE
)

# Get all files from Dye Creek Downloads
dye_files <- list.files(dye_dir, full.names = TRUE, recursive = TRUE)

# Extract satellite name and date (assuming consistent structure)
dye_info <- sapply(strsplit(basename(dye_files), "_"), function(x)
  paste(x[1], x[4], sep = "_"))

# Filter dye_files to include only those matching both satellite and date in a_files OR those that were collected in 2025
matched_dye_files <- dye_files[dye_info %in% a_files | grepl("2025", dye_info)]

# Length
length(matched_dye_files)

# Create folders for each landsat date and place files in those folders
df_dc <- data.frame(file = matched_dye_files, label = 'Dye_Creek')
df_dc$file_base <- basename(df_dc$file)

# Identifying file name for each scene/date
folders <- str_split_fixed(df_dc$file_base, 'T1_', 2)[, 1]
# add T1 to folder names
folders <- paste0(folders, 'T1')
# remove anything with pattern T2
folders <- folders[!grepl('T2', folders)]
folders <- folders[!grepl('desktop', folders)]
folders <- unique(folders)

length(unique(folders)) # 26 dates missing that were labeled clear but not in folders

# create folder for each folder
for (i in 1:length(folders)) {
  dir.create(
    paste0(
      'H:/My Drive/16_Oak_Woodland/1_Data/02_DyeCreek_Landsat_raw/',
      # edit path based on your directory
      folders[i],
      '/'
    )
  )
}

# copy files to folders
for (i in 1:nrow(df_dc)) {
  print(i)
  
  from_path <- file.path('H:/My Drive/16_Oak_Woodland/1_Data', df_dc$file[i])
  
  to_path <- paste0(
    'H:/My Drive/16_Oak_Woodland/1_Data/02_DyeCreek_Landsat_raw/',
    basename(str_split_fixed(df_dc$file_base[i], 'T1_', 2)[, 1]),
    'T1/',
    df_dc$file_base[i]
  )
  
  # Check if file already exists before copying
  if (!file.exists(to_path)) {
    file.copy(from = from_path, to = to_path)
  } else {
    print(paste("File already exists:", to_path))
  }
}

##%######################################################%##
#                                                          #
####           Stack and subset - Dye Creek             ####
#                                                          #
##%######################################################%##

# Dye Creek shapefile with correct projected coordinate system -- always check!
dye_shape <- vect('00_Shapefiles/TNC_Data/Dye_Creek_Preserve_clean_proj.shp') 

# Get unique filenames from A folder
a_files <- list.files(file.path(google_drive_dir, "A_new"), pattern = ".jpg", full.names = FALSE) %>%
  strsplit("_") %>%
  sapply(function(x) paste(x[1], x[4], sep = "_")) %>%
  unique()

folders <- list.dirs('02_DyeCreek_Landsat_raw', full.names = F)
# Remove duplicate files with "(1)" or "(2)" in the name
# dye_files <- list.files('02_DyeCreek_Landsat_raw', full.names = TRUE, recursive = TRUE)
# dup_patterns <- c(" \\(1\\)\\.", " \\(2\\)\\.", " \\(3\\)\\.", " \\(4\\)\\.")
# 
# for (pattern in dup_patterns) {
#   file.remove(dye_files[grepl(pattern, dye_files)])
# }

# Extract satellite and date information from the folder names (assuming they follow the same structure)
folder_info <- sapply(strsplit(folders, "_"), function(x) paste(x[1], x[4], sep = "_"))

# Select folders where the satellite and image date match those in a_files
matching_folders <- folders[folder_info %in% a_files | grepl("2025", folder_info)]

# Print the matching folders
print(matching_folders) # 861

# Start timing
start_time <- Sys.time()

# Initialize an empty vector to store skipped folder names
skipped_folders <- c()

for (i in seq_along(matching_folders)) {
  print(i)
  
  output_path <- paste0('03_DyeCreek_mosaics/', matching_folders[i], '.tif')
  
  # Skip processing if the output file already exists
  if (file.exists(output_path)) {
    message("Output already exists for ", matching_folders[i], " — skipping.")
    next
  }
  
  # List the .TIF files
  files <- list.files(
    paste0('02_DyeCreek_Landsat_raw/', matching_folders[i], '/'),
    pattern = '_SR_B.*\\.TIF$',
    full.names = TRUE
  )
  
  # Try reading the raster files
  r <- tryCatch({
    rast(files)
  }, error = function(e) {
    message("Error reading raster files for folder: ", matching_folders[i], " — skipping.")
    skipped_folders <- c(skipped_folders, matching_folders[i])
    return(NULL)
  })
  
  # If reading failed, continue to next iteration
  if (is.null(r)) next
  
  # Layer count checks
  if (grepl("LC08|LC09", matching_folders[i])) {
    if (nlyr(r) != 7) {
      message("Skipping ", matching_folders[i], ": expected 7 layers, found ", nlyr(r))
      skipped_folders <- c(skipped_folders, matching_folders[i])
      next
    }
  } else if (grepl("LT05|LT04|LE07", matching_folders[i])) {
    if (nlyr(r) != 6) {
      message("Skipping ", matching_folders[i], ": expected 6 layers, found ", nlyr(r))
      skipped_folders <- c(skipped_folders, matching_folders[i])
      next
    }
  }
  
  # Crop and save
  sc <- crop(r, dye_shape)
  writeRaster(sc, output_path, overwrite = TRUE)
}


# After the loop, save the list of skipped folders to a CSV file
write.csv(skipped_folders, 'skipped_folders.csv', row.names = F) # 10


# End timing
end_time <- Sys.time()

# Calculate and print elapsed time
elapsed_time <- end_time - start_time
print(paste("Total time elapsed:", elapsed_time))

##%######################################################%##
#                                                          #
####           Scale Reflectance - Dye Creek            ####
#                                                          #
##%######################################################%##

# 

# At the time of analysis, the band multiplier and additive factor were the same across all bands and satellites for Landsat Level 2 Collection 2 surface reflectance products

# Define scale factors for Landsat Level 2 Collection 2 surface reflectance
b_mult <- 0.0000275
b_add <- -0.2

# List input raster files (ending with T1.tif)
files <- list.files(
  'H:/My Drive/16_Oak_Woodland/1_Data/03_DyeCreek_mosaics/',
  full.names = TRUE,
  pattern = 'T1.tif$'
)

# Start timing
start_time <- Sys.time()

for (i in seq_along(files)) {
  print(i)
  
  # Define output filename (ENVI format uses .dat as default extension)
  outfile <- gsub("\\.tif$", "_Scale", files[i])
  
  # Skip if output already exists
  if (file.exists(outfile)) {
    message("Already processed: ", basename(files[i]), " — skipping.")
    next
  }
  
  # Read the raster and apply scaling
  r <- tryCatch({
    rast(files[i])
  }, error = function(e) {
    message("Error reading: ", files[i])
    return(NULL)
  })
  
  if (is.null(r)) next
  
  # Apply scaling
  rout <- (r * b_mult) + b_add
  
  # Save as ENVI binary (.dat + .hdr)
  writeRaster(rout, outfile, filetype = "ENVI", overwrite = TRUE)
  
  # Free up memory
  gc()
}

# End timing
end_time <- Sys.time()

# Report elapsed time
elapsed_time <- end_time - start_time
print(paste("Total time elapsed:", elapsed_time))

##%######################################################%##
#                                                          #
####             Calendar plot - Dye Creek              ####
#                                                          #
##%######################################################%##

# 853 dates

# Calendar Plot for Landsat Observations at Dye Creek

library(tidyverse)
library(stringr)

# Prepare dataframe with metadata extracted from file names
df2 <- tibble(file = basename(files)) %>%
  mutate(
    satellite = str_extract(file, "^[^_]+"),
    satellite = factor(satellite, levels = c('LT04', 'LT05', 'LE07', 'LC08', 'LC09')),
    pathrow   = str_split_fixed(file, '_', 7)[, 3],
    date      = str_split_fixed(file, '_', 7)[, 4],
    year      = str_sub(date, 1, 4),
    month     = str_sub(date, 5, 6),
    day       = str_sub(date, 7, 8),
    doy       = as.numeric(format(as.Date(date, format = '%Y%m%d'), '%j'))
  )

# Create calendar-style scatter plot
p <- ggplot(df2, aes(x = doy, y = as.numeric(year), color = satellite)) +
  geom_point(size = 4) +
  theme_minimal(base_size = 16) +
  labs(
    title = paste('# of images:', nrow(df2)),
    x = 'Day of Year',
    y = 'Year',
    color = 'Satellite'
  ) +
  scale_color_manual(values = c(
    "#1B9E77", "#7570B3", "#D95F02", "#E7298A", "#66A61E"
  )) +
  scale_y_continuous(breaks = seq(1980, 2024, 10)) +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_line(color = "gray50", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 22, color = "white"),
    axis.title = element_text(size = 24, face = "bold", color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = 'bottom',
    legend.title = element_text(size = 20, face = "bold", color = "white"),
    legend.text = element_text(size = 20, color = "white"),
    legend.key = element_rect(fill = 'black', color = NA),
    legend.background = element_rect(fill = 'black', color = NA)
  )

# Display plot
print(p)

# Save to file
ggsave(
  filename = 'H:/My Drive/16_Oak_Woodland/3_Figures/Manuscript_figures/DyeCreek_Landsat_calendar_plot.png',
  plot = p,
  width = 8.5,
  height = 8,
  dpi = 300
)
