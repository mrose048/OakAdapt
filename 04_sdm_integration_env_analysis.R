##%######################################################%##
#                                                          #
####             Integrate oak cover change             ####
####        products and environmental analysis         ####
#                                                          #
##%######################################################%##

# Load required libraries
library(terra)
library(randomForest)
library(tidyverse)
library(rpart)
library(rpart.plot)
require(elevatr)
require(whitebox)
require(sf)
require(ggnewscale)
require(ggspatial)
require(ggrepel)


setwd("H:/My Drive/16_Oak_Woodland/")

source('2_Scripts/00_functions.R')

# Preserve shapefiles
randall <- vect('1_Data/00_Shapefiles/study_areas/Randall.shp')
dyecreek <- vect('1_Data/00_Shapefiles/study_areas/DyeCreek.shp')

# Blue oak range maps
qd_randall <- vect('1_Data/00_Shapefiles/blue_oak_range/BlueOakRandall.shp')
qd_dyecreek <- vect('1_Data/00_Shapefiles/blue_oak_range/BlueOakDyeCreek.shp')


##%######################################################%##
#                                                          #
####    Load environmental data (baseline and past)     ####
#                                                          #
##%######################################################%##

# baseline environmental variables
# env <- list.files(
#   "H:/My Drive/12_T&E_Plants_SDMs/Models/1_Inputs/2_Predictors/1_Current/",
#   full.names = T
# ) %>%
#   terra::rast() 


# baseline environmental variables
env2 <- list.files(
  "H:/My Drive/UCR Drive/Franklin_grant/project/1-NSF_spatial_and_species_traits/1_Inputs/2_Predictors/1_Current/",
  pattern = '.tif$',
  full.names = T
) %>%
  terra::rast()

# For mac
env2 <- list.files(
  "/Users/brookerose/Library/CloudStorage/GoogleDrive-mbrose@sdsu.edu/My Drive/UCR Drive/Franklin_grant/project/1-NSF_spatial_and_species_traits/1_Inputs/2_Predictors/1_Current/",
  pattern = '.tif$',
  full.names = T
) %>%
  terra::rast()


# Elevation - Randall
# randall_elev <- get_elev_raster(st_as_sf(randall), z = 14, prj = crs(randall))
# randall_elev <- project(rast(randall_elev), crs(randall))

# save
# writeRaster(randall_elev, '1_Data/00_DEM/Randall_dem_10m.tif')
# randall_elev <- rast('1_Data/00_DEM/Randall_dem_10m.tif')
# 
# # Calculate terrain metrics (10m)
# slope <- terrain(randall_elev, v = "slope", unit = "degrees")
# aspect <- terrain(randall_elev, v = "aspect", unit = "radians")
# northness <- cos(aspect)
# names(northness) <- "northness"
# eastness <- sin(aspect)
# names(eastness) <- "eastness"
# tri <- terrain(randall_elev, v = "TRI")
# tpi <- terrain(randall_elev, v = "TPI")
# flow <- terrain(randall_elev, v = "flowdir")
# 
# # Stack and save terrain attributes
# terrain_metrics <- c(randall_elev, slope, aspect, northness, eastness, tri, tpi, flow)
# names(terrain_metrics) <- c(
#   "elevation", "slope", "aspect", "northness", "eastness", "tri", "tpi", "flowdir"
# )
# writeRaster(terrain_metrics, '1_Data/00_DEM/randall_terrain_10m.tif', datatype = 'FLT4S', overwrite = T)

terrain_metrics <- rast('1_Data/00_DEM/randall_terrain_10m.tif')

# Dye Creek
# dyecreek_elev <- get_elev_raster(st_as_sf(dyecreek), z = 14, prj = crs(dyecreek))
# dyecreek_elev <- project(rast(dyecreek_elev), crs(dyecreek))
# # Calculate terrain metrics (10m)
# slope_d <- terrain(dyecreek_elev, v = "slope", unit = "degrees")
# aspect_d <- terrain(dyecreek_elev, v = "aspect", unit = "radians")
# northness_d <- cos(aspect_d)
# names(northness_d) <- "northness"
# eastness_d <- sin(aspect_d)
# names(eastness_d) <- "eastness"
# tri_d <- terrain(dyecreek_elev, v = "TRI")
# tpi_d <- terrain(dyecreek_elev, v = "TPI")
# flow_d <- terrain(dyecreek_elev, v = "flowdir")

# Stack and save terrain attributes
# dye_creek_terrain <- c(dyecreek_elev, slope_d, aspect_d, northness_d, eastness_d, tri_d, tpi_d, flow_d)
# names(dye_creek_terrain) <- c(
#   "elevation", "slope", "aspect", "northness", "eastness", "tri", "tpi", "flowdir"
# )
# writeRaster(dye_creek_terrain, '1_Data/00_DEM/dyecreek_terrain_10m.tif', datatype = 'FLT4S', overwrite = T)

dye_creek_terrain <- rast('1_Data/00_DEM/dyecreek_terrain_10m.tif')

# Historical BCM
# bcm <- list.files(
#   "00_Habitat_maps/Climate_data/30y-selected/",
#   full.names = T,
#   pattern = '.asc$'
# ) %>%
#   terra::rast()
# 
# # set -9999 values to na
# bcm[bcm == -9999] <- NA

##%######################################################%##
#                                                          #
####                  Oak cover change                  ####
#                                                          #
##%######################################################%##

library(terra)
library(dplyr)

# Helper function to project, rasterize, and mask
rasterize_and_mask <- function(raster_path, shapefile, shapefile_name) {
  r <- rast(raster_path)
  shp_proj <- project(shapefile, r)
  shp_proj$presence <- 1
  r_masked <- mask(r, rasterize(shp_proj, r, field = "presence", background = NA, touches = TRUE))
  message(paste("✓", shapefile_name, "-", basename(raster_path), "processed"))
  return(r_masked)
}

# Load shapefiles (assumed already loaded: qd_dyecreek and qd_randall)

# ---- Dye Creek ----
d_lsc <- rasterize_and_mask("1_Data/08_Outputs/DyeCreek_Landsat_change_raw.tif", qd_dyecreek, "DyeCreek Landsat")
d_sc  <- rasterize_and_mask("1_Data/08_Outputs/DyeCreek_S2_change_raw.tif",     qd_dyecreek, "DyeCreek Sentinel")
d_sc2 <- rasterize_and_mask("1_Data/08_Outputs/DyeCreek_S2_change_raw_TMR1.tif", qd_dyecreek, "DyeCreek Sentinel 2")

# ---- Randall ----
r_lsc <- rasterize_and_mask("1_Data/08_Outputs/Randall_Landsat_change_raw.tif", qd_randall, "Randall Landsat")
r_sc  <- rasterize_and_mask("1_Data/08_Outputs/Randall_Sentinel_change_raw.tif", qd_randall, "Randall Sentinel")


qd_randall <- project(qd_randall, r_lsc)
qd_randall$presence <- 1
qd_ran_r <- rasterize(qd_randall, r_lsc, field = "presence", background = NA, touches = TRUE)


##%######################################################%##
#                                                          #
####        Blue oak species distribution models        ####
#                                                          #
##%######################################################%##

# Define the base path and projection
base_path <- "H:/My Drive/UCR Drive/Franklin_grant/project/1-NSF_spatial_and_species_traits/2_Outputs/9_Final_SDM/2_SDM_LU/Quercus douglasii"

# Base path on Mac
base_path <- "/Users/brookerose/Library/CloudStorage/GoogleDrive-mbrose@sdsu.edu/My Drive/UCR Drive/Franklin_grant/project/1-NSF_spatial_and_species_traits/2_Outputs/9_Final_SDM/2_SDM_LU/Quercus douglasii"

crs_proj <- "+proj=longlat +datum=NAD83"

# Helper function to load and reproject raster
load_projected_raster <- function(subfolder, filename) {
  rast(file.path(base_path, subfolder, filename)) %>%
    project(crs_proj)
}

# Load rasters
fut1 <- load_projected_raster("05_hades_rcp85", "Quercus douglasii_2055.asc")
fut2 <- load_projected_raster("03_cnrm_rcp85", "Quercus douglasii_2055.asc")
curr <- load_projected_raster("01_current", "Quercus douglasii_1995.asc")
fut3 <- mean(c(fut1, fut2), na.rm = TRUE)

# Crop fut3 to randall extent
fut3_randall <- fut3 %>% project(randall) %>% crop(randall)
curr_randall <- curr %>% project(randall) %>% crop(randall)

# Crop fut3 to dye creek extent
fut3_dye <- fut3 %>% project(dyecreek) %>% crop(dyecreek)
curr_dye <- curr %>% project(dyecreek) %>% crop(dyecreek)

##%######################################################%##
#                                                          #
####    Resample habitat suitability to 10m and 30m     ####
#                                                          #
##%######################################################%##

# Resample to 10m for Sentinel-2
fut3_randall_s2 <- fut3_randall %>% resample(r_sc, method = "bilinear")
fut3_dye_s2 <- fut3_dye %>% resample(d_sc, method = "bilinear")
cur_randall_s2 <- curr_randall %>% resample(r_sc, method = 'bilinear')
curr_dye_s2 <- curr_dye %>% resample(d_sc, method = 'bilinear')

# Resample to 30m for Landsat
fut3_randall_lsc <- fut3_randall %>% resample(r_lsc, method = "bilinear")
fut3_dye_lsc <- fut3_dye %>% resample(d_lsc, method = "bilinear")


##%######################################################%##
#                                                          #
####     Resample environmental data to 10m and 30m     ####
#                                                          #
##%######################################################%##

randall_env2 <- env2 %>% project(randall) %>% crop(randall)
dyecreek_env2 <- env2 %>% project(dyecreek) %>% crop(dyecreek)

# Resample to 10m for Sentinel-2
randall_envs2_resamp <- randall_env2 %>% resample(r_sc, method = "bilinear")
dyecreek_envs2_resamp <- dyecreek_env2 %>% resample(d_sc, method = "bilinear")

# Resample to 30m for Landsat
randall_envls_resamp <- randall_env2 %>% resample(r_lsc, method = "bilinear")
dyecreek_envls_resamp <- dyecreek_env2 %>% resample(d_lsc, method = "bilinear")

# Terrain data resample to 30 for Landsat
randall_terrain_ls <- resample(terrain_metrics, r_lsc, method = "bilinear")
dye_terrain_ls <- resample(dye_creek_terrain, d_lsc, method = "bilinear")

# Terrain data resample to 10 for Sentinel-2
terrain_metrics_s2 <- resample(terrain_metrics, r_sc, method = "bilinear")
dye_terrain_s2 <- resample(dye_creek_terrain, d_sc, method = "bilinear")

# Randall Predictor stacks
## Landsat
r_lsc_pred <- c(randall_envls_resamp, randall_terrain_ls) %>%
  subset(!names(.) %in% c("ppt_jja", "aet", "slope", "flowdir", "aspect", "terrain", "northness", "eastness"))

## Sentinel-2
r_sc_pred <- c(randall_envs2_resamp, terrain_metrics_s2) %>%
  subset(!names(.) %in% c("ppt_jja", "aet", "slope", "flowdir", "aspect", "terrain", "northness", "eastness"))

# Dye Creek Predictor stacks
## Landsat
d_lsc_pred <- c(dyecreek_envls_resamp, dye_terrain_ls) %>%
  subset(!names(.) %in% c("ppt_jja", "aet", "slope", "flowdir", "aspect", "terrain", "northness", "eastness"))

## Sentinel-2
d_sc_pred <- c(dyecreek_envs2_resamp, dye_terrain_s2) %>%
  subset(!names(.) %in% c("ppt_jja", "aet", "slope", "flowdir", "aspect", "terrain", "northness", "eastness"))

##%######################################################%##
#                                                          #
####        Hillshade rasters for visualization         ####
#                                                          #
##%######################################################%##

make_hillshade_df <- function(slope_rast, aspect_rast, shape, angle = 45, direction = 315, w = 5) {
  hs <- shade(slope = slope_rast, aspect = aspect_rast, angle = angle, direction = direction)
  hs <- hs %>% crop(shape) %>% mask(shape)
 # hs_smooth <- focal(hs, w = matrix(1, w, w), fun = mean, na.policy = "omit")
  df <- as.data.frame(hs, xy = TRUE, na.rm = TRUE)
  colnames(df)[3] <- "hillshade"
  return(df)
}

hillshade_dfs <- list(
  randall_lsc = make_hillshade_df(
    slope = terrain(randall_terrain_ls$elevation, "slope", unit = "radians"),
    aspect = terrain(randall_terrain_ls$elevation, "aspect", unit = "radians"),
    shape = randall
  ),
  randall_s2 = make_hillshade_df(
    slope = terrain(terrain_metrics_s2$elevation, "slope", unit = "radians"),
    aspect = terrain(terrain_metrics_s2$elevation, "aspect", unit = "radians"),
    shape = randall
  ),
  dye_lsc = make_hillshade_df(
    slope = terrain(dye_terrain_ls$elevation, "slope", unit = "radians"),
    aspect = terrain(dye_terrain_ls$elevation, "aspect", unit = "radians"),
    shape = dyecreek
  ),
  dye_s2 = make_hillshade_df(
    slope = terrain(dye_terrain_s2$elevation, "slope", unit = "radians"),
    aspect = terrain(dye_terrain_s2$elevation, "aspect", unit = "radians"),
    shape = dyecreek
  )
)


##%######################################################%##
#                                                          #
####        Select relevant predictor variables         ####
#                                                          #
##%######################################################%##
 
library(ggcorrplot)

# Randall Landsat

# 1. Stack response and predictor variables, mask them
full_stack <- c(r_lsc, fut3_randall_lsc, r_lsc_pred)
names(full_stack)[1:2] <- c("oak_change", "suitability")
full_stack <- mask(full_stack, r_lsc)

# 2. Extract values as dataframe
env_df <- as.data.frame(full_stack, na.rm = TRUE)

# 3. Compute full correlation matrix
cor_mat <- cor(env_df, use = "complete.obs")

# 4. Plot correlations
corr_plot <- ggcorrplot(
  cor_mat,
  lab = TRUE,
  lab_size = 4,
  colors = c("blue", "white", "red"),
  title = "Correlation with Oak Cover Change",
  show.legend = TRUE,
  type = "lower"
)

print(corr_plot)

# Randall Sentinel
# 1. Stack response and predictor variables, mask them
full_stack_s2 <- c(r_sc, fut3_randall_s2, r_sc_pred)
names(full_stack_s2)[1:2] <- c("oak_change", "suitability")
full_stack_s2 <- mask(full_stack_s2, r_sc)

# 2. Extract values as dataframe
env_df_s2 <- as.data.frame(full_stack_s2, na.rm = TRUE)

# 3. Compute full correlation matrix
cor_mat_s2 <- cor(env_df_s2, use = "complete.obs")

# 4. Plot correlations
corr_plot_s2 <- ggcorrplot(
  cor_mat_s2,
  lab = TRUE,
  lab_size = 4,
  colors = c("blue", "white", "red"),
  title = "Correlation with Oak Cover Change",
  show.legend = TRUE,
  type = "lower"
)

print(corr_plot_s2)

# Dye Creek Landsat
# 1. Stack response and predictor variables, mask them
full_stack_dye <- c(d_lsc, fut3_dye_lsc, d_lsc_pred)
names(full_stack_dye)[1:2] <- c("oak_change", "suitability")
full_stack_dye <- mask(full_stack_dye, d_lsc)

# 2. Extract values as dataframe
env_df_dye <- as.data.frame(full_stack_dye, na.rm = TRUE)
# 3. Compute full correlation matrix
cor_mat_dye <- cor(env_df_dye, use = "complete.obs")
# 4. Plot correlations
corr_plot_dye <- ggcorrplot(
  cor_mat_dye,
  lab = TRUE,
  lab_size = 4,
  colors = c("blue", "white", "red"),
  title = "Correlation with Oak Cover Change",
  show.legend = TRUE,
  type = "lower"
)
print(corr_plot_dye)

# Dye Creek Sentinel
# 1. Stack response and predictor variables, mask them
full_stack_dye_s2 <- c(d_sc, fut3_dye_s2, d_sc_pred)
names(full_stack_dye_s2)[1:2] <- c("oak_change", "suitability")
full_stack_dye_s2 <- mask(full_stack_dye_s2, d_sc)

# 2. Extract values as dataframe
env_df_dye_s2 <- as.data.frame(full_stack_dye_s2, na.rm = TRUE)
# 3. Compute full correlation matrix
cor_mat_dye_s2 <- cor(env_df_dye_s2, use = "complete.obs")
# 4. Plot correlations
corr_plot_dye_s2 <- ggcorrplot(
  cor_mat_dye_s2,
  lab = TRUE,
  lab_size = 4,
  colors = c("blue", "white", "red"),
  title = "Correlation with Oak Cover Change",
  show.legend = TRUE,
  type = "lower"
)
print(corr_plot_dye_s2)

##%######################################################%##
#                                                          #
####             K-means clustering of Oak              ####
####     Resilience Scores and Habitat Suitability      ####
#                                                          #
##%######################################################%##

library(terra)
library(tidyverse)
library(cluster)
library(rpart)
library(rpart.plot)
library(gridExtra)
library(sf)
library(ggspatial)
library(ggnewscale)
library(ggrepel)

# Create one viridis color per class
viridis_colors <- viridis(5, option = "D", direction = -1)

# Assign each class a single color wrapped in a list
box_palette <- as.list(viridis_colors)

run_kmeans_analysis <- function(oak_rast, suit_rast, mask_rast, k, seed = 123, pred_rast, shape, hill_df, scale = FALSE, cluster_order = NULL) {
  
  # 1. Mask and stack rasters
  oak_mask <- mask(oak_rast, mask_rast)
  suit_mask <- mask(suit_rast, mask_rast)
  
  # 2. Standardize rasters
  # if scale = TRUE
  if (scale) {
    oak_scaled <- (oak_mask - global(oak_mask, "mean", na.rm=TRUE)[[1]]) / global(oak_mask, "sd", na.rm=TRUE)[[1]]
    suit_scaled <- (suit_mask - global(suit_mask, "mean", na.rm=TRUE)[[1]]) / global(suit_mask, "sd", na.rm=TRUE)[[1]]
  } else {
    oak_scaled <- oak_mask
    suit_scaled <- suit_mask
  }
  
  combo_stack <- c(oak_scaled, suit_scaled)
  names(combo_stack) <- c("oak_change", "suitability")
  
  # 3. Extract training data
  training_df <- as.data.frame(combo_stack, xy = TRUE, na.rm = TRUE, cells = TRUE)
  
  # round oak_change and suitability to 5 decimal points
  training_df <- training_df %>%
    mutate(
      oak_change = round(oak_change, 5),
      suitability = round(suitability, 5)
    )
  
  # 4. K-means clustering (already scaled, so don't re-scale again)
  set.seed(seed)
  km_res <- kmeans(training_df[, c("oak_change", "suitability")], centers = k, nstart = 25)
  training_df$class <- as.factor(km_res$cluster)
  
  # 5. Reorder class labels based on combined mean values
  # 5. Reorder cluster labels
  cluster_stats <- training_df %>%
    group_by(class) %>%
    summarize(
      mean_suit = mean(suitability, na.rm = TRUE),
      mean_oak  = mean(oak_change, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(mean_combo = mean_suit + mean_oak)
  
  if (is.null(cluster_order)) {
    # Default ordering: descending by combined mean
    cluster_stats <- cluster_stats %>%
      arrange(desc(mean_combo)) %>%
      mutate(new_class = row_number())
  } else {
    # User-defined ordering
    if (length(cluster_order) != k) {
      stop("cluster_order must be the same length as k")
    }
    cluster_stats <- cluster_stats %>%
      mutate(new_class = cluster_order[match(class, sort(unique(class)))])
  }
  
  # Map new labels
  class_map <- cluster_stats %>% select(class, new_class)
  
  training_df <- training_df %>%
    left_join(class_map, by = "class") %>%
    mutate(class = as.factor(new_class)) %>%
    select(-new_class)
  
  # 6. Assign cluster labels to raster
  cluster_rast <- rast(combo_stack[[1]])  # create empty raster with same structure
  values(cluster_rast)[training_df$cell] <- as.numeric(training_df$class)
  cluster_rast <- mask(cluster_rast, mask_rast)
  
  # 7. Create ggplot biplot
  biplot <- ggplot(training_df %>%  slice_sample(n = min(10000, nrow(training_df))), aes(x = oak_change, y = suitability, color = class)) +
    geom_point(alpha = 0.5, size = 1.5) +
    stat_ellipse(aes(group = class), level = 0.95, size = 1, alpha = 0.6) +
    scale_fill_viridis_d(option = "D", name = "Class", guide = 'none') +
    scale_color_viridis_d(option = "D", direction = -1, guide = 'none') +
    theme_minimal(base_size = 14) +
    labs(
      x = "TMR 3",
      y = "Future Habitat Suitability"
    ) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20),
      panel.grid.minor = element_blank()
    )
  
  # 8. Create spatial plot
  class_df <- as.data.frame(cluster_rast, xy = TRUE, na.rm = TRUE)
  names(class_df)[3] <- "class"
  class_df$class <- as.factor(class_df$class)
  
  spatial_plot <- ggplot() +
    # Hillshade as background
    geom_raster(data = hill_df, aes(x = x, y = y, fill = hillshade)) +
    scale_fill_gradient(low = "white", high = "gray60", guide = "none") +  # soft shading
    new_scale_fill() +
    
    geom_raster(data = class_df, aes(x = x, y = y, fill = class)) +
    scale_fill_viridis_d(option = "D",
                         name = "Class",
                         direction = -1) +
    geom_sf(
      data = st_as_sf(shape),
      fill = NA,
      color = "black",
      lwd = 1
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
      height = unit(2, "cm"),
      width = unit(1.5, "cm")
    ) +
    coord_sf() +
    theme_minimal() +
    labs(title = NULL, x = NULL, y = NULL) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.position = 'bottom',
      legend.direction = 'horizontal',
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 24)
    )
  
  # 9. PCA on environmental predictors
  
  # Extract environmental predictor values from raster (cells to match class)
  pred_values <- as.data.frame(pred_rast %>% mask(mask_rast), cells = TRUE, na.rm = TRUE)
  
  # Join environmental values with class labels using cell index
  training_df_env <- left_join(
    pred_values,
    training_df %>% select(cell, class),
    by = "cell"
  ) %>%
    drop_na() %>% # Remove rows with missing predictor or class
    slice_sample(n = min(10000, nrow(pred_values)))
  
  # Separate predictors and class
  predictors_scaled <- scale(training_df_env %>% select(-cell, -class))
  
  # Run PCA
  pca_res <- prcomp(predictors_scaled, center = TRUE, scale. = TRUE)
  
  # Get PCA scores (first 2 PCs)
  pca_scores <- as.data.frame(pca_res$x[, 1:2])
  pca_scores$class <- training_df_env$class
  
  # Get PCA loadings for environmental interpretation
  loadings_df <- as.data.frame(pca_res$rotation[, 1:2])
  loadings_df$varname <- rownames(loadings_df)
  
  loadings_df <- loadings_df %>%
    mutate(magnitude = sqrt(PC1^2 + PC2^2)) %>%
    filter(magnitude > 0.45)  # adjust threshold as needed
  
  # Scale factor for arrows and text
  arrow_scale <- 5
  
  max_val <- max(abs(pca_scores$PC1), abs(pca_scores$PC2))
  
  # PCA plot with biplot arrows
  env_pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = class)) +
    geom_point(alpha = 0.1, size = 1.6) +
    stat_ellipse(aes(group = class, fill = class), 
                 geom = "polygon", alpha = 0.15, color = NA) +
    scale_color_viridis_d(option = "D", direction = -1, guid = 'none') +
    scale_fill_viridis_d(option = "D", direction = -1, guide = 'none') +
    
    # Arrows
    geom_segment(
      data = loadings_df,
      inherit.aes = FALSE,
      aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
      arrow = arrow(length = unit(0.25, "cm")),
      color = "gray20",
      linewidth = 01.5
    ) +
    
    # Variable labels with repel
    geom_text_repel(
      data = loadings_df,
      inherit.aes = FALSE,
      aes(x = PC1 * arrow_scale, y = PC2 * arrow_scale, label = varname),
      size = 10,
      box.padding = 0.4,
      segment.color = "gray50",
      max.overlaps = 20,
      min.segment.length = 0
    ) +
    
    labs(
      x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)"),
    ) +
    
   # coord_equal() +
    theme_minimal(base_size = 18) +
    theme(
      axis.title = element_text(size = 25),
      axis.text = element_text(size = 25),
      panel.grid.minor = element_blank()
    ) 
  
  
  # 10. Return results
  list(
    biplot = biplot,
    spatial_plot = spatial_plot,
    env_pca = env_pca_plot,
    raster = cluster_rast,
    data = training_df
  )
}

## Randall

### Landsat
result_k5 <- run_kmeans_analysis(
  oak_rast =  r_lsc,
  suit_rast = fut3_randall_lsc,
  mask_rast =  r_lsc,
  k = 5,
  pred_rast = r_lsc_pred,
  shape = randall,
  hill_df = hillshade_dfs$randall_lsc,
  scale = FALSE,
  cluster_order = NULL
)

# save map and biplot
ggsave(result_k5$spatial_plot, filename = "3_Figures/Randall_Landsat_map.png", width = 8, height = 8, dpi = 300)
ggsave(result_k5$biplot, filename = "3_Figures/Randall_Landsat_biplot.png", width = 8, height = 6, dpi = 300)
ggsave(result_k5$env_pca, filename = "3_Figures/Randall_Landsat_PCA.png", width = 8, height = 8, dpi = 300)

# save map as categorical raster
writeRaster(result_k5$raster, "1_Data/08_Outputs/Randall_Landsat_ecological_classes.tif", overwrite = TRUE)

# summarize proportion of blue oak woodland landscape covered by each class
class_proportions <- result_k5$data %>%
  group_by(class) %>%
  summarize(proportion = n() / nrow(.)) %>%
  arrange(desc(proportion))


### Sentinel
result_k5_s2 <- run_kmeans_analysis(
  oak_rast = r_sc,
  suit_rast = fut3_randall_s2,
  mask_rast = r_sc,
  k = 5,
  pred_rast = r_sc_pred,
  shape = randall,
  hill_df = hillshade_dfs$randall_s2,
  scale = FALSE,
  cluster_order = c(1,3,2,4,5)
)

# save map and biplot
ggsave(result_k5_s2$spatial_plot, filename = "3_Figures/Randall_Sentinel_map.png", width = 8, height = 8, dpi = 300)
ggsave(result_k5_s2$biplot, filename = "3_Figures/Randall_Sentinel_biplot.png", width = 8, height = 6, dpi = 300)
ggsave(result_k5_s2$env_pca, filename = "3_Figures/Randall_Sentinel_PCA.png", width = 8, height = 8, dpi = 300)

# save map as categorical raster
writeRaster(result_k5_s2$raster, "1_Data/08_Outputs/Randall_Sentinel_ecological_classes.tif", overwrite = TRUE)

# summarize proportion of blue oak woodland landscape covered by each class
class_proportions <- result_k5_s2$data %>%
  group_by(class) %>%
  summarize(proportion = n() / nrow(.)) %>%
  arrange(desc(proportion))


## Dye Creek

#### Landsat
result_k5_dye <- run_kmeans_analysis(
  oak_rast = d_lsc,
  suit_rast = fut3_dye_lsc,
  mask_rast = d_lsc,
  k = 5,
  pred_rast = d_lsc_pred,
  hill_df = hillshade_dfs$dye_lsc,
  shape = dyecreek,
  scale = FALSE,
  cluster_order = c(5,4,3,1,2)
)

ggsave(result_k5_dye$spatial_plot, filename = "3_Figures/DyeCreek_Landsat_map.png", width = 8, height = 8, dpi = 300)
ggsave(result_k5_dye$biplot, filename = "3_Figures/DyeCreek_Landsat_biplot.png", width = 8, height = 6, dpi = 300)
ggsave(result_k5_dye$env_pca, filename = "3_Figures/DyeCreek_Landsat_PCA.png", width = 8, height = 8, dpi = 300)

# save map as categorical raster
writeRaster(result_k5_dye$raster, "1_Data/08_Outputs/DyeCreek_Landsat_ecological_classes.tif", overwrite = TRUE)

# summarize proportion of blue oak woodland landscape covered by each class
class_proportions <- result_k5_dye$data %>%
  group_by(class) %>%
  summarize(proportion = n() / nrow(.)) %>%
  arrange(desc(proportion))


#### Sentinel
result_k5_s2_dye <- run_kmeans_analysis(
  oak_rast = d_sc,
  suit_rast = fut3_dye_s2,
  mask_rast = d_sc,
  k = 5,
  pred_rast = d_sc_pred,
  shape = dyecreek,
  hill_df = hillshade_dfs$dye_s2,
  scale = FALSE,
  cluster_order = c(4,5,2,1,3)
)

# save map and biplot
ggsave(result_k5_s2_dye$spatial_plot, filename = "3_Figures/DyeCreek_Sentinel_map.png", width = 8, height = 8, dpi = 300)
ggsave(result_k5_s2_dye$biplot, filename = "3_Figures/DyeCreek_Sentinel_biplot.png", width = 8, height = 6, dpi = 300)
ggsave(result_k5_s2_dye$env_pca, filename = "3_Figures/DyeCreek_Sentinel_PCA.png", width = 8, height = 8, dpi = 300)

# save map as categorical raster
writeRaster(result_k5_s2_dye$raster, "1_Data/08_Outputs/DyeCreek_Sentinel_ecological_classes.tif", overwrite = TRUE)

# summarize proportion of blue oak woodland landscape covered by each class
class_proportions <- result_k5_s2_dye$data %>%
  group_by(class) %>%
  summarize(proportion = n() / nrow(.)) %>%
  arrange(desc(proportion))

## TMR 1
dye_s2_scaled <- run_kmeans_analysis(
  oak_rast = d_sc2,
  suit_rast = fut3_dye_s2,
  mask_rast = d_sc,
  k = 5,
  pred_rast = d_sc_pred,
  shape = dyecreek,
  hill_df = hillshade_dfs$dye_s2,
  scale = TRUE
)

# TMR 1 unscaled
dye_s2_unscaled <- run_kmeans_analysis(
  oak_rast = d_sc2,
  suit_rast = fut3_dye_s2,
  mask_rast = d_sc,
  k = 5,
  pred_rast = d_sc_pred,
  shape = dyecreek,
  hill_df = hillshade_dfs$dye_s2
)



##%######################################################%##
#                                                          #
####            Elbow plots                             ####
#                                                          #
##%######################################################%##
library(factoextra)

plot_elbow <- function(result_obj, title) {
  df <- data.frame(
    oak_cover = result_obj$data$oak_change,
    suitability = result_obj$data$suitability
  )
  
  # Sample 1000 points without replacement if possible
  n_sample <- min(1000, nrow(df))
  df_sample <- df[sample(nrow(df), n_sample, replace = FALSE), ]
  
  fviz_nbclust(df_sample, kmeans, method = "wss") +
    labs(title = title) +
    theme_minimal() +
    theme(text = element_text(size = 20))
}

p3 <- plot_elbow(result_k5_dye, "a) Dye Creek – Landsat")
p4 <- plot_elbow(result_k5_s2_dye, "b) Dye Creek – Sentinel-2")
p1 <- plot_elbow(result_k5, "c) Randall – Landsat")
p2 <- plot_elbow(result_k5_s2, "d) Randall – Sentinel-2")

library(patchwork)
(p1 | p2) / (p3 | p4)

combined_plot <- (p1 | p2) / (p3 | p4) 

ggsave("3_Figures/elbow_plots_combined.png", combined_plot, width = 12, height = 8, dpi = 300)

##%######################################################%##
#                                                          #
####            Comparing to Oak Recruitment            ####
#                                                          #
##%######################################################%##

### Randall

plots <- vect('1_Data/00_Shapefiles/Randall_Blue_Oak_Recruitment/Randall_Blue_Oak_Recruitment.shp') %>%
  project(randall)

#### Landsat
# Extract ecological class from raster at plot locations
result_k5_r <- rast('1_Data/08_Outputs/Randall_Landsat_ecological_classes.tif')
oak_classes_ls <- terra::extract(result_k5_r, plots, bind = TRUE, search_radius = 120)
oak_classes_ls <- terra::extract(r_lsc_pred, oak_classes_ls, bind = TRUE)
oak_class_df_ls <- as.data.frame(oak_classes_ls)

# Ensure recruitment is an ordered factor
# oak_class_df_ls$class <- factor(
#   oak_class_df_ls$class,
#   levels = c("High", "Low/Medium", "None")
# )

# Define your custom Viridis palette (skipping yellow)
custom_viridis <- viridis::viridis(n = 5, option = "D", direction = -1)

# Plot: proportion of ecological classes per recruitment class
ggplot(oak_class_df_ls, aes(x = class, fill = as.factor(oak_change))) +
  geom_bar(position = "stack") +
  scale_fill_manual(
    name = "Ecological Class",
    values = c("#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")
  ) +
  labs(
    x = "Recruitment Class",
    y = "Number of Plots"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.direction = 'horizontal',
    legend.position = "bottom"
  ) 

# save
ggsave("3_Figures/Randall_Landsat_Ecological_classes_recruitment.png", width = 10, height = 8, dpi = 300)

# Environmental characteristics of oak recruitment plots
ggplot(oak_class_df_ls, aes(x = tri, y = count, color = Ranch)) + geom_point(size = 3)
ggplot(oak_class_df_ls, aes(x = elevation, y = count, color = Ranch)) + geom_point(size = 3)
ggplot(oak_class_df_ls, aes(x = cwd, y = count, color = Ranch)) + geom_point(size = 3)
ggplot(oak_class_df_ls, aes(x = tmn, y = count, color = Ranch)) + geom_point(size = 3)
ggplot(oak_class_df_ls, aes(x = tpi, y = count, color = Ranch)) + geom_point(size = 3)
ggplot(oak_class_df_ls, aes(x = pct_clay, y = count, color = Ranch)) + geom_point(size = 3)
ggplot(oak_class_df_ls, aes(x = depth, y = count, color = Ranch)) + geom_point(size = 3)
ggplot(oak_class_df_ls, aes(x = ppt_djf, y = count, color = Ranch)) + geom_point(size = 3)
ggplot(oak_class_df_ls, aes(x = awc, y = count, color = Ranch)) + geom_point(size = 3)
ggplot(oak_class_df_ls, aes(x = oak_change, y = count, color = Ranch)) + geom_point(size = 3)

oak_class_df_ls$outlier <- with(oak_class_df_ls,
                                oak_change %in% c(1,2) &
                                  count <= quantile(count, 0.5, na.rm = TRUE))

oak_class_df_ls %>%
  group_by(outlier) %>%
  summarise(across(c(awc, cwd, depth, pct_clay, ph, ppt_djf, tmn, elevation, tri, tpi),
                   list(mean = mean, sd = sd), na.rm = TRUE))

ggplot(oak_class_df_ls %>% filter(oak_change %in% c(1,2)),
       aes(x = outlier, y = ph)) +
  geom_point() +
  labs(x = "Outlier (High class + Low recruit)", y = "")

# Figure caption: 
# Proportional distribution of ecological classes across recruitment categories. The bar plot shows the relative abundance of each ecological class—derived from k-means clustering of historical Landsat-based oak resilience and future blue oak habitat suitability—within field-assigned recruitment classes ("High", "Medium", "Low", "None"). Each ecological class represents a unique combination of predicted oak cover resilience and future habitat suitability.

#### Sentinel
# Extract ecological class from raster at plot locations
oak_classes_s2 <- terra::extract(result_k5_s2$raster, plots, bind = TRUE, search_radius = 120)
oak_class_df_s2 <- as.data.frame(oak_classes_s2)
# Ensure recruitment is an ordered factor
# oak_class_df_s2$class<- factor(
#   oak_class_df_s2$class,
#   levels = c("High", "Medium", "Low", "None")
# )

# Plot: proportion of ecological classes per recruitment class
ggplot(oak_class_df_s2, aes(x = class, fill = as.factor(oak_change))) +
  geom_bar(position = "stack") +
  scale_fill_manual(
    name = "Ecological Class",
    values = c("#FDE725FF", "#5DC863FF", "#21908CFF","#3B528BFF", "#440154FF")
  ) +
  labs(
    x = "Recruitment Class",
    y = "Number of Plots"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.direction = 'horizontal',
    legend.position = "bottom"
  )

# save
ggsave("3_Figures/Randall_Sentinel_Ecological_classes_recruitment.png", width = 10, height = 8, dpi = 300)

# Map of Randall Preserve, blue oak range map, and plots
# raster = qd_ran_r
# shapefile = randall
# plots = plots (color = class)


# Convert raster to dataframe
qd_ran_df <- as.data.frame(qd_ran_r, xy = TRUE, na.rm = TRUE)
names(qd_ran_df)[3] <- "blue_oak"
qd_ran_df <- qd_ran_df %>% filter(blue_oak == 1)
qd_ran_df$blue_oak <- as.factor(qd_ran_df$blue_oak)

# Convert plots to sf
plots_sf <- st_as_sf(plots, coords = c("x", "y"), crs = crs(qd_ran_r))

# Convert preserve boundary to sf if it's not already
randall_sf <- st_as_sf(randall)

# Build map
randall_map <- ggplot() +
  geom_raster(data = hillshade_dfs$randall_s2, aes(x = x, y = y, fill = hillshade)) +
  scale_fill_gradient(low = "white", high = "gray60", guide = "none") +  # soft shading
  new_scale_fill() +  # Allows for multiple fill scales
  # Raster layer
  geom_raster(data = qd_ran_df, aes(x = x, y = y, fill = blue_oak), alpha = .6) +
  scale_fill_manual(
    values = c("#3B528BFF"),
    breaks = c(1),
    guide = 'none'
    
  ) +
  
  # Preserve boundary
  geom_sf(data = randall_sf, fill = NA, color = "black", size = 1) +
  
  # Plot locations
  geom_sf(data = plots_sf, aes(color = class), size = 6) +
  # scale_color_viridis_d(option = "C", direction = -1, name = 'Recruitment Class') +
  scale_color_manual(values = c(
    "High" = "#00CED1",
    "Low/Medium" = "#FF8C00",
    "None" = "#C71585"
  ), name = 'Recruitment Class:') +
  geom_sf(
    data = st_as_sf(randall),
    fill = NA,
    color = "black",
    lwd = 1
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
    height = unit(2, "cm"),
    width = unit(1.5, "cm")
  ) +
  coord_sf() +
  theme_minimal() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
   # legend.background = element_rect(fill = "white", color = "black"),
    legend.position = 'bottom',
    legend.direction = 'horizontal',
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30)
  )

# Show the map
randall_map

# save
ggsave("3_Figures/Randall_Preserve_Map.png", plot = randall_map, width = 12, height = 8, dpi = 300)

##%######################################################%##
#                                                          #
####           Nice plots of the TMR3 results           ####
#                                                          #
##%######################################################%##

# list of rasters r_sc, r_lsc, d_sc, d_lsc

library(ggnewscale)
library(ggspatial)
library(viridis)

plot_tmr_raster <- function(
    rast,                  # SpatRaster or RasterLayer
    hill_df,               # Hillshade data frame
    boundary_sf,           # Preserve boundary (sf)
    fill_label = "Value",  # Legend title
    fill_limits = NULL,    # Optional limits for fill color scale
    out_file = NULL        # Optional output file path
) {
  # Convert raster to dataframe
  rast_df <- as.data.frame(rast, xy = TRUE, na.rm = TRUE)
  names(rast_df)[3] <- "value"
  
  # Build plot
  tmr_map <- ggplot() +
    # Hillshade
    geom_raster(data = hill_df, aes(x = x, y = y, fill = hillshade)) +
    scale_fill_gradient(low = "white", high = "gray60", guide = "none") +
    new_scale_fill() +
    
    # Main raster
    geom_raster(data = rast_df, aes(x = x, y = y, fill = value), alpha = 0.9) +
    scale_fill_viridis(
      option = "B",
      name = fill_label,
      limits = fill_limits,
      oob = scales::squish,
      guide = guide_colorbar(
        barwidth = unit(8, "cm"),
        barheight = unit(0.6, "cm"),
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    
    # Preserve boundary
    geom_sf(data = boundary_sf, fill = NA, color = "black", size = 1) +
    
    annotation_scale(
      location = "bl",
      width_hint = 0.2,
      height = unit(0.6, "cm"),
      text_cex = 1.5
    ) +
    annotation_north_arrow(
      location = "br",
      which_north = "true",
      height = unit(2, "cm"),
      width = unit(1.5, "cm")
    ) +
    coord_sf() +
    theme_minimal() +
    labs(title = NULL, x = NULL, y = NULL) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
     # legend.background = element_rect(fill = "white", color = "black"),
      legend.position = 'bottom',
      legend.direction = 'horizontal',
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 24)
    )
  
  # Display the map
  print(tmr_map)
  
  # Save if requested
  if (!is.null(out_file)) {
    ggsave(out_file, plot = tmr_map, width = 12, height = 8, dpi = 300)
  }
}

# r_sc: Randall - Sentinel Composite
plot_tmr_raster(
  rast = r_sc,
  hill_df = hillshade_dfs$randall_s2,
  boundary_sf = randall_sf,
  fill_limits = c(-.2,.2),
  fill_label = "TMR3 Score",
  out_file = "3_Figures/Randall_TMR3_Sentinel.png"
)

# r_lsc: Randall - Landsat Composite
plot_tmr_raster(
  rast = r_lsc,
  hill_df = hillshade_dfs$randall_ls,
  boundary_sf = randall_sf,
  fill_limits = c(-.2,.2),
  fill_label = "TMR3 Score",
  out_file = "3_Figures/Randall_TMR3_Landsat.png"
)

# d_sc: Dye Creek - Sentinel Composite
plot_tmr_raster(
  rast = d_sc,
  hill_df = hillshade_dfs$dye_s2,
  boundary_sf = st_as_sf(dyecreek),
  fill_label = "TMR3 Score",
  fill_limits = c(-.2, .2),
  out_file = "3_Figures/Dye_TMR3_Sentinel.png"
)

# d_lsc: Dye Creek - Landsat Composite
plot_tmr_raster(
  rast = d_lsc,
  hill_df = hillshade_dfs$dye_lsc,
  boundary_sf = st_as_sf(dyecreek),
  fill_label = "TMR3 Score",
  fill_limits = c(-.2, .2),
  out_file = "3_Figures/Dye_TMR3_Landsat.png"
)

##%######################################################%##
#                                                          #
####                  Fire perimeters                   ####
#                                                          #
##%######################################################%##

plot_fire_map <- function(hillshade_df, preserve_sf, fire_sf, title = NULL) {
  ggplot() +
    # Hillshade background
    geom_raster(data = hillshade_df, aes(x = x, y = y, fill = hillshade)) +
    scale_fill_gradient(low = "white", high = "gray60", guide = "none") +
    new_scale_fill() +
    
    geom_sf(data = fire_sf, aes(fill = fire_year), color = "black", size = 0.4, alpha = 0.6) +
    scale_fill_viridis_d(
      option = "C", direction = -1,
      name = "Fire Year",
      guide = guide_legend(ncol = 4)  # Adjust number of columns for clarity
    ) +
    
    # Preserve boundary
    geom_sf(data = preserve_sf, fill = NA, color = "black", size = 1) +
    
    # North arrow and scale
    annotation_scale(
      location = "bl", width_hint = 0.2, height = unit(0.6, "cm"), text_cex = 1.5
    ) +
    annotation_north_arrow(
      location = "br", which_north = "true",
      height = unit(2, "cm"), width = unit(1.5, "cm")
    ) +
    
    coord_sf() +
    theme_minimal() +
    labs(title = title, x = NULL, y = NULL) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = 'bottom',
      legend.direction = 'horizontal',
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.background = element_rect(fill = "white", color = "black"),
      plot.title = element_text(size = 20, hjust = 0.5)
    )
}

# Load fire perimeters

### Randall
r_fire <- terra::vect('1_Data/00_Shapefiles/Fire/randall_fires.shp') %>%
  project(randall)

# fires after 1930
r_fire <- fire[fire$YEAR_ >= 1930, ]

r_fire_rx <- terra::vect('1_Data/00_Shapefiles/Fire/randall_fires_rx.shp') %>%
  project(randall)

r_fire_rx <- fire_rx[fire_rx$YEAR_ >= 1930, ]

### Dye Creek
d_fire <- terra::vect('1_Data/00_Shapefiles/Fire/dyecreek_fires.shp') %>%
  project(dyecreek)

# fires after 1930
d_fire <- d_fire[fire$YEAR_ >= 1930, ]

d_fire_rx <- terra::vect('1_Data/00_Shapefiles/Fire/dyecreek_fires_rx.shp') %>%
  project(dyecreek)

d_fire_rx <- d_fire_rx[fire_rx$YEAR_ >= 1930, ]

park_fire <- terra::vect('1_Data/00_Shapefiles/Fire/Park_Fire/Park_Fire.shp') %>%
  project(dyecreek) %>% crop(dyecreek)

park_fire$YEAR_ <- 2024


r_fire$fire_year <- as.factor(r_fire$YEAR_)
r_fire_rx$fire_year <- as.factor(r_fire_rx$YEAR_)

d_fire$fire_year <- as.factor(d_fire$YEAR_)
d_fire_rx$fire_year <- as.factor(d_fire_rx$YEAR_)
park_fire$fire_year <- as.factor(park_fire$YEAR_)


# Randall Wildfire Map
randall_fire_map <- plot_fire_map(
  hillshade_df = hillshade_dfs$randall_lsc,
  preserve_sf = randall_sf,
  fire_sf = st_as_sf(r_fire),
  title = "Randall Preserve Wildfires (Post-1930)"
)

# save
ggsave("3_Figures/Randall_Preserve_Wildfire_Map.png", plot = randall_fire_map, width = 12, height = 8, dpi = 300)

# Randall Rx Map
randall_rx_map <- plot_fire_map(
  hillshade_df = hillshade_dfs$randall_lsc,
  preserve_sf = randall_sf,
  fire_sf = st_as_sf(r_fire_rx),
  title = "Randall Preserve Prescribed Fires (Post-1930)"
)

# save
ggsave("3_Figures/Randall_Preserve_Rx_Map.png", plot = randall_rx_map, width = 12, height = 8, dpi = 300)

# Dye Creek Wildfire Map (including Park Fire)
d_fire_sf <- st_as_sf(d_fire)
park_fire_sf <- st_as_sf(park_fire)

# Select only YEAR_ and geometry (assuming the year column is called "YEAR_")
d_fire_sf <- d_fire_sf[, c("fire_year", "geometry")]
park_fire_sf <- park_fire_sf[, c("fire_year", "geometry")]

# Now combine
dye_fire_combined <- rbind(d_fire_sf, park_fire_sf)
dye_fire_map <- plot_fire_map(
  hillshade_df = hillshade_dfs$dye_lsc,
  preserve_sf = st_as_sf(dyecreek),
  fire_sf = dye_fire_combined,
  title = "Dye Creek Wildfires (Post-1930)"
)

# save
ggsave("3_Figures/DyeCreek_Wildfire_Map.png", plot = dye_fire_map, width = 12, height = 8, dpi = 300)

# Dye Creek Rx Map
dye_rx_map <- plot_fire_map(
  hillshade_df = hillshade_dfs$dye_lsc,
  preserve_sf = st_as_sf(dyecreek),
  fire_sf = st_as_sf(d_fire_rx),
  title = "Dye Creek Prescribed Fires (Post-1930)"
)

# save
ggsave("3_Figures/DyeCreek_Rx_Map.png", plot = dye_rx_map, width = 12, height = 8, dpi = 300)


##%######################################################%##
#                                                          #
####        Future Blue Oak Habitat Suitability         ####
#                                                          #
##%######################################################%##


fut3_dye_df <- as.data.frame(fut3_dye_s2 %>% mask(dyecreek), xy = TRUE, na.rm = TRUE)
colnames(fut3_dye_df)[3] <- "suitability"

fut3_randall_df <- as.data.frame(fut3_randall_s2 %>% mask(randall), xy = TRUE, na.rm = TRUE)
colnames(fut3_randall_df)[3] <- "suitability"

plot_suitability_map <- function(hillshade_df, suitability_df, preserve_sf, title = NULL) {
  ggplot() +
    # Hillshade background
    geom_raster(data = hillshade_df, aes(x = x, y = y, fill = hillshade)) +
    scale_fill_gradient(low = "white", high = "gray60", guide = "none") +
    new_scale_fill() +
    
    # Suitability raster
    geom_raster(data = suitability_df, aes(x = x, y = y, fill = suitability), alpha = .6) +
    scale_fill_viridis_c(
      option = "C", direction = 1,
      name = "Suitability", limits = c(0, 1),
      guide = guide_colorbar(
        barwidth = 20, barheight = 1.5,  # wider and thicker
        title.position = "top", title.hjust = 0.5,
        label.position = "bottom", label.hjust = 0.5
      )
    ) +
    
    # Preserve boundary
    geom_sf(data = preserve_sf, fill = NA, color = "black", size = 1) +
    
    # North arrow and scale bar
    annotation_scale(location = "bl", width_hint = 0.2, height = unit(0.6, "cm"), text_cex = 1.5) +
    annotation_north_arrow(location = "br", which_north = "true", height = unit(2, "cm"), width = unit(1.5, "cm")) +
    
    coord_sf() +
    theme_minimal() +
    labs(title = title, x = NULL, y = NULL) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = 'bottom',
      legend.direction = 'horizontal',
      legend.title = element_text(size = 22, margin = margin(b = 10)),
      legend.text = element_text(size = 22),
      plot.title = element_text(size = 20, hjust = 0)
    )
}



p1 <- plot_suitability_map(
  suitability_df = fut3_dye_df,
  hillshade_df = hillshade_dfs$dye_lsc,
  preserve_sf = st_as_sf(dyecreek),
  title = "c) Blue oak habitat suitability (2055) at Dye Creek Preserve"
)

p2 <- plot_suitability_map(
  suitability_df = fut3_randall_df,
  hillshade_df = hillshade_dfs$randall_lsc,
  preserve_sf = st_as_sf(randall),
  title = "d) Blue oak habitat suitability (2055) at Randall Preserve"
)

library(patchwork)

# Combine with shared legend below
final_plot <- (p2 / p1) +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom',
        panel.background = element_rect(fill = "white"))

# save
ggsave("3_Figures/Future_Blue_Oak_Habitat_Suitability.png", final_plot, width = 15, height = 14, dpi = 300)


# Current suitability maps

curr_dye_df <- as.data.frame(curr_dye_s2 %>% mask(dyecreek), xy = TRUE, na.rm = TRUE)
colnames(curr_dye_df )[3] <- "suitability"

p3 <- plot_suitability_map(
  suitability_df = curr_dye_df ,
  hillshade_df = hillshade_dfs$dye_lsc,
  preserve_sf = st_as_sf(dyecreek),
  title = "a) Blue oak habitat suitability (1995) at Dye Creek Preserve"
)

cur_randall_df <- as.data.frame(cur_randall_s2 %>% mask(randall), xy = TRUE, na.rm = TRUE)
colnames(cur_randall_df)[3] <- "suitability"

p4 <- plot_suitability_map(
  suitability_df = cur_randall_df ,
  hillshade_df =  hillshade_dfs$randall_lsc,
  preserve_sf = st_as_sf(randall),
  title = "b) Blue oak habitat suitability (1995) at Randall Preserve"
)

# Combine with shared legend below
final_plot <- (p3 / p4) +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom',
        panel.background = element_rect(fill = "white"))

# save
ggsave("3_Figures/Current_Blue_Oak_Habitat_Suitability.png", final_plot, width = 15, height = 14, dpi = 300)

##%######################################################%##
#                                                          #
####                 Elevation figures                  ####
#                                                          #
##%######################################################%##

library(ggplot2)
library(viridis)
library(ggspatial)
library(sf)
library(ggnewscale)  # for new_scale_fill()

# Function to plot elevation
plot_elevation_map <- function(hillshade_df, elevation_raster, preserve_sf, title = NULL) {
  
  # Convert raster to data frame if not already
  if (inherits(elevation_raster, "SpatRaster")) {
    elevation_df <- as.data.frame(elevation_raster, xy = TRUE)
    names(elevation_df) <- c("x", "y", "elevation")
  } else {
    elevation_df <- elevation_raster
  }
  
  ggplot() +
    # Hillshade background
    geom_raster(data = hillshade_df, aes(x = x, y = y, fill = hillshade)) +
    scale_fill_gradient(low = "white", high = "gray60", guide = "none") +
    new_scale_fill() +
    
    # Elevation raster
    geom_raster(data = elevation_df, aes(x = x, y = y, fill = elevation), alpha = 0.7) +
    scale_fill_viridis_c(
      option = "D", direction = 1,
      name = "Elevation (m)",
      guide = guide_colorbar(
        barwidth = 20, barheight = 1.5,
        title.position = "top", title.hjust = 0.5,
        label.position = "bottom", label.hjust = 0.5
      )
    ) +
    
    # Preserve boundary
    geom_sf(data = preserve_sf, fill = NA, color = "black", size = 1) +
    
    coord_sf() +
    theme_minimal() +
    labs(title = title, x = NULL, y = NULL) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 22, margin = margin(b = 10)),
      legend.text = element_text(size = 22),
      plot.title = element_text(size = 20, hjust = 0)
    )
}

# Example usage:
# Assuming hillshade_df, preserve_sf, and terrain_metrics_s2$elevation are available
r_elev <- plot_elevation_map(hillshade_dfs$randall_lsc, terrain_metrics_s2$elevation %>% mask(randall), st_as_sf(randall), title = "")
dye_elev <- plot_elevation_map(hillshade_dfs$dye_lsc, dye_terrain_s2$elevation %>% mask(dyecreek), st_as_sf(dyecreek), title = "")

ggsave('3_Figures/Randall_elevation.png', r_elev, width = 8, height = 8, dpi = 500)
ggsave('3_Figures/DyeCreek_elevation.png', dye_elev, width = 10, height = 8, dpi = 500)

##%######################################################%##
#                                                          #
####     Environmental Drivers of Oak Cover Change      ####
#                                                          #
##%######################################################%##

# Prepare data for Randall

### Terrain

##### Randall
#### Sentinel-2 (10m)
terrain_pred <- terrain_metrics %>% resample(r_sc)
terrain_pred <- subset(terrain_pred, !names(terrain_pred) %in% c("aspect", "flowdir", "slope"))

result_s2_terrain <- analyze_oak_change(r_sc,
                                        terrain_pred,
                                        palette = "magma",
                                        top_n = 5)


#### Landsat (30 m)
terrain_pred_lsc <- terrain_metrics %>% resample(r_lsc)
terrain_pred_lsc <- subset(terrain_pred_lsc, !names(terrain_pred_lsc) %in% c("aspect", "flowdir", "slope"))

result_ls_terrain <- analyze_oak_change(r_lsc,
                                        terrain_pred_lsc,
                                        palette = "magma",
                                        top_n = 5)

### Environmental Drivers (270m)



#### Sentinel-2
r_sc_resamp <- resample(r_sc, randall_env2, method = "bilinear")

randall_env2_subset <- subset(randall_env2, !names(randall_env2) %in% c("ppt_jja", "aet"))

result_s2_env <- analyze_oak_change(r_sc_resamp,
                                    randall_env2_subset,
                                    palette = "magma",
                                    top_n = 8)


#### Landsat
r_lsc_resamp <- resample(r_lsc, randall_env2, method = "bilinear")

result_ls_env <- analyze_oak_change(r_lsc_resamp,
                                    randall_env2_subset,
                                    palette = "magma",
                                    top_n = 8)


#### Dye Creek

### Terrain
#### Sentinel-2 (10m)
d_terrain_pred <- dye_creek_terrain %>% resample(d_sc)
d_terrain_pred <- subset(d_terrain_pred, !names(d_terrain_pred) %in% c("aspect", "flowdir", "slope"))
result_s2_terrain_d <- analyze_oak_change(d_sc,
                                          d_terrain_pred,
                                          palette = "magma",
                                          top_n = 5)
#### Landsat
d_terrain_pred_lsc <- dye_creek_terrain %>% resample(d_lsc)
d_terrain_pred_lsc <- subset(d_terrain_pred_lsc, !names(d_terrain_pred_lsc) %in% c("aspect", "flowdir", "slope"))
result_ls_terrain_d <- analyze_oak_change(d_lsc,
                                          d_terrain_pred_lsc,
                                          palette = "magma",
                                          top_n = 5)


### Environment
#### Sentinel-2
d_sc_resamp <- resample(d_sc, dyecreek_env2, method = "bilinear")
dyecreek_env2_subset <- subset(dyecreek_env2, !names(dyecreek_env2) %in% c("ppt_jja", "aet"))
result_s2_env_d <- analyze_oak_change(d_sc_resamp,
                                      dyecreek_env2_subset,
                                      palette = "magma",
                                      top_n = 8)

#### Landsat
d_lsc_resamp <- resample(d_lsc, dyecreek_env2, method = "bilinear")
result_ls_env_d <- analyze_oak_change(d_lsc_resamp,
                                      dyecreek_env2_subset,
                                      palette = "magma",
                                      top_n = 8)
