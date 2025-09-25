##%######################################################%##
#                                                          #
####    Summarizing vegetation trends for TMR scores    ####
####           at both Dye Creek and Randall            ####
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

##%######################################################%##
#                                                          #
####             Define a tmr plot function             ####
#                                                          #
##%######################################################%##

compute_pc_phenology <- function(veg_rast = NULL, veg_df = NULL,
                                 pcs_rast = NULL,
                                 n_sample = 1000,
                                 out_dir = NULL,
                                 save_plots = TRUE,
                                 save_data = TRUE) {
  
  library(terra); library(dplyr); library(tidyr); library(zoo); library(ggplot2); library(lubridate)
  
  # --- 1. Sample points from PC raster ---
  if (!is.null(pcs_rast)) {
    set.seed(2024)
    sample_points <- spatSample(pcs_rast, size = n_sample, method = "random", as.points = TRUE, na.rm = TRUE)
    
    # Extract PC values
    pc_vals <- terra::extract(pcs_rast, sample_points)
    pc_df <- cbind(pixel_id = seq_len(nrow(pc_vals)), pc_vals) %>% dplyr::select(-ID)
    colnames(pc_df)[-1] <- paste0("PC", 1:(ncol(pc_df)-1))
  } else {
    stop("You must provide pcs_rast")
  }
  
  # --- 2. Extract vegetation fractions ---
  if (!is.null(veg_rast)) {
    veg_vals <- terra::extract(veg_rast, sample_points)
    colnames(veg_vals) <- gsub("^X", "", colnames(veg_vals))
    colnames(veg_vals)[1] <- "pixel_id"
    veg_df <- veg_vals %>%
      pivot_longer(cols = -pixel_id, names_to = "date", values_to = "veg_fraction") %>%
      mutate(date = as.numeric(gsub("^X", "", date))) %>%
      filter(!is.na(date), !is.na(veg_fraction))
  } else if (is.null(veg_df)) {
    stop("You must provide either veg_rast or veg_df")
  }
  
  plot_list <- list()
  summary_list <- list()
  annual_list <- list()
  
  # --- 3. Loop over PCs ---
  for (i in 1:(ncol(pc_df)-1)) {
    pc_col <- paste0("PC", i)
    
    # High/Low thresholds
    low_thr <- quantile(pc_df[[pc_col]], 0.05, na.rm = TRUE)
    high_thr <- quantile(pc_df[[pc_col]], 0.95, na.rm = TRUE)
    
    # Assign groups
    groups <- pc_df %>%
      select(pixel_id, !!sym(pc_col)) %>%
      mutate(group = case_when(
        !!sym(pc_col) <= low_thr ~ "Low_5%",
        !!sym(pc_col) >= high_thr ~ "High_95%",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(group))
    
    # Merge with veg_df
    veg_grp <- veg_df %>% inner_join(groups, by = "pixel_id")
    
    # Aggregate mean vegetation fraction per group/date
    agg_ts <- veg_grp %>%
      group_by(group, date) %>%
      summarise(mean_veg = mean(veg_fraction, na.rm = TRUE), .groups = "drop") %>%
      mutate(year = floor(date),
             DOY = floor((date - year) * 365.25) + 1)  # fractional year -> day of year
    
    pheno_annual <- agg_ts %>%
      group_by(group, year) %>%
      arrange(DOY) %>%
      mutate(
        roll_mean = zoo::rollmean(mean_veg, 3, fill = NA),
        slope = c(NA, diff(roll_mean))
      ) %>%
      reframe(
        GreenUp_DOY = if(any(!is.na(slope))) DOY[which.max(slope)] else NA_integer_,
        Senescence_DOY = if(any(!is.na(slope))) {
          idx_gu <- which.max(slope)
          if(length(idx_gu) == 0 || idx_gu >= length(slope)) {
            NA_integer_
          } else {
            idx_sen <- which.min(slope[(idx_gu+1):length(slope)]) + idx_gu
            DOY[idx_sen]
          }
        } else NA_integer_,
        LOS = ifelse(!is.na(GreenUp_DOY) & !is.na(Senescence_DOY), 
                     Senescence_DOY - GreenUp_DOY, NA_integer_),
        Amplitude = if(any(!is.na(mean_veg))) max(mean_veg, na.rm = TRUE) - min(mean_veg, na.rm = TRUE) else NA_real_,
        Integral = if(any(!is.na(mean_veg))) sum(mean_veg, na.rm = TRUE) else NA_real_
      )
    
    # Average metrics across years
    pheno_annual_summary <- pheno_annual %>%
      group_by(group) %>%
      summarise(
        GreenUp_DOY_mean = mean(GreenUp_DOY, na.rm = TRUE),
        Senescence_DOY_mean = mean(Senescence_DOY, na.rm = TRUE),
        LOS_mean = mean(LOS, na.rm = TRUE),
        Amplitude_mean = mean(Amplitude, na.rm = TRUE),
        Integral_mean = mean(Integral, na.rm = TRUE),
        .groups = "drop"
      )
    
    # --- Full time-series metrics ---
    pheno_full <- agg_ts %>%
      mutate(year = floor(date)) %>%
      group_by(group, year) %>%
      summarise(mean_veg_year = mean(mean_veg, na.rm = TRUE), .groups = "drop") %>%
      arrange(year) %>%
      group_by(group) %>%
      summarise(
        Max = max(mean_veg_year, na.rm = TRUE),
        Min = min(mean_veg_year, na.rm = TRUE),
        Amplitude = Max - Min,
        LinearSlope = coef(lm(mean_veg_year ~ year))[2],
        PercentChange = {
          n_years <- n()
          first_avg <- mean(mean_veg_year[1:min(2, n_years)], na.rm = TRUE)
          last_avg  <- mean(mean_veg_year[(n_years-1):n_years], na.rm = TRUE)
          100 * (last_avg - first_avg) / first_avg
        },
        .groups = "drop"
      )
    
    # Combine metrics
    pheno_combined <- pheno_full %>%
      left_join(pheno_annual_summary, by = "group") %>%
      mutate(PC = pc_col)
    
    summary_list[[pc_col]] <- pheno_combined
    annual_list[[pc_col]] <- pheno_annual
    
    # --- Plot ---
    p <- ggplot(agg_ts, aes(x = date, y = mean_veg, color = group)) +
      geom_line(size = 1.5) +
      geom_smooth(se = FALSE, size = 1.7) +
      labs(title = paste0(pc_col, " High vs Low TMR"),
           x = "Year", y = "Mean Vegetation Fraction", color = "Group") +
      scale_color_manual(values = c("Low_5%" = "#000004FF", "High_95%" = "gold")) +
      theme_minimal(base_size = 30) +
      annotate("rect", xmin = 1987, xmax = 1992, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.25) +
      annotate("rect", xmin = 2011.8, xmax = 2017, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.25) +
      annotate("rect", xmin = 2020, xmax = 2022, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.25) +
      annotate("rect", xmin = 2007, xmax = 2009, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.25)
    
    plot_list[[pc_col]] <- p
    
    # Optionally save
    if (!is.null(out_dir)) {
      if (save_plots) ggsave(file.path(out_dir, paste0(pc_col, "_HighLow.png")), p, width = 12, height = 8, dpi = 600)
      if (save_data) write.csv(agg_ts, file.path(out_dir, paste0(pc_col, "_HighLow_TS.csv")), row.names = FALSE)
    }
  }
  
  summary_df <- bind_rows(summary_list)
  annual_df <- bind_rows(annual_list)
  
  return(list(
    plots = plot_list,
    summary = summary_df,
    annual = annual_df
  ))
}



##%######################################################%##
#                                                          #
####                    Data inputs                     ####
#                                                          #
##%######################################################%##

# Helper function to project, rasterize, and mask
rasterize_and_mask <- function(raster_path, shapefile, shapefile_name) {
  r <- rast(raster_path)[[1:3]]
  shp_proj <- project(shapefile, r)
  shp_proj$presence <- 1
  r_masked <- mask(r, rasterize(shp_proj, r, field = "presence", background = NA, touches = TRUE))
  message(paste("✓", shapefile_name, "-", basename(raster_path), "processed"))
  return(r_masked)
}


# Blue oak range maps
qd_randall <- vect('1_Data/00_Shapefiles/blue_oak_range/BlueOakRandall.shp')
qd_dyecreek <- vect('1_Data/00_Shapefiles/blue_oak_range/BlueOakDyeCreek.shp')

# Randall 

## Landsat

rls_veg <- rast('1_Data/Landsat/05_Randall_temporal_stack/Randall_Landsat_1982_2025')
rls_tmr <- rasterize_and_mask(
  raster_path = '1_Data/Landsat/06_Randall_tem_unmixed/Randall_LS_TMR_PCA',
  shapefile = qd_randall,
  shapefile_name = "Randall Landsat"
)

## Sentinel-2

rs2_veg <- rast('1_Data/Sentinel_2/05_Randall_temporal_stack/Randall_S2_stack')
rs2_tmr <- rasterize_and_mask(raster_path = '1_Data/Sentinel_2/06_Randall_tem_unmixed/Randall_S2_TMR_PCA',
                              shapefile = qd_randall,
                              shapefile_name = 'Randall Sentinel')


# Dye Creek

## Landsat

dls_veg <- rast('1_Data/Landsat/05_DyeCreek_temporal_stack/DyeCreek_LS_1982_2025_veg_stack')
dls_tmr <- rasterize_and_mask(raster_path = '1_Data/Landsat/06_DyeCreek_tem_unmixed/full/DyeCreek_LS_1982_2025_TMR_PCA',
                              shapefile = qd_dyecreek,
                              shapefile_name = 'Dye Creek Landsat')


## Sentinel-2

ds2_veg <- rast('1_Data/Sentinel_2/05_DyeCreek_temporal_stack/DyeCreek_S2_stack')
ds2_tmr <- rasterize_and_mask(raster_path = '1_Data/Sentinel_2/06_DyeCreek_tem_unmixed/full/DyeCreek_S2_TMR_full_PCA',
                              shapefile = qd_dyecreek,
                              shapefile_name = 'Dye Creek Sentinel')




##%######################################################%##
#                                                          #
####        Apply vegetation trend plot function        ####
#                                                          #
##%######################################################%##


# --- Randall Landsat ---
res_rls <- compute_pc_phenology(
  veg_rast = rls_veg,
  pcs_rast = rls_tmr[[1:3]],
  n_sample = 1000,
  save_plots = F
)

# Save plots
ggsave(
  res_rls$plots[[3]],
  filename = '3_Figures/Randall_TMR3_temporal_trend_Landsat.png',
  height = 7,
  width = 12,
  dpi = 500
)

# --- Randall Sentinel-2 ---
res_rs2 <- compute_pc_phenology(
  veg_rast = rs2_veg,
  pcs = rs2_tmr[[1:3]],
  n_sample = 1000,
  save_plots = F
)

# Save plots
ggsave(
  res_rs2$plots[[3]],
  filename = '3_Figures/Randall_TMR3_temporal_trend_Sentinel.png',
  height = 7,
  width = 12,
  dpi = 500
)

# --- Dye Creek Landsat ---
res_dls <- compute_pc_phenology(
  veg_rast = dls_veg,
  pcs = dls_tmr[[1:3]],
  n_sample = 1000,
  save_plots = F
)

# Save plots
ggsave(
  res_dls$plots[[3]],
  filename = '3_Figures/DyeCreek_TMR3_temporal_trend_Landsat.png',
  height = 7,
  width = 12,
  dpi = 500
)

# --- Dye Creek Sentinel-2 ---
res_ds2 <- compute_pc_phenology(
  veg_rast = ds2_veg,
  pcs = ds2_tmr[[1:3]],
  n_sample = 1000,
  save_plots = F
)

# Save plots
ggsave(
  res_ds2$plots[[3]],
  filename = '3_Figures/DyeCreek_TMR3_temporal_trend_Sentinel.png',
  height = 7,
  width = 12,
  dpi = 500
)
