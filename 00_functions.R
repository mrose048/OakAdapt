##%######################################################%##
#                                                          #
####  Functions for oak temporal mixture model paper/s  ####
#                                                          #
##%######################################################%##


# Function to stretch contrast (adjusts values between 2nd and 98th percentile)
stretch_contrast <- function(x) {
  q <- quantile(x, probs = c(0.02, 0.98), na.rm = TRUE)  # Get min and max cutoffs
  x <- (x - q[1]) / (q[2] - q[1])  # Scale between 0 and 1
  x[x < 0] <- 0  # Clip values below 0
  x[x > 1] <- 1  # Clip values above 1
  return(x)
}


# Function to plot a raster with custom fill and color palette
plot_map <- function(raster_input, boundary, fill_var, palette, legend_title, legend_position = c(0.1, 0.15)) {
  
  # Convert raster to dataframe for ggplot
  raster_df <- as.data.frame(raster_input, xy = T)
  colnames(raster_df) <- c("x", "y", "value")  # Rename value column
  
  # Convert value column to a factor for categorical mapping
  raster_df$value <- as.factor(raster_df$value)
  
  # Get bounding box of the boundary (this will define the map limits)
  bbox <- st_bbox(boundary)  # Get the bounding box for the boundary
  
  # Ensure only mapped categories get colors
  used_colors <- palette[names(palette) %in% levels(raster_df$value)]
  
  # Create ggplot map
  # Create ggplot map
  p <- ggplot() +
    geom_raster(data = raster_df, aes(x = x, y = y, fill = value)) +
    scale_fill_manual(values = used_colors, name = legend_title) +
    geom_sf(
      data = st_as_sf(boundary),
      fill = NA,
      color = "black",
      lwd = 1
    ) +
    
    # Scale Bar
    annotation_scale(
      location = "bl",
      width_hint = 0.2,
      height = unit(0.6, "cm"),
      text_cex = 1.5
    ) +
    
    # North Arrow
    annotation_north_arrow(
      location = "br",
      which_north = "true",
      height = unit(2, "cm"),
      width = unit(1.5, "cm")
    ) +
    
    # Formatting
    theme_minimal() +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), 
             ylim = c(bbox["ymin"], bbox["ymax"]), 
             expand = FALSE) +
    labs(title = NULL, x = NULL, y = NULL) +
    
    
    # Adjust labels, legend, and axis breaks
    theme(
      axis.text = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.position = legend_position,
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16)
    )
  
  return(p)
}


plot_map_c <- function(raster_input,
                       boundary,
                       fill_var,
                       palette,
                       legend_title,
                       legend_position = c(0.1, 0.15),
                       fire,
                       fire2,
                       direction) {
  # Convert raster to dataframe for ggplot
  raster_df <- as.data.frame(raster_input, xy = T)
  colnames(raster_df) <- c("x", "y", "value")  # Rename value column
  
  # Get bounding box of the boundary (this will define the map limits)
  bbox <- st_bbox(boundary)  # Get the bounding box for the boundary
  
  
  # Create ggplot map
  # Create ggplot map
  p <- ggplot() +
    geom_raster(data = raster_df, aes(x = x, y = y, fill = value)) +
    scale_fill_distiller(
      palette = palette,
      direction = direction,
      limits = c(-0.25, 0.25),
      oob = scales::squish,
      name = legend_title,
    ) +
    geom_sf(
      data = st_as_sf(boundary),
      fill = NA,
      color = "black",
      size = 1
    ) +
    
    # Add the outline of the randall_shape (assuming it is an sf object)
    geom_sf(
      data = st_as_sf(fire),
      aes(color = as.factor(YEAR_)),
      fill = NA,
      lwd = 1.8
    ) +
    geom_sf(data = st_as_sf(fire2), color = 'red', fill = NA, lwd = 1.8)+
    scale_colour_viridis_d(
      option = "C",                # Ensure it's treated as categorical
      name = "Fire Year",              # Add a title for the legend
      end = 0.8
    ) +
    # Scale Bar
    annotation_scale(
      location = "bl",
      width_hint = 0.2,
      height = unit(0.5, "cm"),
      text_cex = 1.5
    ) +
    
    # North Arrow
    annotation_north_arrow(
      location = "br",
      which_north = "true",
      height = unit(2, "cm"),
      width = unit(1.5, "cm")
    ) +
    
    # Formatting
    theme_minimal() +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), 
             ylim = c(bbox["ymin"], bbox["ymax"]), 
             expand = FALSE) +  # This ensures the map stays within boundary limits
    labs(title = NULL, x = NULL, y = NULL) +
    
    
    # Adjust labels, legend, and axis breaks
    theme(
      axis.text = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.position = legend_position,
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 20)
    )
  
  return(p)
}


# Define a bivariate color palette
library(pals)
bivcol = function(pal) {
  tit = substitute(pal)
  pal = pal()
  ncol = length(pal)
  
  # Create the matrix to map colors
  z = matrix(seq_along(pal), nrow = sqrt(ncol), byrow = TRUE)  # Matrix with values from 1 to 9
  
  # Use 'image' to plot the bivariate color palette
  image(
    1:sqrt(ncol),
    1:sqrt(ncol),
    z,
    axes = FALSE,
    col = pal,
    asp = 1
  )
  
  # Add color labels on top of each color
  text(
    expand.grid(1:sqrt(ncol), 1:sqrt(ncol)),
    labels = names(pal),
    col = "black",
    cex = 1.2
  )
  mtext(tit)
}


# Function to load and clip Sentinel images
clip_sentinel_image <- function(files, study_area, date) {
  # Filter files by date
  files_by_date <- files[str_detect(files, date)]
  
  # Load images
  sentinel_stack <- rast(files_by_date)
  
  # Clip to study area
  clipped_stack <- mask(sentinel_stack, study_area)
  return(clipped_stack)
}


library(ranger)
library(pdp)
library(vip)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggcorrplot)

my_theme <- theme_minimal(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

analyze_oak_change <- function(response_raster,
                               predictor_stack,
                               response_name = "oak_change",
                               top_n = 4,
                               palette = "viridis",
                               n_samples = 2000) {
  # Mask + prepare
  mask_r <- !is.na(response_raster)
  response_raster <- mask(response_raster, mask_r)
  predictor_stack <- resample(predictor_stack, mask_r)
  full_stack <- c(response_raster, predictor_stack)
  names(full_stack)[1] <- response_name
  
  # Sample training data
  set.seed(42)
  training_data <- as.data.frame(spatSample(full_stack, method = 'random', size = n_samples, na.rm = TRUE), na.rm = TRUE)
  
  # Full correlation matrix
  cor_matrix <- cor(training_data, use = "pairwise.complete.obs")
  
  # Heatmap of correlation matrix
  cor_plot <- ggcorrplot(
    cor_matrix,
    hc.order = TRUE,
    type = "lower",
    lab = TRUE,
    lab_size = 3,
    show.legend = TRUE,
    colors = viridis(3, option = palette),
    title = paste("Correlation Matrix:", response_name),
    ggtheme = theme_minimal()
  )
  
  
  # Train random forest
  rf_model <- ranger(
    formula = as.formula(paste(response_name, "~ .")),
    data = training_data,
    importance = "permutation",
    num.trees = 500
  )
  
  # Variable importance plot (VIP)
  vip_plot <- vip(
    rf_model,
    num_features = top_n,
    geom = "col",
    aesthetics = list(fill = "#440154FF")
  ) +
    my_theme +
    scale_fill_viridis(discrete = TRUE, option = palette) +
    labs(
      title = paste("Variable Importance:", response_name),
      x = "Predictor",
      y = "Importance"
    )
  
  # Partial Dependence Plots
  top_vars <- names(sort(rf_model$variable.importance, decreasing = TRUE))[1:top_n]
  
  pdp_list <- lapply(seq_along(top_vars), function(i) {
    var <- top_vars[i]
    pd <- partial(
      rf_model,
      pred.var = var,
      train = training_data,
      grid.resolution = 50
    )
    plot <- autoplot(pd, contour = FALSE) +
      scale_color_viridis_d(option = palette) +
      labs(
        x = var,
        y = if (i == 1) "Oak cover change (TMR proxy)" else NULL
      ) +
      my_theme +
      theme(
        axis.title.y = element_text(if (i == 1) element_text() else element_blank())
      )
    return(plot)
  })
  
  pdp_grid <- plot_grid(plotlist = pdp_list, nrow = 1, align = 'hv')
  
  return(
    list(
      training_data = training_data,
      model = rf_model,
      vip_plot = vip_plot,
      pdp_grid = pdp_grid,
      correlations = cor_matrix,
      cor_plot = cor_plot
    )
  )
}

# Run multivariate classificaiton:
run_mrt_classification <- function(pc_stack,
                                   predictors,
                                   n_samples = 2000,
                                   min_split_var_exp = 0.02,
                                   seed = 42) {
  set.seed(seed)
  
  # Merge response (first 3 PCs) and predictors into one stack
  all_data <- c(pc_stack[[1:3]], predictors)
  names(all_data)[1:3] <- c("PC1", "PC2", "PC3")
  
  # Sample from stack
  sample_df <- as.data.frame(spatSample(all_data, method = "random", size = n_samples, na.rm = TRUE))
  
  # Run MRT with first 3 PCs as multivariate response
  mrt_model <- mvpart(cbind(PC1, PC2, PC3) ~ ., 
                      data = sample_df,
                      xv = "pick",   # Cross-validation to choose best size
                      xval = n_samples, 
                      cp = 0)        # Grow full tree
  
  # Prune based on custom criterion (splits explaining at least 2% variance)
  rel_error <- mrt_model$cptable[, "rel error"]
  split_var_exp <- 1 - rel_error
  keep_nodes <- which(split_var_exp >= min_split_var_exp)
  
  if (length(keep_nodes) > 0) {
    best_cp <- mrt_model$cptable[max(keep_nodes), "CP"]
    pruned_model <- prune(mrt_model, cp = best_cp)
  } else {
    warning("No splits met minimum explained variance; returning full tree.")
    pruned_model <- mrt_model
  }
  
  # Plot tree
  rpart.plot(pruned_model, main = "MRT: Ecological Classification from PCs")
  
  return(list(
    model = pruned_model,
    data = sample_df,
    pc_stack = pc_stack
  ))
}

