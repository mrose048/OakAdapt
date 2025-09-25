##%######################################################%##
#                                                          #
####        Functions for handling Landsat data         ####
#                                                          #
##%######################################################%##

library(terra)
library(stringr)
library(dplyr)
library(progress)
library(abind) 
library(lubridate)


remove_duplicate_files <- function(dir_path) {
  dup_patterns <- paste0(" \\(", 1:5, "\\)\\.")
  all_files <- list.files(dir_path, full.names = TRUE, recursive = TRUE)
  for (pattern in dup_patterns) {
    file.remove(all_files[grepl(pattern, all_files)])
  }
}

get_unique_ids <- function(file_list) {
  strsplit(basename(file_list), "_") %>%
    sapply(function(x) paste(x[1], x[4], sep = "_"))
}

compare_sources <- function(dir1, dir2) {
  f1 <- list.files(dir1, full.names = FALSE) %>% get_unique_ids() %>% unique()
  f2 <- list.files(dir2, full.names = FALSE) %>% get_unique_ids() %>% unique()
  list(
    only_in_dir1 = setdiff(f1, f2),
    only_in_dir2 = setdiff(f2, f1),
    shared = intersect(f1, f2)
  )
}

organize_landsat_files <- function(files, out_dir) {
  require(progress)
  require(stringr)
  
  # Skip anything that doesn't match Landsat pattern or contains T2
  files <- files[
    grepl("^L(C08|C09|E07|T05|T04).*\\.TIF$", basename(files), ignore.case = TRUE) &
      !grepl("T2", basename(files), ignore.case = TRUE)
  ]
  
  pb <- progress_bar$new(
    format = "Organizing files [:bar] :current/:total (:percent) eta: :eta",
    total = length(files),
    clear = FALSE,
    width = 60
  )
  
  df <- data.frame(file = files, base = basename(files))
  df$folder <- str_extract(df$base, ".*T1")  # Adjust if your structure varies
  folders <- unique(df$folder)
  
  for (i in seq_len(nrow(df))) {
    folder_path <- file.path(out_dir, df$folder[i])
    dir.create(folder_path, showWarnings = FALSE, recursive = TRUE)
    
    target_file <- file.path(folder_path, df$base[i])
    
    if (!file.exists(target_file)) {
      file.copy(df$file[i], target_file, overwrite = FALSE)
    }
    
    pb$tick()
  }
  
  return(folders)
}

mosaic_by_folder <- function(folders, base_dir, shape, out_dir) {
  require(progress)
  pb <- progress_bar$new(
    format = "Mosaicking [:bar] :current/:total (:percent) eta: :eta",
    total = length(folders),
    clear = FALSE,
    width = 60
  )
  
  for (f in folders) {
    pb$tick()
    files <- list.files(file.path(base_dir, f), pattern = "_SR_B.*\\.TIF$", full.names = TRUE)
    if (length(files) == 0) next
    
    out_file <- file.path(out_dir, paste0(f, ".tif"))
    if (file.exists(out_file)) next
    
    r <- tryCatch(rast(files), error = function(e) NULL)
    if (is.null(r)) next
    
    # Layer check
    expected_layers <- if (grepl("LC08|LC09", f)) 7 else 6
    if (nlyr(r) != expected_layers) next
    
    cropped <- crop(r, shape)
    writeRaster(cropped, out_file, overwrite = TRUE)
  }
}


scale_reflectance <- function(in_dir, out_dir) {
  require(progress)
  files <- list.files(in_dir, pattern = "T1.tif$", full.names = TRUE)
  
  pb <- progress_bar$new(
    format = "Scaling reflectance [:bar] :current/:total (:percent) eta: :eta",
    total = length(files),
    clear = FALSE,
    width = 60
  )
  
  for (i in seq_along(files)) {
    pb$tick()
    out_file <- gsub("\\.tif$", "_Scale", files[i])
    if (file.exists(out_file)) next
    
    r <- tryCatch(rast(files[i]), error = function(e) NULL)
    if (is.null(r)) next
    
    rout <- (r * 0.0000275) - 0.2
    writeRaster(rout, out_file, filetype = "ENVI", overwrite = TRUE)
    gc()
  }
}


library(terra)

spectral_unmixing <- function(files, out_dir, shape, rescale = 1/10000) {
  # Define badbands and endmember files per sensor mission
  sensor_params <- list(
    LT04 = list(
      badbands = '1_Data/00_Badbands/L5bbl.csv',
      endmember = '1_Data/00_Endmembers/EM_MR_L4_txtNH.txt'
    ),
    LT05 = list(
      badbands = '1_Data/00_Badbands/L5bbl.csv',
      endmember = '1_Data/00_Endmembers/EM_MR_L5_txtNH.txt'
    ),
    LE07 = list(
      badbands = '1_Data/00_Badbands/L7bbl.csv',
      endmember = '1_Data/00_Endmembers/EM_MR_L7_txtNH.txt'
    ),
    LC = list(  # Shared for LC08 and LC09
      badbands = '1_Data/00_Badbands/L8bbl.csv',
      endmember = '1_Data/00_Endmembers/EM_MR_L8_txtNH.txt'
    )
  )
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  pb <- txtProgressBar(min = 0, max = length(files), style = 3)
  
  for (i in seq_along(files)) {
    f <- files[i]
    cat(sprintf("\nProcessing file [%d/%d]: %s\n", i, length(files), f))
    
    if (!file.exists(f)) {
      cat(sprintf("❌ File not found: %s\n", f))
      next
    }
    
    filename <- basename(f)
    filename_no_ext <- tools::file_path_sans_ext(filename)
    extension <- tools::file_ext(filename)
    
    parts <- unlist(strsplit(filename_no_ext, "_"))
    
    if (length(parts) >= 8) {
      sensor_mission <- parts[1]
      sensor_key <- if (grepl("^LC", sensor_mission)) "LC" else sensor_mission
      
      if (!(sensor_key %in% names(sensor_params))) {
        cat(sprintf("⚠️ Unrecognized sensor mission code: %s in file %s\n", sensor_mission, filename))
        next
      }
      
      # Try loading badbands and endmembers
      bands <- tryCatch({
        read.csv(sensor_params[[sensor_key]]$badbands, header = FALSE, stringsAsFactors = FALSE)
      }, error = function(e) {
        cat(sprintf("❌ Failed to read badbands for sensor %s: %s\n", sensor_key, e$message))
        next
      })
      
      G <- tryCatch({
        read.csv(sensor_params[[sensor_key]]$endmember, header = FALSE, sep = "", stringsAsFactors = FALSE)
      }, error = function(e) {
        cat(sprintf("❌ Failed to read endmember file for sensor %s: %s\n", sensor_key, e$message))
        next
      })
      
      product <- parts[2]
      path_row <- parts[3]
      date1 <- parts[4]
      date2 <- parts[5]
      version <- parts[6]
      tier <- parts[7]
      scale <- paste(parts[8:length(parts)], collapse = "_")
      unmixed <- "Unmixed"
      
      new_filename <- paste(product, date1, date2, version, tier, sensor_mission, path_row, scale, unmixed, sep = "_")
    } else {
      cat(sprintf("⚠️ Filename format issue, skipping: %s\n", filename))
      next
    }
    
    outfile_path <- file.path(out_dir, new_filename)
    
    if (file.exists(outfile_path)) {
      setTxtProgressBar(pb, i)
      next
    }
    
    # Load raster with error handling
    D1 <- tryCatch({
      rast(f)
    }, error = function(e) {
      cat(sprintf("❌ Failed to load raster %s: %s\n", f, e$message))
      next
    })
    
    crs_D1 <- crs(D1)
    
    selected_bands <- which(bands == 1)
    G2 <- G[, -1]
    G3 <- G2[selected_bands, , drop = FALSE] * rescale
    
    D2 <- tryCatch({
      as.array(D1)[, , selected_bands, drop = FALSE]
    }, error = function(e) {
      cat(sprintf("❌ Failed to extract bands from raster: %s\n", e$message))
      next
    })
    
    dim_D2 <- dim(D2)
    D3 <- matrix(D2, nrow = dim_D2[1]*dim_D2[2], ncol = dim_D2[3])
    D4 <- t(cbind(D3, rep(1, nrow(D3))))
    G5 <- rbind(as.matrix(G3), rep(1, ncol(G3)))
    
    u <- tryCatch({
      solve(t(G5) %*% G5) %*% t(G5) %*% D4
    }, error = function(e) {
      cat(sprintf("❌ Failed to solve unmixing equation for %s: %s\n", filename, e$message))
      next
    })
    
    u2 <- t(u)
    u_reshape <- array(u2, dim = c(dim_D2[1], dim_D2[2], ncol(G3)))
    u_raster <- rast(u_reshape, ext = ext(D1), crs = crs_D1)
    
    u_crop <- tryCatch({
      crop(u_raster, shape)
    }, error = function(e) {
      cat(sprintf("❌ Failed to crop raster: %s\n", e$message))
      next
    })
    
    u_mask <- tryCatch({
      mask(u_crop, shape)
    }, error = function(e) {
      cat(sprintf("❌ Failed to mask raster: %s\n", e$message))
      next
    })
    
    tryCatch({
      writeRaster(u_mask, filename = outfile_path, overwrite = TRUE, datatype = "FLT4S", filetype = 'ENVI')
    }, error = function(e) {
      cat(sprintf("❌ Failed to write raster: %s\n", e$message))
      next
    })
    
    gc()
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  cat("\n✅ Spectral unmixing completed for", length(files), "files.\n")
}




stack_veg_fractions <- function(site_name,
                                unmixed_dir,
                                clear_img_dir,
                                out_dir,
                                out_file_base) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # STEP 1: Build list of clear image identifiers
  clear_list <- list.files(clear_img_dir, pattern = "\\.jpg$", full.names = FALSE)
  clear_list <- gsub("\\.jpg$", "", clear_list)
  clear_list <- gsub("L1TP|L1GS", "L2SP", clear_list)
  clear_list <- clear_list[!grepl("T2| \\(1\\)$", clear_list)]
  clear_ids <- apply(str_split_fixed(clear_list, "_", 5)[, 1:4], 1, paste, collapse = "_")
  clear_dates <- sub(".*_(\\d{8})$", "\\1", clear_ids)
  
  # STEP 2: Get matching unmixed files
  unmixed_files <- list.files(unmixed_dir, pattern = "Unmixed$", full.names = TRUE)
  selected_files <- unmixed_files[sapply(unmixed_files, function(f) any(str_detect(f, paste(clear_dates, collapse = "|"))))]
  
  if (length(selected_files) == 0) {
    stop(paste("No unmixed files matched for", site_name))
  }
  
  # STEP 3: Calculate fractional year
  calc_frac_year <- function(filename) {
    match <- str_extract(filename, "\\d{8}")
    if (is.na(match)) return(NA)
    date <- ymd(match)
    year(date) + (yday(date) - 1) / ifelse(leap_year(date), 366, 365)
  }
  frac_years <- sapply(selected_files, calc_frac_year)
  
  # STEP 4: Stack vegetation fraction (assumed to be 2nd band)
  veg_rasters <- lapply(selected_files, function(f) rast(f)[[2]])
  veg_stack <- rast(veg_rasters)
  names(veg_stack) <- frac_years
  
  # STEP 5: Create and apply mask
  mask_raster <- any(is.na(veg_stack))
  # writeRaster(mask_raster, file.path(out_dir, paste0(out_file_base, "_mask")), 
  #             overwrite = TRUE, filetype = "ENVI")
  
  masked_stack <- mask(veg_stack, mask_raster, maskvalues = 1)
  
  # STEP 6: Write output ENVI file
  final_envi_file <- file.path(out_dir, paste0(out_file_base))
  writeRaster(masked_stack, final_envi_file, datatype = "FLT4S", overwrite = TRUE, filetype = "ENVI")
  
  # STEP 7: Modify header to include wavelengths
  hdr_file <- paste0(final_envi_file, ".hdr")
  if (file.exists(hdr_file)) {
    hdr <- readLines(hdr_file)
    idx <- grep("data ignore value = nan", hdr)
    
    # Get band names as wavelengths
    wavelengths <- names(masked_stack)
   # wavelength_str <- paste(wavelengths, collapse = ", ")
    
    # Build block and insert
    wavelength_block <- c("wavelength = {", paste0("  ", wavelengths), "}")
    hdr <- append(hdr, wavelength_block, after = idx)
    
    # Write updated header
    writeLines(hdr, hdr_file)
  }
  
  message(paste("Finished stacking for", site_name, "with", length(selected_files), "images."))
}


perform_pca_on_stack <- function(stack, out_prefix, out_dir) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Get the values of the raster stack as a matrix (rows = pixels, columns = time)
  stack_values <- as.matrix(stack)
  
  # Remove any rows (pixels) with NA across time
  valid_idx <- complete.cases(stack_values)
  valid_data <- stack_values[valid_idx, ]
  
  # Perform PCA
  pca <- prcomp(valid_data, center = TRUE, scale. = TRUE)
  
  # Get variance explained
  variance <- pca$sdev^2
  prop_var <- variance / sum(variance)
  var_df <- data.frame(
    PC = paste0("PC", seq_along(prop_var)),
    Eigenvalue = variance,
    ProportionExplained = prop_var,
    CumulativeExplained = cumsum(prop_var)
  )
  
  # Save variance explained
  write.csv(var_df, file = file.path(out_dir, paste0(out_prefix, "_PCA_variance_explained.csv")), row.names = FALSE)
  
  # Take first 3 components and insert into empty raster
  pcs <- matrix(NA, nrow = ncell(stack), ncol = 3)
  pcs[valid_idx, ] <- pca$x[, 1:3]
  
  # Create raster stack of PC1–3
  pc_rasters <- rast(pcs, nrows = nrow(stack), ncols = ncol(stack),
                     ext = ext(stack), crs = crs(stack))
  names(pc_rasters) <- paste0("PC", 1:3)
  
  # Write raster stack
  writeRaster(pc_rasters,
              filename = file.path(out_dir, paste0(out_prefix, "_PCA_stack")),
              overwrite = TRUE,
              filetype = "ENVI",
              datatype = "FLT4S")
  
  return(invisible(var_df))
}



