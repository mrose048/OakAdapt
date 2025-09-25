##%######################################################%##
#                                                          #
####      Spectral/Temporal Mixture Model Function      ####
#                                                          #
##%######################################################%##

# w = weight unit contraint (set to 1)
# n_layers = number of layers in the raster stack that will be unmixed
# rescale = 1 (temporal) or 1/10000 (spectral) # rescale factor for endmembers

unmix <- function(endmember_file,
                  D1_file,
                  n_layers,
                  w = 1,
                  rescale = 1) {
  #Multiply EMs by a scale factor?
  # rescale = 1/10000
  
  bands = data.frame(matrix(1, nrow = 1, ncol = n_layers))
  
  #Read endmember file
  G = read.csv(endmember_file, header = F, sep = "", skip = 6)
  
  # Read spatial data file
  D1 <- rast(D1_file)  # Read vegetation fraction stack as a raster
  crs_D1 <- crs(D1)  # Store the original CRS
  
  # Get wavelengths from the first column of the endmember file
  wavelengths <- G[, 1]
  G2 <- G[, 2:ncol(G)]
  
  # Rescale endmembers
  # G3 <- G2 * rescale
  G3 <- G2
  
  # Subset using bad bands
  D2 <- as.array(D1)[, , bands == 1]  # Extract valid bands
  gc()
  G4 <- G3[bands == 1, ]
  
  # Reshape arrays
  D3 <- array_reshape(D2, c(dim(D2)[1] * dim(D2)[2], dim(D2)[3]), order = "C")
  gc()
  
  # Add unit sum constraint
  D4 <- t(as.matrix(cbind(D3, rep(1, nrow(
    D3
  )))))  # Transpose to add unit sum constraint
  gc()
  G5 <- as.matrix(rbind(G4, rep(1, ncol(G4))))
  
  # Compute residual of time series - you will transform this with PCA and evaluate the first three PC
  r <- D4 - (G5 %*% solve(t(G5) %*% G5) %*% t(G5) %*% D4)
  
  gc()
  
  # Transpose to matrix to match D3
  #r2 <- t(r)
  r2 <- t(r[-nrow(r), ])  # Exclude the last row (unit sum constraint)
  
  gc()
  
  # reshape back to original dimensions
  r_reshape <- array_reshape(r2, c(dim(D1)[1], dim(D1)[2], dim(r2)[2]), order = "C")
  
  # RMSE
  # Square residuals, mean across bands, take square root
  rmse_map <- sqrt(apply(r_reshape^2, c(1, 2), mean, na.rm = TRUE))
  
  # Convert to raster
  rmse_raster <- rast(rmse_map, ext = ext(D1), crs = crs_D1)
  # plot(rmse_raster)
  
  # Compute RMSE from residual matrix
  rmse_all <- sqrt(mean(r2^2, na.rm = TRUE))
  
  gc()
  
  # Convert back to raster while preserving CRS
  r_raster <- rast(r_reshape, ext = ext(D1), crs = crs_D1)
  
  # Compute mixture fractions
  u <- solve(t(G5) %*% G5) %*% t(G5) %*% D4
  
  gc()
  
  # Transpose the matrix to match D3 format
  u2 <- t(u)
  
  gc()
  
  # Reshape back to original dimensions
  u_reshape <- array_reshape(u2, c(dim(D1)[1], dim(D1)[2], dim(u2)[2]), order = "C")
  
  # Convert back to raster while preserving CRS
  u_raster <- rast(u_reshape, ext = ext(D1), crs = crs_D1)
  
  names(r_raster) <- names(D1)
  
  s <- list(r_raster, u_raster)
  names(s) <- c("residuals", "unmixing_results")
  
  return(list(s, rmse_all))
  
}
