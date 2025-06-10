# Calculates the buoyancy work required to bring a parcel of water from different depths (zint) 
# to a reference depth (z2=-10 m). It is considered a resting ocean and that the plot in zint is in hydrostatic 
# equilibrium. The vertical discretization of the vertical profile can be non-equidistant.

buoyancy_potential_work <- function(rho, z, zint=-10) {
	# Parameters: 
  # Depth (z) and density (rho) data should go from the deepest to the shallowest and be column vectors
	# - rho: Sigma-0 Potential Density Anomaly profile (kg*m^-3)
	# - z: Heights (m) with negative units. The z-vector can be non-equidistant.
  # 
	# Return:
	# - WB: Work done by buoyancy profile (J*m^-3)
	# - z: Heights (m)
  # - zint: Reference height (default = -10 m)
  
	# A shallower depth than that of interest is required
  if (z[1] < zint) {
    return(list(WB=numeric(0), z=numeric(0)))
  }
  
	# If there is no data at zint, it is interpolated
  if (!(zint %in% z)) {
    rho_int <- approx(z, rho, xout=zint)$y
    insert_pos <- which(z <= zint)[1]
    z <- append(z, zint, after=insert_pos - 1)
    rho <- append(rho, rho_int, after=insert_pos - 1)
  }

  #### 1. Define variables ####
  nz <- length(rho) # Amount of data in the vertical
  WB <- rep(NA, nz) # Buoyancy work (J*m^-3)
  g <- 9.81 # Acceleration of gravity (m*s^-2)
  
  izint <- which(z==zint) # Index of zint in the depth vector z
  dz <- diff(z) # Spacing of the depth vector (m) with positive units

  for (it in 1:nz) { # counter for the cycle, which runs through all the depths of z
    
    # 2.1 Given the depth of interest zint, find its corresponding density
    zint <- z[it] # Depth of interest
    if (is.na(zint)) {
      # If the depth of interest is NA, the work is also NA
      WB[it] <- NA
    } else {
      rho_int <- rho[it] # Density in the depth of interest
      
      # 2.2 Calculation of the buoyancy potential for zint
      S <- rep(NA, nz)
      S[it] <- 0
      
      for (i in 1:(it-1)) {
        sr <- rho[i:(it-1)] + rho[(i+1):it]
        sdz <- dz[i:(it-1)]
        S[i] <- -0.5 * sum(sr * sdz, na.rm = TRUE)
      }
      
      for (i in (it+1):nz) {
        sr <- rho[it:(i-1)] + rho[(it+1):i]
        sdz <- dz[it:(i-1)]
        S[i] <- 0.5 * sum(sr * sdz, na.rm = TRUE)
      }
      
      V <- -(g * (z - zint) * rho_int) + (g * S[1:nz]) # Buoyancy potential (J*m^-3)
      V <- -V # Changing the sign of potential to make it more intuitive
      
      # 2.3 Calculating the work from z1 to z2
      WB[it] <- V[izint] # Buoyancy work (J*m^-3)
    }
  }

  return(list(WB=WB, z=z))
}
