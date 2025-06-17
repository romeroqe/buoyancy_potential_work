# Function to calculate the WORK DONE BY THE BUOYANCY FORCE (WB) in
# vertically displacing a water parcel from different depths to a reference
# depth zref.
buoyancy_potential_work <- function(rho, z, zref=-10) {
	# WB quantifies the water column's vertical homogeneity in terms of the work
	# done by the buoyancy force in vertically displacing a water parcel under 
	# static instability conditions.
	#
	# Inputs:
	# * rho: Sigma-0. Potential Density Anomaly profile referred to 0 dbar (kg m^-3)
	# * z: Depths in negative values (m)
	# * zref: Reference depth (default -10 m). WB(zref)=0 by definition.
	# The data must be a column vector and be arranged from the greatest to the shallowest depth.
	#
	# Outputs:
	# * WB: Vertical profile of work done by buoyancy (J m^-3)
	# * z: Depths (m)

	# If there is no data at zref, it is interpolated
	if (!(zref %in% z)) {
		# Interpolate rho at zref
		rho_zref <- approx(z, rho, xout=zref, rule=2)$y

		# Construct the vectors z and rho, adding zref and rho_zref
		z_extended <- c(zref, z)
		rho_extended <- c(rho_zref, rho)

		# Sort vectors z and rho
		order_idx <- order(z_extended)
		z <- z_extended[order_idx]
		rho <- rho_extended[order_idx]
	}

	# Defining variables
	nz <- length(rho) # Number of data in the vertical
	WB <- rep(NA, nz) # Work done by buoyancy (J m^-3)
	g <- 9.81 # Acceleration of gravity (m s^-2)
	izref <- which(z==zref) # Index of zref in the vector z
	dz <- diff(z) # Differences between adjacent elements of z (m), with positive units

	# Calculation of the vertical profile of WB
	for (it in 1:nz) { # Loop for every depth
		zref <- z[it] # Depth of interest
		if (is.na(zref)) { # If the depth of analysis is nan, WB is also nan
			WB[it] <- NA
		} else {
			rho_int <- rho[it]

			# Computation of WB
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

			V <- -(g * (z - zref) * rho_int) + (g * S[1:nz]) # Buoyancy potential (J m^-3)
			V <- -V # Change the sign of the buoyancy potential
			WB[it] <- V[izref] # WB at the depth of analysis (J m^-3)
		}
	}
	
	return(list(WB=WB, z=z))
}
