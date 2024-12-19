# A script that calculates the buoyancy work required to bring a parcel of water from different depths (zint) 
# to a reference depth (z2=-10 m). It is considered a resting ocean and that the plot in zint is in hydrostatic 
# equilibrium. The vertical discretization of the vertical profile can be non-equidistant.

library(R.matlab)
library(ggplot2)

#### 1. Upload data ####
# Argo Profiles
# These Argo profiles were interpolated to a value of -10 m.

iargo <- 3

if (iargo == 1) {
  filename <- 'D5903264_520.mat' # every 2 m
} else if (iargo == 2) {
  filename <- 'D1900039_112.mat' # every 10 m
} else {
  filename <- 'D1900044_079.mat' # varies across the profile
}

mat <- readMat(filename)

# Depth and density data should go from the deepest to the shallowest and be column vectors
# z: Depth file (m) with negative units
# rho: Sigma-0 Potential Density Anomaly File (kg*m^-3)

z <- as.vector(mat$z)
rho <- as.vector(mat$rho)

#### 2. Define variables ####
g <- 9.81 # Acceleration of gravity (m*s^-2)
nz <- length(rho) # Amount of data in the vertical
WV <- rep(NA, nz) # Buoyancy work (J*m^-3)
z2 <- -10 # Depth (m) for job calculation from different depths of interest to z2
iz2 <- which(z == z2) # Index of z2 in the depth vector z
dz <- diff(z) # Spacing of the depth vector (m) with positive units

#### 3. Calculate work per buoyancy for various depths of interest zint ####
for (it in 1:nz) { # counter for the cycle, which runs through all the depths of z
  
  # 3.1 Given the depth of interest zint, find its corresponding density
  zint <- z[it] # Depth of interest
  rho_int <- rho[it] # Density in the depth of interest
  
  # 3.2 Buoyancy strength (N*m^-3)
  F <- g * (rho - rho_int)
  
  # 3.3 Calculation of the buoyancy potential for zint
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
  
  # 3.4 Calculating the work from z1 to z2
  WV[it] <- V[iz2] # Buoyancy work (J*m^-3)
}

## --- 4. Graphics ---
ggplot() +
  geom_line(aes(x = rho, y = z), color = 'blue', size = 1) +
  labs(x = expression(Density ~ (kg ~ m^{-3})), y = "Depth (m)") +
  theme_minimal() +
  theme(axis.title.x = element_text(color = 'blue'),
        axis.text.x = element_text(color = 'blue'),
        axis.line.x.bottom = element_line(color = 'blue'))

ggplot() +
  geom_line(aes(x = WV, y = z), color = 'red', size = 1) +
  labs(x = expression(Work ~ (N ~ m^{-3})), y = "Depth (m)") +
  theme_minimal() +
  theme(axis.title.x = element_text(color = 'red'),
        axis.text.x = element_text(color = 'red'),
        axis.line.x.top = element_line(color = 'red'))
