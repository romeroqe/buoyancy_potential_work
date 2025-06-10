library(R.matlab)

source("../scripts/buoyancy_potential_work.R")

#### 1. Upload data ####
# These Argo profiles were interpolated to a value of -10 m.

iargo <- 3

if (iargo == 1) {
  filename <- '../demo/data/D5903264_520.mat' # every 2 m
} else if (iargo == 2) {
  filename <- '../demo/data/D1900039_112.mat' # every 10 m
} else {
  filename <- '../demo/data/D1900044_079.mat' # varies across the profile
}

#### 2. Define variables ####
mat <- readMat(filename)
z <- rev(as.vector(mat$z))
rho <- rev(as.vector(mat$rho))

#### 3. Compute the Work done by buoyancy ####
data <- buoyancy_potential_work(rho, z)

#### 4. Figures ####
plot(rho, z, type = "l", col = "blue", lwd = 2,
     xlab = expression(Density ~ (kg ~ m^{-3})), ylab = "Depth (m)")
axis(1, col.axis = "blue", col = "blue")  # eje x en azul
title(xlab = expression(Density ~ (kg ~ m^{-3})), col.lab = "blue")

plot(data$WB, data$z, type = "l", col = "red", lwd = 2,
     xlab = expression(Work ~ (J ~ m^{-3})), ylab = "Depth (m)")
axis(1, col.axis = "red", col = "red")  # eje x rojo
title(xlab = expression(Work ~ (J ~ m^{-3})), col.lab = "red")

