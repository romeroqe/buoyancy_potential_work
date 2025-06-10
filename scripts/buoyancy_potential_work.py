import numpy as np
from scipy import interpolate


"""
Function to calculate the WORK DONE BY BUOYANCY (WB)
WB quantifies the water column's vertical homogeneity in terms of the work done by the buoyancy force in vertically displacing a water parcel under static instability conditions.
"""
def buoyancy_potential_work(rho, z, zint=-10):
	"""
	Parameters:
    Depth (z) and density (rho) data should go from the deepest to the shallowest and be column vectors
	- rho: Sigma-0 Potential Density Anomaly profile (kg*m^-3)
	- z: Heights (m) with negative units. The z-vector can be non-equidistant.
	- zint: Reference height (default = -10 m)

	Return:
	- WB: Work done by buoyancy profile (J*m^-3)
	- z: Heights (m)
	"""
	
	# A shallower depth than that of interest is required
	if z[0] < zint: 
		return [], []

	# If there is no data at zint, it is interpolated
	if np.nansum(z[z == zint]) == 0:
		f = interpolate.interp1d(z, rho)
		i10 = np.argwhere(z <= zint)[0]
		z = np.insert(z, i10, zint)
		rho = np.insert(rho, i10, f(zint))
	
	WB = np.zeros_like(z)*np.NaN
	nz  = rho.shape[0] # Amount of data in the vertical
	g = 9.81 # Acceleration of gravity
	
	izint = np.nonzero(z==zint)
	dz = np.diff(z) # Depth vector spacing (m) with positive units
	
	# Calculate buoyancy work for various depths of interest
	for it in range(nz):
		zint = z[it]
		if np.isnan(zint): 
			#If the depth of interest is nan, the work is also nan
			WB[it] = np.nan
		else:
			rho_int = rho[it]
			
			# Computation of WB
			S = np.nan*np.zeros_like(z)
			S[it] = 0
			
			for i in range(0,it):
				sr=rho[i:it]+rho[i+1:it+1]
				sdz=dz[i:it]
				S[i] = -0.5*np.nansum(sr*sdz)
			for i in range(it+1,nz):
				sr=rho[it:i]+rho[it+1:i+1]
				sdz=dz[it:i]
				S[i] = 0.5*np.nansum(sr*sdz)

			V = -(g*(z-zint)*rho_int)+(g*S) # Buoyancy potential (J*m^-3)
			V = -V # Change the sign of the potential to make it more intuitive
			WB[it] = V[izint] # Work done by buoyancy (J*m^-3)

	return WB, z
