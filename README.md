# model_pol
polarization modeling 
To use this IDL based packgae, you need 
1) the output density and temperature files from the three dimensional radiation transfer code RADMC3D (http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/).
2) the integration part is called from an external fortran code. gfortran is required. 
3) the dust properties (absorption & emissive efficiencies) are calculated from the package DDSCAT (http://ddscat.wikidot.com/).  
STEPS
1. write the magnetic field stuctures in cartesian cooridnates (ms_b.pro)
2. regrid the density and temperature distribution from spherical coordinate to cartesian coordinate (pm_coord_regrid.pro)
3. use pm_integrate to do the simple ray-tracing and calculate the resultant stokes parameters I, U, and Q along the line of sight (pm_integrate.pro). 
4. reduce the I, U, and Q images and obtain the spatial polarization percentage p and polarization position angle theta (pm_image.pro).  
