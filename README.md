# model_pol
polarization modeling 
To use this IDL based packgae, you need 1) the output density and temperature files from the three dimensional radiation transfer code RADMC3D http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/
2) the integration part is called from an external fortran code. fortran 77 is required. 
3) the dust properties (absorption & emissive efficiencies) are calculated from the package DDSCAT http://ddscat.wikidot.com/.  
STEPS
1. write the magnetic field stuctures in cartesian cooridnates (ms_b.pro)
2. regrid the 
