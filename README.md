GRACE tools 

Description

This set of tools covers a workflow for obtaining GRACE Terrestrial Water Storage Anomalies time series
The idea is to minimize user input (i.e. only a basin/aquifer/region of interest shape file or text file containing the points for the polygon). 

The tools helps the user download, filter and convert GRACE geoid coefficients into equivalent water height (mm).

Acknowledgements / Credit

The spherical harmonic analysis portion of the tools (i.e. transforming the shape file to spherical harmonic coefficients) was largely based on the original FORTRAN program dragon.F, available at 
www.iamg.org/documents/oldftp/VOL32/v32-10-11.zip

The original dragon.F was written by Wang, H. et al. (2006) for the journal paper:
Wang, H., P. Wu, and Z. Wang (2006), An approach for spherical harmonic analysis of non-smooth data, Computers & Geosciences, 32(10), 1654 ? 1668, doi: http://dx.doi.org/10.1016/j.cageo.2006.03.004.

Requirements
- R (https://www.r-project.org/), along with the following packages:
ncdf4, raster, rgdal, maptools, geosphere
- netCDF library must be installed using Intel compilers

Installation
Most of the processing is done in R scripts, so no installation is needed
The first line of the R scripts (#!/usr/.../Rscript) must be modified accordingly (path to Rscript)

The FORTRAN program must be built using Intel compilers
The Makefiles located in ./snake_src/ must be modified to point to the netcdf library, for example: 
CINC = -I/usr/local/include
LDFLAGS = -L/usr/local/lib
The script INSTALL.sh calls the Makefiles and moves the executables to ./bin, this is optional 

Using the program

1. The first thing to do is to download and process the GRACE coefficients from NASA (currently CSR Level 2 RL 05 product):

./getCSR.sh							#Creates folder ./CSR and downloads the coefficients (.gz files)

./extractCSR.sh 						#Extracts the .gz files inside the CSR folder

./texttoncdf.R CSRcoef_RAW.nc ./CSR/ 		#Assembles a netCDF file containing the coefficients. 
This script also creates a file GRACE_timevector.txt, which is useful later on (contains the time bounds for each time step). 
The first argument (name of the output netCDF file) can be named differently, but must be consistent with the following steps.
The second argument is the path to the downloaded coefficients (e.g. ./CSR/ <- notice the last slash must be present, otherwise the script will not work)

./coefdelta.R CSRcoef_RAW.nc CSRcoef_delta.nc 		#Calculate the geoid changes (e.g. delta C, delta S), removes the mean for the whole period

./runDS.R CSRcoef_delta.nc CSRcoef_DS.nc 			#Correlated-error filter (or de-striping filter). See: Swenson and Wahr (2006). This process takes a couple of minutes.

./C11.R C11.txt ./CSR/							#Downloads C10 and C11 coefficients to file C11.txt. See: Swenson, Chambers, and Wahr (2008)

./C20.R C20.txt ./CSR/ 							#Downloads C20 (SLR) coefficients to file C20.txt. See: Cheng et al. (2013)

./multiplycoef.R CSRcoef_DS.nc CSR.DSG300km.SH.nc C20.txt C11.txt 300	#  Replace correct C10, C11, and S11 coefficients in netCDF file, convert geoid changes to mm equivalent water height, apply a Gaussian filter (e.g. 300km radius)

2. The second part of the process is to convert the shape file (e.g. basin, aquifer, or region of interest) into spherical harmonic coefficients:
Assuming we are in the folder ./workspace

../basin2netcdf.R Saq_aquifer_shape.dat mask_saq.nc	#Creates a rasterized mask of the polygon. The first argument is the shape file (.shp or text file). The second argument is the name of the output netCDF file

The following line performs the spherical harmonic analysis using the FORTRAN program "snake" 
../bin/snakeSHAS mask mask_saq.nc ../AvKernels/saq.nc ../AvKernels/saq_check.nc 60 
OR
../bin/snakeSHA mask mask_saq.nc ../AvKernels/saq.nc  60 

The first argument (e.g. mask) is the name of the variable within the input netCDF file (second argument, e.g. mask_saq.nc) 
The third argument is the name of the output netCDF file containing the SH coefficients 
The fourth argument in the first option (e.g. ../AvKernels/saq_check.nc) is the name of the output file containing the SH synthesis. This is useful to visualize the averaging kernel, but is not strictly necessary for the rest of the process. 
The last argument is the maximum degree to perform the analysis (lmax). 

3. Calculate the region-averaged TWSA time series:
./calcTWSAts.R ./workspace/Saq_aquifer_shape.dat CSR.DSG300km.SH.nc  ./AvKernels/saq.nc out.txt

This last step requires the following arguments:
- the name of the shape file
- name of the GRACE netCDF file (the last one created in step 1)
- name of the basin SH netCDF file (last one created in step 2)
- name of output file

A scaling factor (see e.g. Long et al., 2015) can be supplied as an optional fifth argument.

EXTRA

The program ./bin/snakeVGRACE can be used to visualize the global TWS anomalies:
./bin/snakeVGRACE CSR.DSG300km.SH.nc CSR.DSG300km.nc 

It could also be used to visualize the effect of the de-striping filter, for example, 
./multiplycoef.R CSRcoef_delta.nc CSR.noDSG300km.SH.nc C20.txt C11.txt 300  #	notice the difference with step 1 above, here we use the coefficients before de-striping
./bin/snakeVGRACE CSR.noDSG300km.SH.nc CSR.noDSG300km.nc		


REFERENCES

Cheng, M., B. D. Tapley, and J. C. Ries (2013). Deceleration in the Earth's oblateness, J. Geophys. Res. Solid Earth, 118, 740-747, doi:10.1002/jgrb.50058.

Long, D., Longuevergne, L., and Scanlon, B. R.: Global analysis of approaches for deriving total water storage changes from GRACE satellites, Water Resources Research, doi:10.1002/2014WR016853, URL http://dx.doi.org/10.1002/2014WR016853, 2015.

Swenson, S. and Wahr, J.: Methods for inferring regional surface-mass anomalies from Gravity Recovery and Climate Experiment (GRACE) measurements of time-variable gravity, Journal of Geophysical Research, 107, 2002.

Swenson, S. and Wahr, J.: Post-processing removal of correlated errors in GRACE data, Geophysical Research Letters, 33, L08 402, 2006.

Swenson, S., D. Chambers, and J. Wahr (2008).  Estimating geocenter variations from a combination of GRACE and ocean model output, J. Geophys. Res., 113, B08410, doi:10.1029/2007JB005338.

Wang, H., P. Wu, and Z. Wang (2006), An approach for spherical harmonic analysis of non-smooth data, Computers & Geosciences, 32(10), 1654 ? 1668, doi: http://dx.doi.org/10.1016/j.cageo.2006.03.004.


