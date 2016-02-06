GRACE tools 

This set of tools covers a workflow for obtaining GRACE Terrestrial Water Storage Anomalies time series. 

The user only needs a set of points (or shape file) for the region of interest.

The tools helps the user download, filter and convert GRACE geoid coefficients into equivalent water height (mm).

Requirements:
R needs to be installed (https://www.r-project.org/), along with the following packages:
- ncdf4
- (to do: complete this list)
Netcdf library must be installed and its module loaded
Intel Fortran compiler (ifort)

Installation:
The FORTRAN program 'snake' must be compiled using the Makefile
Must be built with ifort, as building with gfortran gives inaccurate results

For all R executables, change the first line accordingly to provide the correct path to Rscript (e.g. #!/usr/bin/Rscript)


The step-by-step process is as follows:

1. Download geoid coefficients (currently CSR RL 05 product) 
getCSR.sh     This script creates a folder called CSR on the working directory
extractCSR.sh This script extracts the files inside the folder CSR

2. Assemble a netCDF file
./texttoncdf.R CSRcoef_RAW.nc ./CSR/  
This script reads the downloaded files and creates a netCDF file. The first argument is the name of the output file (e.g. CSRcoef_RAW.nc). The second argument is the path to the downloaded coefficients (e.g. ./CSR/ <- notice the last slash must be present, otherwise the script will not work)


3. Calculate the geoid changes (delta C, delta S) and remove the mean for the whole period
./coefdelta.R CSRcoef_RAW.nc CSRcoef_delta.nc

4. Remove correlated errors using the de-striping filter:
./runDS.R CSRcoef_delta.nc CSRcoef_DS.nc 

5. Download correct C10, C11, and S11 coefficients:
./C11.R C11.txt ./CSR/
./C20.R C20.txt ./CSR/

6. Replace correct C10, C11, and S11 coefficients in netCDF file, convert geoid changes to mm equivalent water height, apply a Gaussian filter (e.g. 300km radius):
./multiplycoef.R CSRcoef_DS.nc CSR.DSG300km.coef.nc C20.txt C11.txt 300;


--- After this point we need the FORTRAN program to be built

8. Visualize global TWSA:
To do.. 

9. Read shape file (e.g. basin, aquifer, or region of interest) and calculate averaging kernel
9.1 Read shape file and create a masked netcdf file
./basin2netcdf.R ./WS/SAQ.dat 1 mask_saq.nc
mv mask_saq.nc ./WS/mask_saq.nc
9.2 Transform to SH 
./bin/snakeSHAS mask ./WS/mask_saq.nc ./AvKernels/saq.nc ./AvKernels/saq_check.nc 60

10. Calculate region-averaged TWSA time series:
./sumcoefs.R source TWSA300km.nc saq_avkernel.nc saq_TWSA.txt
