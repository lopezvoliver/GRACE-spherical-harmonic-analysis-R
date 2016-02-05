Workflow for generating basin-averaged Terrestrial Water Storage (TWS) anomalies from GRACE CSR data:

1. Download geoid coefficients from ftp://podaac.jpl.nasa.gov/allData/grace/L2/CSR/RL05/
getCSR.sh
extractCSR.sh

2. Assemble a netCDF file
./texttoncdf.R CSRcoef_RAW.nc 

3. Calculate the geoid changes (delta C, delta S) and remove the mean for the whole period
./coefdelta.R CSRcoef_RAW.nc CSRcoef_delta.nc

4. Remove correlated errors using the de-striping filter:
./runDS.R 

5. Download correct C10, C11, and S11 coefficients:
./C11.R C11.txt
./C20.R C20.txt

6. Replace correct C10, C11, and S11 coefficients in netCDF file, convert geoid changes to mm equivalent water height, apply a Gaussian filter:
./multiplycoef.R CSRcoef_delta.nc TWSA300km.nc C20.txt C11.txt 200;

7. Read shape file (e.g. basin, aquifer, or region of interest) and calculate averaging kernel


8. Calculate region-averaged TWSA time series:
./sumcoefs.R source TWSA300km.nc saq_avkernel.nc saq_TWSA.txt
