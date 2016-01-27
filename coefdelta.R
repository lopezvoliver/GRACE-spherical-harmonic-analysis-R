#!/usr/bin/env Rscript 	
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<2){
   cat('Usage: ./coefdelta.R input.nc output.nc\n')
   quit()
}
fnamei = args[1]
fnameo = args[2]
library(ncdf4)
#################################################################################
#										#
# GRACE GSM files represent the full gravity field as Stokes coefficients	# 
#							   Clm(t), Slm(t)	#
# Here we apply the following post-processing needed to obtain the signal in 	#
# equivalent water height:							#
#										#
# 1) Remove time-mean in order to get Stokes coefficients anomalies: 		#
#							  ΔClm(t), ΔSlm(t)	#
nci=nc_open(fnamei)
day1<-ncvar_get(nci,'day1')
day2<-ncvar_get(nci,'day2')
time<-ncvar_get(nci,'time')
Clm<-ncvar_get(nci,'clm')
Slm<-ncvar_get(nci,'slm')
nc_close(nci)

# Step 1: Remove temporal mean 
mClm=array(rep(colMeans(aperm(Clm,c(3,1,2))),times=dim(Clm)[3]),dim=dim(Clm))
mSlm=array(rep(colMeans(aperm(Slm,c(3,1,2))),times=dim(Slm)[3]),dim=dim(Slm))
dClm=Clm-mClm
dSlm=Slm-mSlm

#we will store the coefficients small netcdf file with dimensions l,m,t
#Define dimensions	
dimt<-ncdim_def('time','GRACE month',1:length(time))
diml<-ncdim_def('degree','l',0:60)
dimm<-ncdim_def('order','m',0:60)
	
#Define variables 
ncday1=ncvar_def('day1','days since 2000-01-01',list(dimt),-9999,longname='GRACE day start',prec='integer')
ncday2=ncvar_def('day2','days since 2000-01-01',list(dimt),-9999,longname='GRACE day end',prec='integer')	
ncoef1=ncvar_def('clm','geoid height anomaly',list(diml,dimm,dimt),-9999,longname='delta Clm Temporal mean from 2002120-2015243 removed',prec='float')
ncoef2=ncvar_def('slm','geoid height anomaly',list(diml,dimm,dimt),-9999,longname='delta Slm Temporal mean from 2002120-2015243 removed',prec='float')
	
#Define nc file
ncout<-nc_create(fnameo,list(ncday1,ncday2,ncoef1,ncoef2),force_v4=TRUE)

#save day periods and de-trended coefficients
ncvar_put(ncout,'day1',day1)
ncvar_put(ncout,'day2',day2)
ncvar_put(ncout,'clm',dClm)
ncvar_put(ncout,'slm',dSlm)

#Close file
nc_close(ncout)

command<-sprintf('ncdump -h %s',fnameo)
system(command)

