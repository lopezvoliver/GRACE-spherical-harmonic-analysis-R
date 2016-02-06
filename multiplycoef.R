#!/usr/bin/env Rscript 	
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<5){
   cat('./multiplycoef.R fnamei fnameo c20f c11f gausr')
   cat('Provide file name of input and output\n')
   cat('Also provide: C20 file (txt) and C11 file (txt)\n')
   quit()
}
fnamei = args[1]
fnameo = args[2]
c20f = args[3]
c11f = args[4]
gausr = as.numeric(args[5])

onlyone=0
if(length(args)==6){
   nts = args[6]  #hidden feature: only save 1 point
   if (nts==1){
      onlyone=1
   }
}
library(ncdf4)
#################################################################################
#										#
# GRACE GSM files represent the full gravity field as Stokes coefficients	# 
#							   Clm(t), Slm(t)	#
# Here we apply the following post-processing needed to obtain the signal in 	#
# equivalent water height:							#
#										#
# Read DS coefficients, replace C20, C11, S11, and C10
# Multiply by Kl and gaussian filter
# Save to output file

nci=nc_open(fnamei)
time<-ncvar_get(nci,'time')
Clm<-ncvar_get(nci,'clm')
Slm<-ncvar_get(nci,'slm')
nc_close(nci)

# Step 1: Read and replace coefficients:
C20<-read.table(c20f)
C11<-read.table(c11f) #C10 C11 S11

tC20<-C20$V1
tC11<-C11$V1

#replace coefficients in C20 (index 3,1)
#plot(time,Clm[3,1,],type='l')
Clm[3,1,tC20]=C20$V2
#lines(time,Clm[3,1,],col='red')

#replace coefficients in C10 (index 2,1)
#plot(time,Clm[2,1,],type='l')  #all zeroes
Clm[2,1,tC11]=C11$V2
#plot(time,Clm[2,1,],type='l')  # zeros do not affect that much

#replace coefficients in C11 (index 2,2)
#plot(time,Clm[2,2,],type='l') #all zeros
Clm[2,2,tC11]=C11$V3
#plot(time,Clm[2,2,],type='l')  #zeros do not affect that much

#replace coefficients in S11 (index 2,2)
#plot(time,Slm[2,2,],type='l') #all zeros
Slm[2,2,tC11]=C11$V3
#plot(time,Slm[2,2,],type='l')  #zeros do not affect that much

#since zeros do not affect the signal significantly, no cropping in time will be done

#Now, we multiply by Kl
a = 6371 #km            Mean equatorial radius of the Earth     
rhoe = 5517 #kg/m3      Average density of the Earth
rhow = 1000 #kg/m3      Density of fresh water

ltrunc=40 #lmax is 40 because of de-striping filter applied
if(gausr==0) ltrunc=60  #if no gaussian 

#Define load love numbers (notice the - sign!!) obtained from Wahr et al 1998
kl_wahr <- -c(0,-.027,0.303,0.194,0.132,0.104, 0.089, 0.081, 0.076, 0.072, 0.069, 0.064, 0.058, 0.051, 0.040, 0.033, 0.027, 0.020, 0.014, 0.01, 0.007)
lk<-c(0,1,2,3,4,5,6,7,8,9,10,12,15,20,30,40,50,70,100,150,200)
out<-approx(lk,kl_wahr,seq(0,ltrunc))   
kl=out$y #load love numbers interpolated
#However, for the geocenter:
kl[2]=0.021 #load love number that needs to be applied to geocenter

#Generate gaussian filter:
if(gausr>0) W = exp(-((0:ltrunc)*gausr/a)^2/(4*log(2)))

#And now, calculate capital K, units will be in mm:
Kl = a*rhoe*(2*seq(0,ltrunc,1)+1)/(3*(1+kl))*1000
if(gausr>0) Kl = Kl*W

CLM=Clm
SLM=Slm
Clm=Clm[1:(ltrunc+1),1:(ltrunc+1),]
Slm=Slm[1:(ltrunc+1),1:(ltrunc+1),]
for (l in 0:ltrunc){
         Clm[l+1,,]=Clm[l+1,,]*Kl[l+1]
         Slm[l+1,,]=Slm[l+1,,]*Kl[l+1]
}

if (onlyone==0){
#we will store the coefficients small netcdf file with dimensions l,m,t
#Define dimensions	
dimt<-ncdim_def('time','GRACE month',1:length(time))
diml<-ncdim_def('degree','l',0:ltrunc)
dimm<-ncdim_def('order','m',0:ltrunc)
	
#Define variables 
ncoef1=ncvar_def('clm','mm',list(diml,dimm,dimt),-9999,longname='delta Clm Temporal mean from 2002120-2015243 removed',prec='float')
ncoef2=ncvar_def('slm','mm',list(diml,dimm,dimt),-9999,longname='delta Slm Temporal mean from 2002120-2015243 removed',prec='float')
	
#Define nc file
ncout<-nc_create(fnameo,list(ncoef1,ncoef2),force_v4=TRUE)

#save day periods and de-trended coefficients
ncvar_put(ncout,'clm',Clm)
ncvar_put(ncout,'slm',Slm)

#Close file
nc_close(ncout)

command<-sprintf('ncdump -h %s',fnameo)
system(command)
}else{
print('Saving only first time step')
#Define dimensions	
dimt<-ncdim_def('time','GRACE month',1)
diml<-ncdim_def('degree','l',0:ltrunc)
dimm<-ncdim_def('order','m',0:ltrunc)
	
#Define variables 
ncoef1=ncvar_def('clm','mm',list(diml,dimm,dimt),-9999,longname='delta Clm Temporal mean from 2002120-2015243 removed',prec='float')
ncoef2=ncvar_def('slm','mm',list(diml,dimm,dimt),-9999,longname='delta Slm Temporal mean from 2002120-2015243 removed',prec='float')
	
#Define nc file
ncout<-nc_create(fnameo,list(ncoef1,ncoef2),force_v4=TRUE)

#save day periods and de-trended coefficients
ncvar_put(ncout,'clm',Clm[,,1])
ncvar_put(ncout,'slm',Slm[,,1])

#Close file
nc_close(ncout)

command<-sprintf('ncdump -h %s',fnameo)
system(command)
}
