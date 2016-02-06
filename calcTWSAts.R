#!/usr/bin/env Rscript
#Usage: ./calcTWSAts.R basin_file GRACE_SHfile.nc basin_SHfile.nc output.txt scale_factor
args<-commandArgs(trailingOnly = TRUE)
library(ncdf4)

if(length(args)<4){
   cat('Usage: \n')
   cat('./calcTWSAts.R shapefile GRACE_SHfile.nc shape_SHfile.nc output.txt [scale_factor]\n')
   cat('shapefile: Shape file (.shp) or text file (.dat or .txt) containing region polygon\n')
   cat('GRACE_SHfile.nc: name of GRACE TWSA SH coefficients file \n')
   cat('shape_SHfile.nc: name of basin coef file \n')
   cat('output.txt: name of output file (txt) \n')
   cat('scale_factor: Optional scaling factor (Defaut =1)\n')
      quit()
}

scale=1
gia=0
if(length(args)>=5){
   scale=as.numeric(args[5])
   #Secret option: remove GIA trend
   if(length(args)>5){
      gia=as.numeric(args[6])
   }
}

basname = args[1] # name or source file of polygon 
fnamei = args[2] #name of file containing coefficients
fnameb = args[3] #name of basin coefficients file
fnameo = args[4] #name of output file

#Load basin file 
source('loadbasin.R')
bpol<-basin_load(basname)

#Calculate angular region of polygon
ang = areaPolygon(bpol,a=1)

# Read coefficients from GRACE file (clm, slm)
nci<-nc_open(fnamei)
clm<-ncvar_get(nci,'clm')
slm<-ncvar_get(nci,'slm')
time<-ncvar_get(nci,'time')
nc_close(nci)

# Read coefficients from file (basin coefficients)
nci<-nc_open(fnameb)
bflmr<-ncvar_get(nci,'clm')  #We use these ones because they are geodetically normalized
bflmi<-ncvar_get(nci,'slm')  #see swenson and wahr. In this case, they ARE the same as the ones in dragon
nc_close(nci)

lmax=dim(clm)[1]-1

bClm<-bflmr
bSlm<-bflmi

#for(l in 0:lmax){
#   for (m in 0:l){
#      if(m==0){
#      bClm[l+1,m+1] = bflmr[l+1,m+1]
#      bSlm[l+1,m+1] = 0
#      }else{
#      bClm[l+1,m+1] = bflmr[l+1,m+1]*(-1)^(m)
#      bSlm[l+1,m+1] = bflmi[l+1,m+1]*(-1)^(m+1)
#      }
#   }
#}
# Get GRACE data information
gdate<-read.table('GRACE_timevector.txt')

if(length(time)==1){
   cat(signal)
   quit()
}
timedate<-as.Date(gdate$V1)
deltaT = as.numeric(timedate - timedate[1]) #TIme date since reference, in days
timedate<-as.character(timedate)


#Signal reconstructed from global * basin coefficients:
signal = array(NA,dim=length(time))
for (k in 1:length(time)){
if(length(time)==1){
signal[k] = sum(bClm[1:(lmax+1),1:(lmax+1)]*clm+bSlm[1:(lmax+1),1:(lmax+1)]*slm,na.rm=TRUE)/ang
}else{
signal[k] = sum(bClm[1:(lmax+1),1:(lmax+1)]*clm[,,k]+bSlm[1:(lmax+1),1:(lmax+1)]*slm[,,k],na.rm=TRUE)/ang
}
}
signal=signal*4*pi
signal = signal - (gia/365)*deltaT# Remove GIA trend
signal=signal*scale

#Re-center signal in y (remove mean from whole period)
signal = signal - mean(signal)

# Finally we output the data (three columns)
#OUT=cbind(timedate,signal)
OUT=list(date=timedate,y=signal)
write.table(file=fnameo,OUT,row.names=FALSE,col.names=TRUE,quote=FALSE)
#plot(time,osignal,type='l')
#lines(time,signal,col='blue')
#lines(time,signal2,col='red')
