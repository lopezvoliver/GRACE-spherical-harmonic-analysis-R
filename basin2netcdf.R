#!/usr/bin/env Rscript
#Usage: ./basin2netcdf.R shape_file out.nc
#Load basin file 
global = 0 #if 1, global, if not, crop to extent of basin
args<-commandArgs(trailingOnly=TRUE)
if(length(args)==0){
 print('Usage: ./basin2netcdf.R shape_file out.nc')
 print('shape_file: name of shape file (.shp) or text file containing polygon coordinates (.txt or .dat)')
 print('out.nc: name of output netCDF file')
 quit()
}
name=args[1]
fnameo=args[2]

source('loadbasin.R')
out<-basin_load(name)

# Extent of rasterization and resolution defined here:
lon<-seq(-179.75,179.75,0.5)
lat<-seq(-89.75,89.75,0.5)

library(raster) #raster functions
library(rgdal)  #readOGR

# The masking function:
maskregion<-function(lon,lat,maskdf){

     if (class(maskdf)=='SpatialPolygonsDataFrame'| class(maskdf)=='SpatialPolygons'){
            #lon must be in -180:180 format
            if(max(lon)>180) lon=lon-180
    r<-raster(ncols=length(lon),nrows=length(lat))
        extent(r)<-extent(lon[1],lon[length(lon)],lat[1],lat[length(lat)])
        land_raster<-rasterize(maskdf,r)
            lmask<-matrix(land_raster@data@values,ncol=length(lat),nrow=length(lon),byrow=FALSE)
            lmask[!is.na(lmask)]=1
                lmask=lmask[,length(lat):1]
                return(lmask)
                  }else{
                         ttemp=expand.grid(lon,lat)
                    t=array(NA,dim=c(length(lon)*length(lat),2))
                        t[,1]=ttemp$Var1
                        t[,2]=ttemp$Var2
                            inbas=pnt.in.poly(matrix(t,ncol=2),maskdf)$pip
                            infoo<-matrix(inbas,ncol=length(lat),nrow=length(lon))
                                infoo[infoo==0]=NA
                                return(infoo)
                                  }

}
ang<-areaPolygon(out,a=1) #needed for calculating time series
basin_raster<-maskregion(lon,lat,out)

#Transform to 0:360 (only if global) and flip lat to be consistent with dragon
if (global==1){
res=0.5
lon<-lon+180
lat<-lat[length(lat):1]
temp<-basin_raster[1:(180/res),]
basin_raster[1:(180/res),]=basin_raster[(1+180/res):(360/res),]
basin_raster[(1+180/res):(360/res),]=temp
basin_raster=t(basin_raster)
basin_raster=basin_raster[length(lat):1,]
}else{
#crop at extent of basin
   ext<-extent(out)
   print(ext)

   x1=which(lon<xmin(ext)); x1=x1[length(x1)]; #print(x1);
   x2=which(lon>xmax(ext));x2=x2[1]; #print(x2);
   y1=which(lat<ymin(ext)); y1=y1[length(y1)]; #print(y1);
   y2=which(lat>ymax(ext)); y2=y2[1];#print(y2);

   lon<-lon[x1:x2]
   lat<-lat[y1:y2]
   basin_raster=basin_raster[x1:x2,y1:y2]
} 

#Output to netcdf file:
library(ncdf4)
dlon<-ncdim_def("lon","degrees east",lon)
dlat<-ncdim_def("lat","degrees north",lat)
dimt<-ncdim_def("time","time foo",1)

#Define variables
var1out=ncvar_def('mask','Basin mask',list(dlon,dlat,dimt),-9999,longname='Basin mask (raster)',prec='float')
var2out=ncvar_def('ang','Angular region',list(),-9999,longname='Angular region of basin', prec='float')

#Create file
ncout<-nc_create(fnameo,list(var1out,var2out),force_v4=TRUE)
#save variables
ncvar_put(ncout,'mask',basin_raster)
ncvar_put(ncout,'ang',ang)
nc_close(ncout)
