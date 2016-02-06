library(maptools)
library(raster)
library(geosphere)
#Function to load either a shape file or text file and pass it as a Spatial Polygon:
basin_load<-function(fname){
#Check fname for .shp:
extension=unlist(strsplit(fname, "[.]"))
extension=extension[length(extension)]
if(extension=='shp'){
#Assuming shape file
cat(sprintf('Opening shape file %s\n',fname));
out<-readShapeSpatial(dirGRDC)
}else{
#Assuming two-column formatted text file
cat(sprintf('Loading text file %s\n',fname));
pfoo<-unname(as.matrix(read.table(fname),ncol=2))
pfoo<-Polygon(pfoo)
pfoo<-Polygons(list(pfoo),1)
out<-SpatialPolygons(list(pfoo))
}
}

# Calculate angular region (area):
#angular_area = areaPolygon(out,a=1)
#return(list(ang=angular_area,ext=extent(out)))
