# SM & SJSW 2019-20
#Make Figures for the distribution data

###########################################################################################################################

rm(list=ls(all=TRUE))

library(raster)
library(rgeos)
library(sp)
library(rgdal)
library(maptools)
library(ncdf4)
library(viridis)

###########################################################################################################################
setwd("~/Desktop/UCalgary/Main publication/Coding")

#Load raster of combined tolerance and erratic behaviour
COMBOTE<-raster("Processed_files/COMBOTE_Jan13_rcp4_noevol.asc")
COMBOTEWARM<-raster("Processed_files/COMBO_TOL_E_W_Jan13_rcp4_noevol.asc")

#Load bathymetry to get coastlines (nc file)
Bathy<-nc_open("GEBCO_2014_2D_-179.7777_45.2207_-120.8539_76.3978.nc")

#Load shapefile of bathymetry/sea ice availibility Current World
suit<-readOGR("Bathy&SeaIce-prefered")

#Find lat and lon of cells
lon<-ncvar_get(Bathy,"lon")
lat<-ncvar_get(Bathy,"lat")

#Get depth/elevation data
Depths<-ncvar_get(Bathy,"elevation")
#Transpose Depth matrix so map will be right way up
Depths <- t(Depths)
#Create raster
D<-raster(Depths,xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))
#Flip raster so map will be right way up
D<-flip(D,direction='y')
#Create raster for bathy where stickelback cannot persist
B<-D
#Alter raster for plotting
D[D<0]<-NA #Remove all cells below sea level (i.e. land remains)
D[D>=0]<-1 #Give all cells above sea level the same value

#Create cells for bathy where stickelback cannot persist
B[B>0]<-NA #Remove all cells above sea level (i.e. on sea remains)
B[B<=-200]<-NA #Remove cell greater than 200m depth
B[B<=0 & B>-200]<-1 #Set values in raster to all the same value

#Set Northern extent to lat of Wales, Alaska (65.6 N)
B[1:1253,]<-NA #Set cells above the Northern extent to 0
A<-raster("SaraUseThis.tif")

#For current bathy and sea ice
#Dissolve polygons together, to form one polygon of possible habitats
#new<-unionSpatialPolygons(suit,rep(1,length(suit$Id)))

#Plot-Current World
viridis <- viridis_pal(direction = 1, option = "C")
cols <- viridisLite::viridis(5)
ltext<-c("Outside Physiol Limits","Within Physiol Limits","Normal Behav") #legend text
pdf("Figs/rcp4_noevol_current.pdf")
plot(D, col="grey",axes=F,legend=F)
plot(COMBOTE, add=T, legend=F, at=c(0,2,3),col=cols[1:3])
legend("bottomleft",legend=ltext,fill=cols,bg="white")
dev.off()

#Plot-Warmer World
pdf("Figs/rcp4_noevol_warmer.pdf")
plot(D, col="grey",axes=F,legend=F)
plot(COMBOTEWARM, add=T, legend=F, at=c(2,3),col=cols[2:3])
legend("bottomleft",legend=ltext,fill=cols,bg="white")
dev.off()

# Plot - Bathy/sea ice current + warmer
lstext <- c("Suitable Habitat Warmer", "Suitable Habitat Current")
pdf("Figs/suitable_combo.pdf")
plot(D,col="grey",axes=F,legend=F)
plot(B,add=T,col=cols[4],axes=F,legend=F)
plot(A,add=T,col=cols[5],axes=F)
legend("bottomleft",legend=lstext,fill=cols[4:5],bg="white")
dev.off()


## Create suitable habitat rasters for area calculations
# check that the rasters are properly plotted without the NA continent
pdf("suitable_current.pdf")
plot(A)
dev.off()
pdf("suitable_warmer.pdf")
plot(B)
dev.off()

# Extents need to be identical, make the smaller extent match the larger
A <- extend(A,B)

# Change projection to an equal area projection (useable for areas N of 45 lat) 
proj4string(B) <- CRS("+init=epsg:3572")
proj4string(A) <- CRS("+init=epsg:3572")

# Convert to polygons
poly.A <- rasterToPolygons(A,na.rm=TRUE,dissolve=TRUE)
poly.B <- rasterToPolygons(B,na.rm=TRUE,dissolve=TRUE)

# Calculate the proportional difference in area
area(poly.B)/area(poly.A)

