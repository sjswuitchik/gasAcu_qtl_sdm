# SM & SJSW 2019-21
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

#Load raster of combined tolerance and erratic behaviour
COMBOTE<-raster("outputs/Processed_files/COMBOTE_jan2022_rcp4_noEvol.asc")
COMBOTEWARM<-raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp4_noEvol.asc")

#Load bathymetry to get coastlines (nc file)
Bathy<-nc_open("required_files/GEBCO_2014_2D_-179.7777_45.2207_-120.8539_76.3978.nc")

#Load shapefile of bathymetry/sea ice availability Current World
suit<-readOGR("required_files/Bathy&SeaIce-prefered")

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
A<-raster("required_files/SaraUseThis.tif")

#For current bathy and sea ice
#Dissolve polygons together, to form one polygon of possible habitats
#new<-unionSpatialPolygons(suit,rep(1,length(suit$Id)))

#Plot-Current World
viridis <- viridis_pal(direction = 1, option = "C")
cols <- viridisLite::viridis(3)
ltext<-c("Outside Physiol Limits","Within Physiol Limits","Normal Behav") #legend text
pdf("outputs/Figs/rcp4_noevol_current.pdf")
plot(D, col="grey",axes=F,legend=F)
plot(COMBOTE, add=T, legend=F, at=c(0,2,3),col=cols[1:3])
legend("bottomleft",legend=ltext,fill=cols,bg="white")
dev.off()

#Plot-Warmer World
pdf("outputs/Figs/rcp4_noevol_warmer.pdf")
plot(D, col="grey",axes=F,legend=F)
plot(COMBOTEWARM, add=T, legend=F, at=c(2,3),col=cols[2:3])
legend("bottomleft",legend=ltext,fill=cols,bg="white")
dev.off()

# Plot - Bathy/sea ice current + warmer
cols <- viridisLite::viridis(2)
lstext <- c("Suitable Habitat Warmer", "Suitable Habitat Current")
pdf("outputs/Figs/suitable_combo.pdf")
plot(D,col="grey",axes=F,legend=F)
plot(B,add=T,col=cols[1],axes=F,legend=F)
plot(A,add=T,col=cols[2],axes=F)
legend("bottomleft",legend=lstext,fill=cols[1:2],bg="white")
dev.off()
