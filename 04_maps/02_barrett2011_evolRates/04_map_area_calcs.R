# SM & SJSW 2020-22
# Calculate areas from the map rasters 

###########################################################################################################################

rm(list=ls(all=TRUE))

library(raster)
library(rgeos)
library(sp)
library(rgdal)
library(maptools)
library(ncdf4)
library(viridis)

############################################################################################################
# Load rasters, add projection
r4NC<-raster("outputs/Processed_files/COMBOTE_jan2022_rcp4_noEvol.asc")
proj4string(r4NC)<-CRS("+init=epsg:4326")
r4NW<-raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp4_noEvol.asc")
proj4string(r4NW)<-CRS("+init=epsg:4326")
r4EW<-raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp4_lowerEvol.asc")
proj4string(r4EW)<-CRS("+init=epsg:4326")
r8NW<-raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp8_noEvol.asc")
proj4string(r8NW)<-CRS("+init=epsg:4326")
r8EW<-raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp8_lowerEvol.asc")
proj4string(r8EW)<-CRS("+init=epsg:4326")
r4min3<-raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp4_lowerEvol.asc") # rasters are identical
proj4string(r4min3)<-CRS("+init=epsg:4326")
r4min21<-raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp4_lowerEvol.asc") # rasters are identical
proj4string(r4min21)<-CRS("+init=epsg:4326")
r4max20<-raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp4_adjPVE_maxLG20.asc")
proj4string(r4max20)<-CRS("+init=epsg:4326")
r4max12<-raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp4_adjPVE_maxLG12.asc")
proj4string(r4max12)<-CRS("+init=epsg:4326")
r8min3<-raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp8_lowerEvol.asc") # rasters are identical
proj4string(r8min3)<-CRS("+init=epsg:4326")
r8min21<-raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp8_lowerEvol.asc") # rasters are identical
proj4string(r8min21)<-CRS("+init=epsg:4326")
r8max20<-raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp8_adjPVE_maxLG12.asc") # rasters are identical
proj4string(r8max20)<-CRS("+init=epsg:4326")
r8max12<-raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp8_adjPVE_maxLG12.asc") # rasters are identical 
proj4string(r8max12)<-CRS("+init=epsg:4326")
r4max<-raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp4_upperEvol.asc")
proj4string(r4max)<-CRS("+init=epsg:4326")
r8max<-raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp8_upperEvol.asc")
proj4string(r8max)<-CRS("+init=epsg:4326")

## New rasters for Dec revisions
r4minMaxEW <- raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp4_lowerUpperEvol.asc")
proj4string(r4minMaxEW) <- CRS("+init=epsg:4326")
r8minMaxEW <- raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp8_lowerUpperEvol.asc")
proj4string(r8minMaxEW) <- CRS("+init=epsg:4326")

# Find the raster values
unique(r4NC)

# Verify what the raster values correspond to in this version (should be consistent but check it anyways)
pdf("outputs/Figs/raw_r4NC.pdf")
plot(r4NC) 
dev.off()

# Values: 
# 0 = outside of physiological limits
# 2 = within physiological limits
# 3 = normal behaviour

# NB: RCP 4.5 current and RCP 8.5 current are identical (because the current day model is not effected by forecasting scenarios) so for ease of coding, only RCP 4.5 current day is used here

## order of operations
# Replace the raster components we don't care about with NA (ie/ isolating the 'normal behaviour' envelope)
# change to equal area projection (usable for areas N of 45 lat)
# convert to polygon
# calculate area of polygon

r4NW[r4NW!=3]<-NA
pdf("outputs/Figs/r4NW.pdf")
plot(r4NW)
dev.off()

r4EW[r4EW!=3]<-NA
pdf("outputs/Figs/r4EW.pdf")
plot(r4EW)
dev.off()

# Change projection to an equal area projection (useable for areas N of 45 lat)
proj4string(r4NW)<-CRS("+init=epsg:3572") 
proj4string(r4EW)<-CRS("+init=epsg:3572")

# Convert to polygon
poly.r4NW<-rasterToPolygons(r4NW,na.rm=TRUE,dissolve=TRUE)
poly.r4EW<-rasterToPolygons(r4EW,na.rm=TRUE,dissolve=TRUE)

# Areas 
area(poly.r4EW) # 193.9766
area(poly.r4NW) # 152.5273

r8NW[r8NW!=3]<-NA
pdf("outputs/Figs/r8NW.pdf")
plot(r8NW)
dev.off()

r8EW[r8EW!=3]<-NA
pdf("outputs/Figs/r8EW.pdf")
plot(r8EW)
dev.off()

proj4string(r8NW)<-CRS("+init=epsg:3572") 
proj4string(r8EW)<-CRS("+init=epsg:3572")

poly.r8NW<-rasterToPolygons(r8NW,na.rm=TRUE,dissolve=TRUE)
poly.r8EW<-rasterToPolygons(r8EW,na.rm=TRUE,dissolve=TRUE)

area(poly.r8EW) # 193.9766
area(poly.r8NW) # 173.2928


r4NC[r4NC!=3]<-NA
pdf("outputs/Figs/r4NC.pdf")
plot(r4NC)
dev.off()

#Change projection to an equal area projection (useable for areas N of 45 lat)
proj4string(r4NC)<-CRS("+init=epsg:3572") 

# Convert to polygon
poly.r4NC<-rasterToPolygons(r4NC,na.rm=TRUE,dissolve=TRUE)

# Area
area(poly.r4NC) # 25.10572

#### Areas with new QTL data - July 2021 ####
r4min3[r4min3!=3]<-NA
pdf("outputs/Figs/r4min3.pdf")
plot(r4min3)
dev.off()

r4min21[r4min21!=3]<-NA
pdf("outputs/Figs/r4min21.pdf")
plot(r4min21)
dev.off()

r4max20[r4max20!=3]<-NA
pdf("outputs/Figs/r4max20.pdf")
plot(r4max20)
dev.off()

r4max12[r4max12!=3]<-NA
pdf("outputs/Figs/r4max12.pdf")
plot(r4max12)
dev.off()

r8min3[r8min3!=3]<-NA
pdf("outputs/Figs/r8min3.pdf")
plot(r8min3)
dev.off()

r8min21[r8min21!=3]<-NA
pdf("outputs/Figs/r8min21.pdf")
plot(r8min21)
dev.off()

r8max20[r8max20!=3]<-NA
pdf("outputs/Figs/r8max20.pdf")
plot(r8max20)
dev.off()

r8max12[r8max12!=3]<-NA
pdf("r8max12.pdf")
plot(r8max12)
dev.off()

# Change projection to an equal area projection (useable for areas N of 45 lat)
proj4string(r4min3)<-CRS("+init=epsg:3572") 
proj4string(r4min21)<-CRS("+init=epsg:3572")
proj4string(r4max20)<-CRS("+init=epsg:3572")
proj4string(r4max12)<-CRS("+init=epsg:3572")
proj4string(r8min3)<-CRS("+init=epsg:3572") 
proj4string(r8min21)<-CRS("+init=epsg:3572")
proj4string(r8max20)<-CRS("+init=epsg:3572")
proj4string(r8max12)<-CRS("+init=epsg:3572")

# Convert to polygon
poly.r4min3<-rasterToPolygons(r4min3,na.rm=TRUE,dissolve=TRUE)
poly.r4min21<-rasterToPolygons(r4min21,na.rm=TRUE,dissolve=TRUE)
poly.r4max20<-rasterToPolygons(r4max20,na.rm=TRUE,dissolve=TRUE)
poly.r4max12<-rasterToPolygons(r4max12,na.rm=TRUE,dissolve=TRUE)
poly.r8min3<-rasterToPolygons(r8min3,na.rm=TRUE,dissolve=TRUE)
poly.r8min21<-rasterToPolygons(r8min21,na.rm=TRUE,dissolve=TRUE)
poly.r8max20<-rasterToPolygons(r8max20,na.rm=TRUE,dissolve=TRUE)
poly.r8max12<-rasterToPolygons(r8max12,na.rm=TRUE,dissolve=TRUE)

# Area calcs
area(poly.r4min3) # 193.9766
area(poly.r4min21) # 193.9766
area(poly.r4max20) # 152.5273
area(poly.r4max12) # 152.5273
area(poly.r8min3) # 193.9766
area(poly.r8min21) # 193.9766
area(poly.r8max20) # 173.2928
area(poly.r8max12) # 173.2928

## New models for Dec 2021 revisions
# r4minMaxEW
# r8minMaxEW

r4minMaxEW[r4minMaxEW!=3]<-NA
pdf("outputs/Figs/r4minMaxEW.pdf")
plot(r4minMaxEW)
dev.off()

proj4string(r4minMaxEW)<-CRS("+init=epsg:3572") 

poly.r4minMaxEW<-rasterToPolygons(r4minMaxEW,na.rm=TRUE,dissolve=TRUE)

r8minMaxEW[r8minMaxEW!=3]<-NA
pdf("outputs/Figs/r8minMaxEW.pdf")
plot(r8minMaxEW)
dev.off()

proj4string(r8minMaxEW)<-CRS("+init=epsg:3572") 

poly.r8minMaxEW<-rasterToPolygons(r8minMaxEW,na.rm=TRUE,dissolve=TRUE)

area(poly.r4minMaxEW) # 193.9766
area(poly.r8minMaxEW) # 193.9766

## New rasters for Jan 2022 revisions
r4maxW <- raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp4_UpperEvol.asc")
proj4string(r4maxW) <- CRS("+init=epsg:4326")
r8maxW <- raster("outputs/Processed_files/COMBO_TOL_E_W_jan2022_rcp8_UpperEvol.asc")
proj4string(r8maxW) <- CRS("+init=epsg:4326")

r4maxW[r4maxW!=3]<-NA
r8maxW[r8maxW!=3]<-NA

pdf("outputs/Figs/r4maxW.pdf")
plot(r4maxW)
dev.off()
pdf("outputs/Figs/r8maxW.pdf")
plot(r8maxW)
dev.off()

proj4string(r4maxW)<-CRS("+init=epsg:3572") 
proj4string(r8maxW)<-CRS("+init=epsg:3572") 

poly.r4maxW<-rasterToPolygons(r4maxW,na.rm=TRUE,dissolve=TRUE)
poly.r8maxW<-rasterToPolygons(r8maxW,na.rm=TRUE,dissolve=TRUE)

area(poly.r4maxW) # 152.5273
area(poly.r8maxW) # 173.2928