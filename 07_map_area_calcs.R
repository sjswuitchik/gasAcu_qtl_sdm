# SM & SJSW 2020
# Calculate areas from the map rasters (in proportions)

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
# Load rasters, add projection
r4NC<-raster("COMBOTE_Nov26_rcp4_noEvol.asc")
proj4string(r4NC)<-CRS("+init=epsg:4326")
r4NW<-raster("COMBO_TOL_E_W_Nov26_rcp4_noEvol.asc")
proj4string(r4NW)<-CRS("+init=epsg:4326")
r4EW<-raster("COMBO_TOL_E_W_Nov3_rcp4_evol.asc")
proj4string(r4EW)<-CRS("+init=epsg:4326")
r8NW<-raster("COMBO_TOL_E_W_Nov3_rcp8_noEvol.asc")
proj4string(r8NW)<-CRS("+init=epsg:4326")
r8EW<-raster("COMBO_TOL_E_W_Nov3_rcp8_evol.asc")
proj4string(r8EW)<-CRS("+init=epsg:4326")
r4EWadjPVE<-raster("COMBO_TOL_E_W_Nov30_rcp4_adjPVE_minOnly.asc")
proj4string(r4EWadjPVE)<-CRS("+init=epsg:4326")
r8EWadjPVE<-raster("COMBO_TOL_E_W_Nov30_rcp8_adjPVE_minOnly.asc")
proj4string(r8EWadjPVE)<-CRS("+init=epsg:4326")
r4EWadjPVEv6<-raster("COMBO_TOL_E_W_Feb9_rcp4_tolerrat_adj.asc")
proj4string(r4EWadjPVEv6)<-CRS("+init=epsg:4326")
r8EWadjPVEv6<-raster("COMBO_TOL_E_W_Feb9_rcp8_tolerrat_adj.asc")
proj4string(r8EWadjPVEv6)<-CRS("+init=epsg:4326") 

# Find the raster values
unique(r4NC)

# Verify what the raster values correspond to in this version (should be consistent but check it anyways)
pdf("raw_r4NC.pdf")
plot(r4NC) 
dev.off()

# Values: 
# 0 = outside of physiological limits
# 2 = within physiological limits
# 3 = normal behaviour

# NB: RCP 4.5 current and RCP 8.5 current are identical (because the current day model is not effected by forecasting scenarios) so for ease of coding, only RCP 4.5 current day is used here

#### Comparison 1: RCP 4.5 no evol vs evol in warmer world (r4EW/r4NW) ####

# Replace the raster components we don't care about with NA (ie/ isolating the 'normal behaviour' envelope)
# NB: extents are the same between rasters, so don't need to change anything
r4NW[r4NW!=3]<-NA
pdf("r4NW.pdf")
plot(r4NW)
dev.off()

r4EW[r4EW!=3]<-NA
pdf("r4EW.pdf")
plot(r4EW)
dev.off()

# Change projection to an equal area projection (useable for areas N of 45 lat)
proj4string(r4NW)<-CRS("+init=epsg:3572") 
proj4string(r4EW)<-CRS("+init=epsg:3572")

# Convert to polygon
poly.r4NW<-rasterToPolygons(r4NW,na.rm=TRUE,dissolve=TRUE)
poly.r4EW<-rasterToPolygons(r4EW,na.rm=TRUE,dissolve=TRUE)

#Difference in the area as a proportion
area(poly.r4EW)/area(poly.r4NW)

#### Comparison 2: RCP 8.5 no evol warmer vs evol warmer (r8EW/r8NW) ####
# NB: extents are the same
r8NW[r8NW!=3]<-NA
pdf("r8NW.pdf")
plot(r8NW)
dev.off()

r8EW[r8EW!=3]<-NA
pdf("r8EW.pdf")
plot(r8EW)
dev.off()

# Change projection to an equal area projection (useable for areas N of 45 lat)
proj4string(r8NW)<-CRS("+init=epsg:3572") 
proj4string(r8EW)<-CRS("+init=epsg:3572")

# Convert to polygon
poly.r8NW<-rasterToPolygons(r8NW,na.rm=TRUE,dissolve=TRUE)
poly.r8EW<-rasterToPolygons(r8EW,na.rm=TRUE,dissolve=TRUE)

# Difference in the area
area(poly.r8EW)/area(poly.r8NW)

#### Comparison 3: RCP 4.5 current vs evol warmer (r4EW/r4NC) #### 

#Extents need to be identical, make the smaller extent match the larger
r4NC<-extend(r4NC,r4NW)

#In the smaller raster, replace the raster components we don't care about with NA
r4NC[r4NC!=3]<-NA
pdf("r4NC.pdf")
plot(r4NC)
dev.off()

#Change projection to an equal area projection (useable for areas N of 45 lat)
proj4string(r4NC)<-CRS("+init=epsg:3572") 

# Convert to polygon
poly.r4NC<-rasterToPolygons(r4NC,na.rm=TRUE,dissolve=TRUE)

#Difference in the area
area(poly.r4EW)/area(poly.r4NC)

#### Comparison 4: RCP 8.5 current vs evol warmer (r8EW/r4NC) ####
# NB: all polys have been created for this comparison

# Difference in the area
area(poly.r8EW)/area(poly.r4NC)

#### Comparison 5: RCP 4.5 current vs no evol warmer (r4NW/r4NC) ####
area(poly.r4NW)/area(poly.r4NC)

#### Comparison 6: RCP 8.5 current vs no evol warmer (r8NW/r4NC) ####
area(poly.r8NW)/area(poly.r4NC)

#### Comparison 7: RCP 4.5 evol warmer vs adj PVE evol warmer (r4EW/r4EWadjPVE) ####

# process adjPVE rasters
r4EWadjPVE[r4EWadjPVE!=3]<-NA
pdf("r4EWadjPVE.pdf")
plot(r4EWadjPVE)
dev.off()

r8EWadjPVE[r8EWadjPVE!=3]<-NA
pdf("r8EWadjPVE.pdf")
plot(r8EWadjPVE)
dev.off()

# Change projection to an equal area projection (useable for areas N of 45 lat)
proj4string(r4EWadjPVE)<-CRS("+init=epsg:3572") 
proj4string(r8EWadjPVE)<-CRS("+init=epsg:3572")

# Convert to polygon
poly.r4EWadjPVE<-rasterToPolygons(r4EWadjPVE,na.rm=TRUE,dissolve=TRUE)
poly.r8EWadjPVE<-rasterToPolygons(r8EWadjPVE,na.rm=TRUE,dissolve=TRUE)

#Difference in the area as a proportion
area(poly.r4EW)/area(poly.r4EWadjPVE)

#### Comparsion 8: RCP 8.5 evol warmer vs adj PVE evol warmer (r8EW/r8EWadjPVE) ####
area(poly.r8EW)/area(poly.r8EWadjPVE)

#### Comparison 9: RCP 4.5 current vs adj PVE evol warmer (r4EWadjPVE/r4NC) ####
# NB: need to run whole script to get the properly extent-matched r4NC here
area(poly.r4EWadjPVE)/area(poly.r4NC)

#### Comparison 10: RCP 8.5 current vs adj PVE evol warmer (r8EWadjPVE/r4NC) ####
area(poly.r8EWadjPVE)/area(poly.r4NC)

#### Comparison 11: RCP 4.5 current vs adj PVE evol warmer v6 (CTmin+CTmax evol) (r4EWadjPVEv6/r4NC) ####
r4EWadjPVEv6[r4EWadjPVEv6!=3]<-NA
pdf("r4EWadjPVEv6.pdf")
plot(r4EWadjPVEv6)
dev.off()

r8EWadjPVEv6[r8EWadjPVEv6!=3]<-NA
pdf("r8EWadjPVEv6.pdf")
plot(r8EWadjPVEv6)
dev.off()

# Change projection to an equal area projection (useable for areas N of 45 lat)
proj4string(r4EWadjPVEv6)<-CRS("+init=epsg:3572") 
proj4string(r8EWadjPVEv6)<-CRS("+init=epsg:3572")

# Convert to polygon
poly.r4EWadjPVEv6<-rasterToPolygons(r4EWadjPVEv6,na.rm=TRUE,dissolve=TRUE)
poly.r8EWadjPVEv6<-rasterToPolygons(r8EWadjPVEv6,na.rm=TRUE,dissolve=TRUE)

area(poly.r4EWadjPVEv6)/area(poly.r4NC)

#### Comparison 12: RCP 8.5 current vs adj PVE evol warmer v6 (CTmin+CTmax evol) (r8EWadjPVEv6/r4NC) ####
area(poly.r8EWadjPVEv6)/area(poly.r4NC)

