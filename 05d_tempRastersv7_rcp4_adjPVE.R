library(ncdf4)
library(raster)
library(rgeos)
library(sp)
library(rgdal)
library(maptools)

SST<-nc_open("Temp2014.nc")

#Save file details to a text file (comment in to view)
#{
 #sink('Temp2014.nc.txt')
 #print(SST)
 #sink()
#}

lon<-ncvar_get(SST,"lon")
lat<-ncvar_get(SST,"lat")
t<-ncvar_get(SST,"time")

Temp.array<-ncvar_get(SST,"sst")
dim(Temp.array)

#Set temp tolerance, erratic behaviour, and preferences

tol.lb<-0.9
tol.ub<-31.9
errat.lb<-5
errat.ub<-25
pref.lb<-13.6
pref.ub<-19.9

#Warmer world - evolution of traits
evol.tol.lb <- tol.lb - 0.9
evol.tol.ub <- tol.ub 
evol.errat.lb <- errat.lb - 1.4
evol.errat.ub <- errat.ub 
evol.pref.lb <- pref.lb
evol.pref.ub <- pref.ub

#Process Temperature data, to get to min temp, max temp, min of 25% percentile, max of 25% percentile, and median
#Set up output matrices

MIN<-matrix(-1000,nrow=dim(Temp.array)[1],ncol=dim(Temp.array)[2])
MAX<-matrix(-1000,nrow=dim(Temp.array)[1],ncol=dim(Temp.array)[2])
MIN25<-matrix(-1000,nrow=dim(Temp.array)[1],ncol=dim(Temp.array)[2])
MAX25<-matrix(-1000,nrow=dim(Temp.array)[1],ncol=dim(Temp.array)[2])
MEDIAN<-matrix(-1000,nrow=dim(Temp.array)[1],ncol=dim(Temp.array)[2])

#Process the temp data: Find min, max, median value in each location

for (jj in 1:dim(Temp.array)[1])
{
	for (ii in 1:dim(Temp.array)[2])	
	{
	X<-lon[jj]
	Y<-lat[ii]	
	MIN[jj,ii]<-min(Temp.array[jj,ii,])
	MAX[jj,ii]<-max(Temp.array[jj,ii,])
	MEDIAN[jj,ii]<-median(Temp.array[jj,ii,])
	}
}

#Find average of upper and lower quartile in each location (364 days of the year mean lower 25 percentile is from days 1-92, upper 25 percentile is from days 273-364)
for (jj in 1:dim(Temp.array)[1])
{
	for (ii in 1:dim(Temp.array)[2])	
	{
	MIN25[jj,ii]<-mean(Temp.array[jj,ii,1:92])
	MAX25[jj,ii]<-mean(Temp.array[jj,ii,273:364])
	}
}

### Create raster
r.min<-raster(t(MIN),xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))
r.min <- flip(r.min, direction='y')
r.min<-rotate(r.min) #Change lat to go go from -180 to 180, not 0 to 360

r.max<-raster(t(MAX),xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))
r.max <- flip(r.max, direction='y')
r.max<-rotate(r.max) #Change lat to go go from -180 to 180, not 0 to 360

r.med<-raster(t(MEDIAN),xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))
r.med <- flip(r.med, direction='y')
r.med<-rotate(r.med) #Change lat to go go from -180 to 180, not 0 to 360

r.min25<-raster(t(MIN25),xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))
r.min25 <- flip(r.min25, direction='y')
r.min25<-rotate(r.min25) #Change lat to go go from -180 to 180, not 0 to 360

r.max25<-raster(t(MAX25),xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))
r.max25 <- flip(r.max25, direction='y')
r.max25<-rotate(r.max25) #Change lat to go go from -180 to 180, not 0 to 360

#Load shapefile of bathymetry/sea ice availibility
suit<-readOGR("Bathy&SeaIce-prefered")

#Dissolve polygons together, to form one polygon of possible habitats
new<-suit

#Reduce size of temperature raster cells, to improve coverage after masking
#Factor by which to increase the number of cells in the raster (x and y dimensions)
ff<-c(20,20)
r.minD<-disaggregate(r.min,fact=ff)
r.maxD<-disaggregate(r.max,fact=ff)
r.medD<-disaggregate(r.med,fact=ff)
r.min25D<-disaggregate(r.min25,fact=ff)
r.max25D<-disaggregate(r.max25,fact=ff)

#Clip temperature rasters based on Bathy and Sea Ice extent
#Mask min water temp, with bathymetry and sea ice extent
temp.min<-mask(r.minD,new)
temp.min<-crop(temp.min,extent(new))

temp.max<-mask(r.maxD,new) 
temp.max<-crop(temp.max,extent(new))

temp.med<-mask(r.medD,new) 
temp.med<-crop(temp.med,extent(new))

temp.min25<-mask(r.min25D,new)
temp.min25<-crop(temp.min25,extent(new))

temp.max25<-mask(r.max25D,new) 
temp.max25<-crop(temp.max25,extent(new))

########################################################################################################################

#Make the warmer world

#Load bathymetry (nc file)
Bathy<-nc_open("GEBCO_2014_2D_-179.7777_45.2207_-120.8539_76.3978.nc")

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
#Remove cells for depths where stickelback cannot persist
#Remove all cells above sea level (i.e. on sea remains)
D[D>0]<-NA 
#Remove cell greater than 200m depth
D[D<=-200]<-NA 
#Set values in raster to all the same value
D[D<=0 & D>-200]<-1 

#Set Northern extent to lat of Wales, Alaska (65.6 N)
#Set cells above the Northern extent to 0
D[1:1253,]<-NA 
#Plot to confirm this worked
plot(D)

#For saving, commented out to avoid accidental overwrites
#writeRaster(D,"path/Processed Files/BathyWarm.asc",format="ascii")
#This file is converted into a polygon in arcgis. 

#Read the polygon in from ArcGis
bathywarm<-readOGR("Warmer Bathy")

#Mask temperature rasters with bathymetry extent for warmer world
temp.min25W<-mask(r.max25D,bathywarm) 
temp.min25W<-crop(temp.min25W,extent(bathywarm))

temp.max25W<-mask(r.max25D,bathywarm) 
temp.max25W<-crop(temp.max25W,extent(bathywarm))


#Warm the waters, based on location

#Winter (based on temp.min25)-need to set up with lat/lon

exwin1<-c(extent(temp.min25W)[1],extent(temp.min25W)[2],62.22,extent(temp.min25W)[4]) #extent of section 1 warming (Northern)
exwin2<-c(extent(temp.min25W)[1],extent(temp.min25W)[2],60.33,62.22) #extent of section 2 warming
exwin3<-c(extent(temp.min25W)[1],extent(temp.min25W)[2],48.88,60.33) #extent of section 3 warming 
exwin4<-c(extent(temp.min25W)[1],extent(temp.min25W)[2],extent(temp.min25W)[3],48.88) #extent of section 4 warming (S of 48.88)
CellsWin1<-cellsFromExtent(temp.min25W,exwin1) #ID cells for each change
CellsWin2<-cellsFromExtent(temp.min25W,exwin2) 
CellsWin3<-cellsFromExtent(temp.min25W,exwin3) 
CellsWin4<-cellsFromExtent(temp.min25W,exwin4) 

exwin5<-c(extent(temp.max25W)[1],extent(temp.max25W)[2],62.22,extent(temp.min25W)[4]) #extent of section 1 warming (Northern)
exwin6<-c(extent(temp.max25W)[1],extent(temp.max25W)[2],60.33,62.22) #extent of section 2 warming
exwin7<-c(extent(temp.max25W)[1],extent(temp.max25W)[2],48.88,60.33) #extent of section 3 warming 
exwin8<-c(extent(temp.max25W)[1],extent(temp.max25W)[2],extent(temp.max25W)[3],48.88) #extent of section 4 warming (S of 48.88)
CellsWin5<-cellsFromExtent(temp.max25W,exwin1) #ID cells for each change
CellsWin6<-cellsFromExtent(temp.max25W,exwin2) 
CellsWin7<-cellsFromExtent(temp.max25W,exwin3) 
CellsWin8<-cellsFromExtent(temp.max25W,exwin4) 

#Warm the min temps
#Create new raster for warmed
temp.min25WARMED<-temp.min25W
#Warm it by relevant amounts for RCP 4.5, for the diff. extents
temp.min25WARMED[CellsWin1]<-temp.min25W[CellsWin1]+2 
temp.min25WARMED[CellsWin2]<-temp.min25W[CellsWin2]+2.4 
temp.min25WARMED[CellsWin3]<-temp.min25W[CellsWin3]+1.8 
temp.min25WARMED[CellsWin4]<-temp.min25W[CellsWin4]+1.6 

#Warm the max temps
#Create new raster
temp.max25WARMED<-temp.max25W 
#Warm it by relevant amounts for RCP 4.5, for the diff. extents
temp.max25WARMED[CellsWin5]<-temp.max25W[CellsWin5]+2 
temp.max25WARMED[CellsWin6]<-temp.max25W[CellsWin6]+2.4 
temp.max25WARMED[CellsWin7]<-temp.max25W[CellsWin7]+1.8 
temp.max25WARMED[CellsWin8]<-temp.max25W[CellsWin8]+1.6 

########################################################################################################################
#Determine presence/absence based on temperature tolereance, erratic behaviour, preference
#Values from stickleback parameter summary file

#Create, tolerance, errat, pref

tol<-errat<-pref<-matrix(NA,nrow=dim(temp.max)[1],ncol=dim(temp.max)[2])
m<-1

for (i in 1:dim(temp.max)[1])
	{
		n<-1
		for (p in 1:dim(temp.max25)[2])
		{
		if (is.na(temp.min25[i,p])=='TRUE') 
			{
			n<-n+1
			next		
			}
		if (is.na(temp.min25[i,p])=='TRUE') 
			{
			n<-n+1
			next		
			}
		if (is.na(temp.min[i,p])=='TRUE') next
		if (temp.min25[i,p]<tol.lb | temp.max25[i,p]>tol.ub) tol[m,n]<-0
		if (temp.min25[i,p]>=tol.lb & temp.max25[i,p]<=tol.ub) tol[m,n]<-1
		if (temp.min25[i,p]<errat.lb | temp.max25[i,p]>errat.ub) errat[m,n]<-0
		if (temp.min25[i,p]>=errat.lb & temp.max25[i,p]<=errat.ub) errat[m,n]<-1 
		if (temp.med[i,p]<pref.lb | temp.med[i,p]>pref.ub) pref[m,n]<-0
		if (temp.med[i,p]>=pref.lb & temp.med[i,p]<=pref.ub) pref[m,n]<-1	
		n<-n+1
		}
		m<-m+1
		}

combo.t.e<-matrix(NA,nrow=dim(temp.max)[1],ncol=dim(temp.max)[2])

#Create combined raster of tol and erratic behav
for (i in 1:dim(tol)[1])
	{
		n<-1
		for (p in 1:dim(tol)[2])
		{
		if (is.na(tol[i,p])=='TRUE' & is.na(errat[i,p])=='TRUE') 
		{
			combo.t.e[i,p]<-NA 	
			next
			}
		if (tol[i,p]==0 & errat[i,p]==0) combo.t.e[i,p]<-0
		if (tol[i,p]==0 & errat[i,p]==1) combo.t.e[i,p]<-1
		if (tol[i,p]==1 & errat[i,p]==0) combo.t.e[i,p]<-2
		if (tol[i,p]==1 & errat[i,p]==1) combo.t.e[i,p]<-3
			}
		}

tolW<-erratW<-matrix(NA,nrow=dim(temp.max25WARMED)[1],ncol=dim(temp.max25WARMED)[2])
m<-1

########################################################################################################################

#For the warmed world

#Create, tolerance, errat, pref
for (i in 1:dim(temp.min25WARMED)[1])
	{
		n<-1
		for (p in 1:dim(temp.min25WARMED)[2])
		{
		if (is.na(temp.min25WARMED[i,p])=='TRUE') 
			{
			n<-n+1
			next		
			}
		if (is.na(temp.min25WARMED[i,p])=='TRUE') 
			{
			n<-n+1
			next		
			}
		if (is.na(temp.min25WARMED[i,p])=='TRUE') next
		if (temp.min25WARMED[i,p]<evol.tol.lb | temp.max25WARMED[i,p]>evol.tol.ub) tolW[m,n]<-0
		if (temp.min25WARMED[i,p]>=evol.tol.lb & temp.max25WARMED[i,p]<=evol.tol.ub) tolW[m,n]<-1
		if (temp.min25WARMED[i,p]<evol.errat.lb | temp.max25WARMED[i,p]>evol.errat.ub) erratW[m,n]<-0
		if (temp.min25WARMED[i,p]>=evol.errat.lb & temp.max25WARMED[i,p]<=evol.errat.ub) erratW[m,n]<-1
		n<-n+1
		}
		m<-m+1
		}

combo.t.eW<-matrix(NA,nrow=dim(temp.max25WARMED)[1],ncol=dim(temp.max25WARMED)[2])

#Create combined raster of tol and erratic behav
for (i in 1:dim(tolW)[1])
	{
		n<-1
		for (p in 1:dim(tolW)[2])
		{
		if (is.na(tolW[i,p])=='TRUE' & is.na(erratW[i,p])=='TRUE') 
		{
			combo.t.eW[i,p]<-NA 	
			next
			}
		if (tolW[i,p]==0 & erratW[i,p]==0) combo.t.eW[i,p]<-0
		if (tolW[i,p]==0 & erratW[i,p]==1) combo.t.eW[i,p]<-1
		if (tolW[i,p]==1 & erratW[i,p]==0) combo.t.eW[i,p]<-2
		if (tolW[i,p]==1 & erratW[i,p]==1) combo.t.eW[i,p]<-3
			}
		}

#Change outputs into rasters
TOL<-raster(tol)
ERRAT<-raster(errat)
PREF<-raster(pref)
COMBOTE<-raster(combo.t.e)
TOLW<-raster(tolW)
ERRATW<-raster(erratW)
COMBOTEW<-raster(combo.t.eW)

#Set projection
proj4string(TOL)<-CRS("+init=epsg:4326")
proj4string(ERRAT)<-CRS("+init=epsg:4326")
proj4string(PREF)<-CRS("+init=epsg:4326")
proj4string(COMBOTE)<-CRS("+init=epsg:4326")
proj4string(TOLW)<-CRS("+init=epsg:4326")
proj4string(ERRATW)<-CRS("+init=epsg:4326")
proj4string(COMBOTEW)<-CRS("+init=epsg:4326")

#Set extent
extent(TOL)<-extent(ERRAT)<-extent(PREF)<-extent(COMBOTE)<-extent(temp.min)
extent(TOLW)<-extent(ERRATW)<-extent(COMBOTEW)<-extent(temp.min25WARMED)

#Save rasters (use ",overwrite=TRUE" to overwrite files, not in code to avoid accidents)
writeRaster(TOL,"Processed_files/TOL_Nov30_rcp4_adjPVE_minOnly.asc",format="ascii")
writeRaster(ERRAT,"Processed_files/ERRAT_Nov30_rcp4_adjPVE_minOnly.asc",format="ascii")
writeRaster(PREF,"Processed_files/PREF_Nov30_rcp4_adjPVE_minOnly.asc",format="ascii")
writeRaster(COMBOTE,"Processed_files/COMBOTE_Nov30_rcp4_adjPVE_minOnly.asc",format="ascii")
writeRaster(TOLW,"Processed_files/TOLW_Nov30_rcp4_adjPVE_minOnly.asc",format="ascii")
writeRaster(ERRATW,"Processed_files/ERRATW_Nov30_rcp4_adjPVE_minOnly.asc",format="ascii")
writeRaster(COMBOTEW,"Processed_files/COMBO_TOL_E_W_Nov30_rcp4_adjPVE_minOnly.asc",format="ascii")

