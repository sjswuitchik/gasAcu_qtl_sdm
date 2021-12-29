#Let's make a salinity map! This map will have:

#Only the bathymetry where stickleback can live (see model parameterization)
#Takes into account sea ice (but you can change that by using a diff Bathy file)

rm(list=ls(all=TRUE))
library (ncdf4)
library(raster)
library(rgeos)
library(sp)
library(rgdal)
library(maptools)


###########################################################################################################################
#Load shapefile of bathymetry/sea ice availability 
suit<-readOGR("~/Desktop/UCalgary_PhD/Main_publication/Coding/req_files/Bathy&SeaIce-prefered")

#Check data loaded correctly
#head(suit@data)

#Load bathymetry to get coastlines (nc file)
Bathy<-nc_open("~/Desktop/UCalgary_PhD/Main_publication/Coding/req_files/GEBCO_2014_2D_-179.7777_45.2207_-120.8539_76.3978.nc")

#Find lat and lon of cells
lon<-ncvar_get(Bathy,"lon")
lat<-ncvar_get(Bathy,"lat")

#Load salinity
Sal.csv<-read.csv("~/Desktop/UCalgary_PhD/Main_publication/Coding/req_files/01 Jan.csv")
names(Sal.csv)

#Create raster for one depth
sal <- rasterFromXYZ(Sal.csv[, c(2,1,3)]) #If want a diff depth change the final column to some other value (col 2,1 give lon and lat)

#Mask salinity based on suitable habitat (bathy and sea ice)
extent(sal)<-extent(suit)
sal.rast<-mask(sal,suit) 

#Get depth/elevation data
Depths<-ncvar_get(Bathy,"elevation")
#Transpose Depth matrix so map will be right way up
Depths <- t(Depths)
#Create raster
D<-raster(Depths,xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))
#Flip raster so map will be right way up
D<-flip(D,direction='y')
#Alter raster for plotting
D[D<0]<-NA #Remove all cells below sea level (i.e. land remains)
D[D>=0]<-1 #Give all cells above sea level the same value

#Create a plot
plot(D, col="grey",axes=F,legend=F)
plot(sal.rast, add = T)
