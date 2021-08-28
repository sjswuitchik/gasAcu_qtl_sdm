#Let's make a salinity map! This map will have:

#Only the bathymetry where stickleback can live (see model parameterization)
#Takes into account Sea ice (but you can change that by using a diff Bathy file)

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
head(suit@data)

#Load salinity
Sal.csv<-read.csv("~/Desktop/UCalgary_PhD/Main_publication/Coding/req_files/01 Jan.csv")
names(Sal.csv)

#Create raster for one depth
sal <- rasterFromXYZ(Sal.csv[, c(2,1,3)]) #If want a diff depth change the final column to some other value (col 2,1 give lon and lat)

#Mask salinity based on suitable habitat (bathy and sea ice)
extent(sal)<-extent(suit)
sal.rast<-mask(sal,suit) 

#Create a plot
plot(sal.rast)
