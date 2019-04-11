# Characterising geographic data and save to locations.Rdata for future use
# Makes a dataframe that includes site specific colors and additional set of geographic variables are 
# created with an axis that captures distance to shoreline and an orthogonal eigenvector (parallel to shore line, roughly)
# Written by Cynthia Riginos 2017-2019

library(plotrix)
library(maps)
library(rgeos)
library(maptools)
library(raster)
library(rgdal)
library(vegan)


locations<-read.csv("ORIGDATA/ReefID_centroids_tenuis.csv", header =TRUE)  #these are lat and long for centers of each reef location
#locations<-read.csv("ORIGDATA/ReefID_centroids_millepora.csv", header =TRUE)  #these are lat and long for centers of each reef location

#ASSIGN COLORS BY LATITUDE
loc_color<-color.scale(1:100, c(0,1,1), c(0,1,0),c(1,1,0),color.spec="rgb")  #blue to red with white in middle
latvalue<-(locations$Lat+24)/16
latvalue2<-(100*round(latvalue, digits = 2))
locations$loc_color<-loc_color[latvalue2]

#quick check that colors and locations look ok
plot(locations$Long, locations$Lat, ylab = "Latitude", xlab = "Longitude", type ="n") 
map(database = "world2", add= TRUE, interior=FALSE)  #world2 is Pacific centered
points(locations$Long, locations$Lat, col = locations$loc_color, pch=19, cex = 1.5) 
text(locations$Long, locations$Lat, locations$ReefID, cex=0.6, pos=3)

##GET DISTANCE TO COASTLINE
#Coordinate systems
wgs.84 <- CRS("+init=epsg:4326")
#from Geoscience Australia - http://spatialreference.org/ref/epsg/3112/
epsg.3112<-"+proj=lcc +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
coast <- readShapeLines("ORIGDATA/ne_50m_coastline/ne_50m_coastline.shp")  #http://www.naturalearthdata.com
proj4string(coast)<-wgs.84

#Create locations spatial points dataframe
locations.sp<-locations[,c(3,2)]
coordinates(locations.sp)<- c("Long", "Lat")
proj4string(locations.sp)<-wgs.84

#Transform
locations.proj.sp   <- spTransform(locations.sp,CRS(epsg.3112))
coast.proj <- spTransform(coast,CRS(epsg.3112))

#Extract distance to coastline
for (p in 1:length(locations.proj.sp)) {
  locations[p,"DistToCoastline"]<-round(gDistance(locations.proj.sp[p],coast.proj)/1000, digits=0)
}

cor(locations$Lat,locations$DistToCoastline) #AT: -0.66; AM: -0.41
cor(locations$Long,locations$DistToCoastline) #AT: 0.83; AM: 0.61

#Use RDA to find eigenvector and eigenvalues perpendicular to shoreline distance eigenvector
long.lat<-locations[,c(3,2)]
geog.pca<-rda(long.lat~DistToCoastline, data=locations)
summary(geog.pca) #AT: RDA1=54.95% variance, PC1 = 44.00% (total 98.9%): AM: RDA1=26%, PC1 = 73 (total: 98.9)
newscores<-as.data.frame(scores(geog.pca, scaling=1, choices=1:2)$sites)

plot(newscores$RDA1, -newscores$PC1, pch=20, col=locations$loc_color) #quick check

locations$DistToCoastPC<-newscores$RDA1
locations$ParallelToCoastPC<- -newscores$PC1

save(locations, file="CLEANDATA/tenuis/locations.aten.RData" )
#save(locations, file="CLEANDATA/millepora/locations.mill.RData" )


rm(list=ls())

