library(spdep)
library(rgdal)
library(maptools)
library(sp)
library(RColorBrewer)
library(classInt)
library(GISTools)
library(maps)

# Set working directory
setwd("...")

# reading maps with rgdal::
cal<-readOGR(".", "CA_State_TIGER2016") 
county<-readOGR(".", "CA_Counties_TIGER2016")

# changing projections
cal<-spTransform(cal, CRS("+proj=longlat +datum=NAD83"))
county<-spTransform(county, CRS("+proj=longlat +datum=NAD83"))

plot(cal)
plot(county)


# reading point data as data.frame
houses<-read.csv("California_Houses.csv", header=TRUE, dec=".", sep=",")
summary(houses)

houses.sp <- houses
# class change from data.frame to sp
coordinates(houses.sp)<-c("Longitude","Latitude") # conversion to sp

# adding a projection to point data
proj4string(houses.sp)<-CRS("+proj=longlat +datum=NAD83") 
houses.sp<-spTransform(houses.sp, CRS("+proj=longlat +datum=NAD83")) 


# conversion from sp to sf
# houses<-st_as_sf(houses)

plot(county)
points(houses$Longitude, houses$Latitude, pch=".")

