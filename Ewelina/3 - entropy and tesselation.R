library(spdep)
library(rgdal)
library(maptools)
library(sp)
library(RColorBrewer)
library(classInt)
library(GISTools)
library(maps)

############################## Tesselated houses locations ##############################

# Read California shapefile
setwd('/Users/ewelinaplachimowicz/Documents/GitHub/housing-prices-spatial-econometrics/model/')

county<-readOGR(".", "CA_Counties_TIGER2016")
county<-spTransform(county, CRS("+proj=merc +datum=NAD83"))
plot(county)

# Read data
setwd('/Users/ewelinaplachimowicz/Documents/GitHub/housing-prices-spatial-econometrics/data/')

houses<-read.csv("California_Houses.csv", header=TRUE, dec=".", sep=",")
summary(houses)

houses.sp <- houses
# class change from data.frame to sp
coordinates(houses.sp)<-c("Longitude","Latitude") # conversion to sp

# adding a projection to point data
proj4string(houses.sp)<-CRS("+proj=longlat +datum=NAD83") 
houses.sp<-spTransform(houses.sp, CRS("+proj=longlat +datum=NAD83")) 

# Tesselated houses locations in California
houses.sp1<-spTransform(houses.sp, CRS("+proj=merc +datum=NAD83")) 

county.owin<-as(county, "owin")
county.ppp<-ppp(x=houses.sp1@coords[,1], y=houses.sp1@coords[,2], window=county.owin)
county.tes<-dirichlet(county.ppp) # Dirichlet tessellation

a_cal<-tile.areas(county.tes)
a1_cal<-a_cal/sum(a_cal)
ent1_cal<-sum(-1*a1_cal*log(a1_cal))
ent1_cal
n_cal<-length(a_cal)
ent.ref_cal<-log(1/n_cal)*(-1)
ent.ref_cal
ent.rel_cal<-ent1_cal/ent.ref_cal
ent.rel_cal

plot(county.tes, main="Tesselated houses locations in California")
plot(county.ppp, add=TRUE, pch=".", col="darkblue", cex=2.8)
text(-13857275, 3809956, paste("Shannon entropy=", round(ent1_cal,2)), cex=0.8)
text(-13857275, 3889956, paste("Relative H entropy=", round(ent.rel_cal,2)), cex=0.8)
text(-13857275, 3969956, paste("Number of points =", round(n_cal,2)), cex=0.8)

# Random locations for comparison
countyR.sp<-spsample(county, n=1000, type="random")

plot(county)
points(countyR.sp, pch=".", cex=2, col="darkblue" )

countyR.ppp<-ppp(x=countyR.sp@coords[,1], y=countyR.sp@coords[,2], window=county.owin)
countyR.tes<-dirichlet(countyR.ppp) # Dirichlet tessellation

a_cal_ran<-tile.areas(countyR.tes)
a1_cal_ran<-a_cal_ran/sum(a_cal_ran)
ent1_cal_ran<-sum(-1*a1_cal_ran*log(a1_cal_ran))
n_cal_ran<-length(a_cal_ran)
ent.ref_cal_ran<-log(1/n_cal_ran)*(-1)
ent.ref_cal_ran
ent.rel_cal_ran<-ent1_cal_ran/ent.ref_cal_ran
ent.rel_cal_ran

plot(countyR.tes, main="Tesselated houses locations in California - theoretical (random) distribution of points")
plot(countyR.ppp, add=TRUE, pch=".", col="darkblue", cex=2.8)
text(-13857275, 3809956, paste("Shannon entropy=", round(ent1_cal_ran,2)), cex=0.8)
text(-13857275, 3889956, paste("Relative H entropy=", round(ent.rel_cal_ran,2)), cex=0.8)
text(-13857275, 3969956, paste("Number of points =", round(n_cal_ran,2)), cex=0.8)


####################### Analysis for 2 biggest cities in California - Los Angeles and San Diego ####################### 

# Crate file for each city
la <- county[county@data$NAME=="Los Angeles",]
sd <- county[county@data$NAME=="San Diego",]

# Plot cities
plot(la)
points(houses.sp1, pch=".")

plot(sd)
points(houses.sp1, pch=".")

# Los Angeles
la.owin<-as(la, "owin")
la.ppp<-ppp(x=houses.sp1@coords[,1], y=houses.sp1@coords[,2], window=la.owin)
la.tes<-dirichlet(la.ppp) # Dirichlet tessellation

a_la<-tile.areas(la.tes)
a1_la<-a_la/sum(a_la)
ent1_la<-sum(-1*a1_la*log(a1_la))
n_la<-length(a_la)
ent.ref_la<-log(1/n_la)*(-1)
ent.ref_la
ent.rel_la<-ent1_la/ent.ref_la # Relative H
ent.rel_la

plot(la.tes, main="Tesselated houses locations in Los Angeles")
plot(la.ppp, add=TRUE, pch=".", col="darkblue", cex=2.8)
text(-13271646, 3839060, paste("Shannon entropy=", round(ent1_la,2)), cex=0.8)
text(-13271646, 3859060, paste("Relative H entropy=", round(ent.rel_la,2)), cex=0.8)
text(-13271646, 3879060, paste("Number of points =", round(n_la,2)), cex=0.8)

# San Diego
sd.owin<-as(sd, "owin")
sd.ppp<-ppp(x=houses.sp1@coords[,1], y=houses.sp1@coords[,2], window=sd.owin)
sd.tes<-dirichlet(sd.ppp) # Dirichlet tessellation

a_sd<-tile.areas(sd.tes)
a1_sd<-a_sd/sum(a_sd)
ent1_sd<-sum(-1*a1_sd*log(a1_sd))
n_sd<-length(a_sd)
ent.ref_sd<-log(1/n_sd)*(-1)
ent.ref_sd
ent.rel_sd<-ent1_sd/ent.ref_sd # Relative H
ent.rel_sd

plot(sd.tes, main="Tesselated houses locations in San Diego")
plot(sd.ppp, add=TRUE, pch=".", col="darkblue", cex=2.8)
text(-13092407, 3809956, paste("Shannon entropy=", round(ent1_sd,2)), cex=0.8)
text(-13092407, 3819956, paste("Relative H entropy=", round(ent.rel_sd,2)), cex=0.8)
text(-13092407, 3829956, paste("Number of points =", round(n_sd,2)), cex=0.8)