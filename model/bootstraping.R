# loading packages
library(spdep)
library(rgdal)
library(maptools)
library(sp)
library(RColorBrewer)
library(doBy)
library(cluster)
library(e1071)
library(spatstat)
library(sampling)
library(spatialreg)
library(stargazer)


setwd("...")
getwd()

# reading maps with rgdal::
cal<-readOGR(".", "CA_State_TIGER2016") 
county<-readOGR(".", "CA_Counties_TIGER2016")


# changing projections
cal<-spTransform(cal, CRS("+proj=longlat +datum=NAD83"))
county<-spTransform(county, CRS("+proj=longlat +datum=NAD83"))


# reading point data as data.frame
houses<-read.csv("California_Houses.csv", header=TRUE, dec=".", sep=",")


library(dplyr)
# group data to polygons names
houses <- houses %>% group_by(houses$Latitude, houses$Longitude) %>% summarise(Median_House_Value = median(Median_House_Value),
                                                                               Median_Income = median(Median_Income),
                                                                               Median_Age = median(Median_Age),
                                                                               Tot_Rooms = median(Tot_Rooms),
                                                                               Tot_Bedrooms = median(Tot_Bedrooms),
                                                                               Population = median(Population),
                                                                               Households = median(Households),
                                                                               Distance_to_coast = median(Distance_to_coast),
                                                                               Distance_to_LA = median(Distance_to_LA),
                                                                               Distance_to_SanDiego = median(Distance_to_SanDiego),
                                                                               Distance_to_SanJose = median(Distance_to_SanJose),
                                                                               Distance_to_SanFrancisco = median(Distance_to_SanFrancisco),
                                                                               Latitude = min(Latitude),
                                                                               Longitude = min(Longitude))

# drop NAs
houses <- houses %>% drop_na()

houses <- as.data.frame(houses)


summary(houses)
dim(houses)

# random sorting of dataset
houses$los<-runif(dim(houses)[1], 0,1)
houses<-orderBy(~los, data=houses)

# parameters of simulation 
n.col<-50 # number of iterations
n.row<-1000 # number of obs in a sample
nnew<-20 # number of new points in forecast


############################################
# full procedure – bootstrapped models
###################	#################

# split of data and correction of location by epsilon
# in in-sample and out-of-sample datasets
epsilon.x<-rnorm(dim(houses)[1], mean=0, sd=0.015)
epsilon.y<-rnorm(dim(houses)[1], mean=0, sd=0.015)
houses$xxe<-houses$Longitude + epsilon.x
houses$yye<-houses$Latitude + epsilon.y
houses.in<-houses[1:12000,]
houses.out<-houses[12001:12590,]
houses.out$los<-runif(dim(houses.out)[1], 0,1)
houses.out<-orderBy(~los, data=houses.out)

# selector matrix – with k-means clustering
# assigning points to clusters with k-means – to get irregular shapes for sampling
#remove(groups)
groups<-kmeans(houses.in[,c("Longitude", "Latitude")], n.row/100)
houses.in$kmean<-groups$cluster

# selector selects the rows of observations for every iteration and saves in matrix
# later it will be used to recover on which data the best model was estimated
# strata() from sampling:: package samples observations from groups (clusters)
selector<-matrix(0, nrow=n.row, ncol=n.col)
for(i in 1:n.col){
  vec<-sample(1:dim(houses.in)[1], n.row, replace=FALSE)
  x<-strata(houses.in, "kmean", size=rep(100, times=n.row/100), method="srswor") # from sampling::
  selector[,i]<-x$ID_unit}

#####################################
# all models – objects to store results of iterations
coef.ols<-matrix(0, nrow=n.col, ncol=12)
coef.sem<-matrix(0, nrow=n.col, ncol=12)
coef.sar<-matrix(0, nrow=n.col, ncol=12)
coef.sdm<-matrix(0, nrow=n.col, ncol=23)

error.ols<-matrix(0, nrow=n.col, ncol=12)
error.sem<-matrix(0, nrow=n.col, ncol=12)
error.sar<-matrix(0, nrow=n.col, ncol=12)
error.sdm<-matrix(0, nrow=n.col, ncol=23)
sig.ols<-matrix(0, nrow=n.col, ncol=12)

fitted.ols<-matrix(0, nrow=n.row, ncol=n.col)
fitted.sem<-matrix(0, nrow=n.row, ncol=n.col)
fitted.sar<-matrix(0, nrow=n.row, ncol=n.col)
fitted.sdm<-matrix(0, nrow=n.row, ncol=n.col)

quality<-matrix(0, nrow=n.row, ncol=9) # R2, AIC.ols, AIC.sem, AIC.sar, AIC.sdm, BIC.ols, BIC.sem, BIC.sar, BIC.sdm
time<-matrix(0, nrow=n.row, ncol=4) # time.ols, time.sem, time.sar, time.sdm
price<-matrix(0, nrow=n.row, ncol=n.col)
spatial<-matrix(0, nrow=n.col, ncol=3) # lambda.sem, rho.sar, rho.sdm

eq<-Median_House_Value ~ Median_Income + Median_Age +  Tot_Rooms + Tot_Bedrooms + Population + Households + Distance_to_coast + Distance_to_LA + Distance_to_SanDiego + Distance_to_SanJose + Distance_to_SanFrancisco


# big loop – estimation of all models (OLS, SEM, SAR, SDM) on the same subsets with time measurement
for(i in 1:n.col){
  houses_sample<-houses.in[selector[,i],]
  price[,i]<-houses_sample$Median_House_Value # saving y values for each iteration in price matrix
  
  crds<-as.matrix(houses.in[selector[,i],c("Longitude", "Latitude")])  # spatial weights matrix W 
  pkt.knn<-knearneigh(crds, k=5, longlat = NULL) # object in knn class k=5 planar coordinates
  pkt.k.nb<-knn2nb(pkt.knn) 
  pkt.k.sym.nb<-make.sym.nb(pkt.k.nb)
  pkt.k.sym.listw<-nb2listw(pkt.k.sym.nb)
  
  # model OLS
  start.time <- Sys.time()
  model.ols<-lm(eq, data=houses_sample)
  end.time <- Sys.time()
  time.ols<- difftime(end.time, start.time, units="secs")
  # model SEM
  start.time <- Sys.time()
  model.sem<-errorsarlm(eq, data=houses_sample, pkt.k.sym.listw, method="LU")
  end.time <- Sys.time()
  time.sem<- difftime(end.time, start.time, units="secs")
  # model SAR
  start.time <- Sys.time()
  model.sar<-lagsarlm(eq, data=houses_sample, pkt.k.sym.listw, method="LU")
  end.time <- Sys.time()
  time.sar<- difftime(end.time, start.time, units="secs")
  # model SDM
  start.time <- Sys.time()
  model.sdm<-lagsarlm(eq, data=houses_sample, pkt.k.sym.listw, method="LU", type="mixed")
  end.time <- Sys.time()
  time.sdm<- difftime(end.time, start.time, units="secs")
  
  # saving the results
  coef.ols[i,]<-summary(model.ols)$coefficients[,1]
  error.ols[i,]<-summary(model.ols)$coefficients[,2]
  sig.ols[i,]<-summary(model.ols)$coefficients[,4]
  fitted.ols[,i]<-model.ols$fitted.values
  
  coef.sem[i,]<-model.sem$coefficients
  error.sem[i,]<-model.sem$rest.se
  fitted.sem[,i]<-model.sem$fitted.values
  
  coef.sar[i,]<-model.sar$coefficients
  error.sar[i,]<-model.sar$rest.se
  fitted.sar[,i]<-model.sar$fitted.values
  
  coef.sdm[i,]<-model.sdm$coefficients
  error.sdm[i,]<-model.sdm$rest.se
  fitted.sdm[,i]<-model.sdm$fitted.values
  
  quality[i,1]<-summary(model.ols)$r.squared
  quality[i,2]<-AIC(model.ols) # AIC.ols
  quality[i,3]<-AIC(model.sem) # AIC.sem
  quality[i,4]<-AIC(model.sar) # AIC.sar
  quality[i,5]<-AIC(model.sdm) # AIC.sdm
  quality[i,6]<-BIC(model.ols) # BIC.ols
  quality[i,7]<-BIC(model.sem) # BIC.sem
  quality[i,8]<-BIC(model.sar) # BIC.sar
  quality[i,9]<-BIC(model.sdm) # BIC.sdm
  
  time[i,1]<-time.ols
  time[i,2]<-time.sem
  time[i,3]<-time.sar
  time[i,4]<-time.sdm
  spatial[i,1]<-model.sem$lambda
  spatial[i,2]<-model.sar$rho
  spatial[i,3]<-model.sdm$rho}


############################
# selection of the best model with PAM
############################

## PAM for OLS
c1.ols<-pam(coef.ols,1) #cluster::pam(), works for n<65536
summary(c1.ols)
plot(c1.ols)

c1.ols$clustering #  clustering vector – only 1 values
c1.ols$medoids # coefficients of selected best model
c1.ols$id.med # which iteration (model) is most representative
#hopkins(coef.ols, n=nrow(coef.ols)-1) # takes long…

#############################
## PAM for SEM
c1.sem<-pam(cbind(coef.sem, spatial[,1]),1) #cluster::pam(), works for n<65536
summary(c1.sem)
#plot(c1.sem)

c1.sem$clustering # clustering vector – only 1 values
c1.sem$medoids # # coefficients of selected best model 
c1.sem$id.med # which iteration (model) is most representative
#hopkins(coef.sem, n=nrow(coef.sem)-1) # takes long…

#############################
## PAM for SAR
c1.sar<-pam(cbind(coef.sar, spatial[,2]),1) #cluster::pam(), works for n<65536
summary(c1.sar)
#plot(c1.sar)

c1.sar$clustering # 
c1.sar$medoids # 
c1.sar$id.med # 
#hopkins(coef.sar, n=nrow(coef.sar)-1) # 

#############################
## PAM for SDM
c1.sdm<-pam(cbind(coef.sdm, spatial[,3]),1) #cluster::pam(), works for n<65536
summary(c1.sdm)
#plot(c1.sdm)

c1.sdm$clustering # 
c1.sdm$medoids # 
c1.sdm$id.med # 
#hopkins(coef.sem, n=nrow(coef.sem)-1) # 

#############################
# full SEM, SAR & SDM models
eq<-Median_House_Value ~ Median_Income + Median_Age +  Tot_Rooms + Tot_Bedrooms + Population + Households + Distance_to_coast + Distance_to_LA + Distance_to_SanDiego + Distance_to_SanJose + Distance_to_SanFrancisco

houses<-houses.in
crds<-as.matrix(houses.in[ ,c("Longitude", "Latitude")])
pkt.knn<-knearneigh(crds, k=5, longlat = NULL) # knn k=5 planar coordinates
pkt.k.nb<-knn2nb(pkt.knn) 
pkt.k.sym.nb<-make.sym.nb(pkt.k.nb)
pkt.k.sym.listw<-nb2listw(pkt.k.sym.nb)
model.ols.full<-lm(eq, data=houses)
model.sem.full<-errorsarlm(eq, data=houses, pkt.k.sym.listw, method="LU")
model.sar.full<-lagsarlm(eq, data=houses, pkt.k.sym.listw, method="LU")
model.sdm.full<-lagsarlm(eq, data=houses, pkt.k.sym.listw, method="LU", type="mixed")

# install.packages("stargazer")

#########################################
# validation for out-of-sample
#########################################

# Randal maps with rgdal::
cal<-readOGR(".", "CA_State_TIGER2016") 
county<-readOGR(".", "CA_Counties_TIGER2016")


# changing projections
cal<-spTransform(cal, CRS("+proj=longlat +datum=NAD83"))
cal<-spTransform(cal, CRS("+proj=merc +datum=NAD83"))

# preparing objects for spatstat
# library(rgdal)
region.owin<-as.owin(cal) # rgdal:: requires planar coordinates


# OLS
# tessellation of observations from the selected best model
points<-data.frame(x=houses.in[selector[,c1.ols$id.med],c("Longitude")], y=houses.in[selector[,c1.ols$id.med],c("Latitude")]) 
points.sp<-SpatialPoints(points) # new points in sp class - spherical 
proj4string(points.sp)<-CRS("+proj=longlat +datum=NAD83") # spherical 
points.sp<-spTransform(points.sp, CRS("+proj=merc +datum=NAD83")) # planar 
region.ppp<-ppp(x=points.sp@coords[,1], y=points.sp@coords[,2], window=region.owin) # points of ppp class 
region.tes<-dirichlet(region.ppp) # Dirichlet tesselation 

tes.poly<-as(region.tes, "SpatialPolygons") 
proj4string(tes.poly)<-CRS("+proj=merc +datum=NAD83") 
tes.poly<-spTransform(tes.poly, CRS("+proj=merc +datum=NAD83")) #planar 

plot(region.tes, main=" ") # tessellation plot
plot(region.ppp, add=TRUE, pch=".", col="darkblue", cex=2)

nnew<-100 # number of new points in the forecast 
forecasts1<-matrix(0, nrow=nnew, ncol=5) 
colnames(forecasts1)<-c("predicted y","real y","crds x","crds y", "diff")

# restoring the medoid model
dane.x<-houses.in[selector[,c1.ols$id.med],] # we choose the data which were used in estimation of the best model
crds<-as.matrix(dane.x[, c("Longitude", "Latitude")]) # we create W
RAMSE.med.ols<-(sum((houses.in[selector[,c1.ols$id.med], c('Median_House_Value')]-fitted.ols[,c1.ols$id.med])^2)/n.row)^(0.5) 
RAMSE.med.ols

pkt.knn<-knearneigh(crds, k=5) # knn object # to use in Moran test in case of OLS
pkt.k.nb<-knn2nb(pkt.knn) 
pkt.k.sym.nb<-make.sym.nb(pkt.k.nb)
pkt.k.sym.listw<-nb2listw(pkt.k.sym.nb)

eq<-Median_House_Value ~ Median_Income + Median_Age +  Tot_Rooms + Tot_Bedrooms + Population + Households + Distance_to_coast + Distance_to_LA + Distance_to_SanDiego + Distance_to_SanJose + Distance_to_SanFrancisco


model.ols<-lm(eq, data=dane.x)
moran.test(model.ols$residuals, pkt.k.sym.listw)

# selection of out-of-sample points for predictions
points.pred<-SpatialPoints(houses.out[1:nnew, c("Longitude", "Latitude")]) # selection of a new point
proj4string(points.pred)<-CRS("+proj=longlat +datum=NAD83") # spherical projection
points.pred<-spTransform(points.pred, CRS("+proj=merc +datum=NAD83")) 
a1<-over(points.pred, tes.poly) # assigning points to tessellation tiles

# completing the draw when NA occurs # determining the number of new points to be drawn (as from-to) 
a2<-nnew+1 # from … 
a3<-which(is.na(a1)) 
a4<-a2+length(a3)-1 # to … 
points.pred2<-SpatialPoints(houses.out[a2:a4, c("Longitude", "Latitude")]) # new points 
proj4string(points.pred2)<-CRS("+proj=longlat +datum=NAD83") # spherical 
points.pred2<-spTransform(points.pred2, CRS("+proj=merc +datum=NAD83")) 
points.pred2 
a5<-over(points.pred2, tes.poly) # putting new points on the tile 
a5 
a1[which(is.na(a1))]<-a5 # overwriting with new points

# loop for forecasts for new points # there is a separate match for each point 
for(i in 1:nnew){ 
  # point by point - assigning new data to the old data set 
  dane.x.new<-dane.x 
  xxx<-houses.out[i, ] 
  dane.x.new[a1[i],]<-xxx 
  rownames(dane.x.new)<-1:dim(dane.x.new)[1] 
  
  # prediction for out-of-sample calibrated SDM model
  pred<-predict(model.ols, newdata=dane.x.new) 
  
  pred[a1[i]] # prediction for a new point 
  xxx[,c('Median_House_Value')] # empirical y value of the new point 
  forecasts1[i,1]<- pred[a1[i]] # predicted value of y 
  forecasts1[i,2]<- xxx[, c('Median_House_Value')] # empirical value of y 
  forecasts1[i,3]<-xxx[, c("Latitude")] # x coordinates 
  forecasts1[i,4]<-xxx[, c("Longitude")] # y coordinates 
} 
forecasts1[,5]<-(forecasts1[,1]-forecasts1[,2])^2 
RAMSE.ols<-(mean(forecasts1[,5]))^0.5 
head(forecasts1)
RAMSE.ols

##########################
# SEM
# tessellation of observations from the selected best model
points<-data.frame(x=houses.in[selector[,c1.sem$id.med],c("Longitude")], y=houses.in[selector[,c1.sem$id.med],15]) 
points.sp<-SpatialPoints(points) # new points in sp class - spherical 
proj4string(points.sp)<-CRS("+proj=longlat +datum=NAD83") # spherical 
points.sp<-spTransform(points.sp, CRS("+proj=merc +datum=NAD83")) # planar 
region.ppp<-ppp(x=points.sp@coords[,1], y=points.sp@coords[,2], window=region.owin) # points of ppp class 
region.tes<-dirichlet(region.ppp) # Dirichlet tesselation 

tes.poly<-as(region.tes, "SpatialPolygons") 
proj4string(tes.poly)<-CRS("+proj=merc +datum=NAD83") 
tes.poly<-spTransform(tes.poly, CRS("+proj=merc +datum=NAD83")) #planar 

plot(region.tes, main=" ") # tessellation plot
plot(region.ppp, add=TRUE, pch=".", col="darkblue", cex=2)

nnew<-100 # number of new points in the forecast 
forecasts1<-matrix(0, nrow=nnew, ncol=5) 
colnames(forecasts1)<-c("predicted y","real y","crds x","crds y", "diff")

# restoring the medoid model
dane.x<-houses.in[selector[,c1.sem$id.med],] # we choose the data which were used in estimation of the best model
crds<-as.matrix(dane.x[,c("Longitude", "Latitude")]) # we create W
RAMSE.med.sem<-(sum((houses.in[selector[,c1.sem$id.med], c("Median_House_Value")]-fitted.sem[,c1.sem$id.med])^2)/n.row)^(0.5) 
RAMSE.med.sem

pkt.knn<-knearneigh(crds, k=5) 
pkt.k.nb<-knn2nb(pkt.knn) 
pkt.k.sym.nb<-make.sym.nb(pkt.k.nb)
pkt.k.sym.listw<-nb2listw(pkt.k.sym.nb)

eq<-Median_House_Value ~ Median_Income + Median_Age +  Tot_Rooms + Tot_Bedrooms + Population + Households + Distance_to_coast + Distance_to_LA + Distance_to_SanDiego + Distance_to_SanJose + Distance_to_SanFrancisco

model.sem<-errorsarlm(eq, data=dane.x, pkt.k.sym.listw, method="LU")
moran.test(model.sem$residuals, pkt.k.sym.listw)

# selection of out-of-sample points for predictions
points.pred<-SpatialPoints(houses.out[1:nnew, c("Longitude","Latitude")]) # selection of a new point
proj4string(points.pred)<-CRS("+proj=longlat +datum=NAD83") # spherical projection
points.pred<-spTransform(points.pred, CRS("+proj=merc +datum=NAD83")) 
a1<-over(points.pred, tes.poly) # assigning points to tessellation tiles

# completing the draw when NA occurs # determining the number of new points to be drawn (as from-to) 
a2<-nnew+1 # from … 
a3<-which(is.na(a1)) 
a4<-a2+length(a3)-1 # to … 
points.pred2<-SpatialPoints(houses.out[a2:a4, c("Longitude","Latitude")]) # new points 
proj4string(points.pred2)<-CRS("+proj=longlat +datum=NAD83") # spherical 
points.pred2<-spTransform(points.pred2, CRS("+proj=merc +datum=NAD83")) 
points.pred2 
a5<-over(points.pred2, tes.poly) # putting new points on the tile 
a5 
a1[which(is.na(a1))]<-a5 # overwriting with new points

# loop for forecasts for new points # there is a separate match for each point 
for(i in 1:nnew){ 
  # point by point - assigning new data to the old data set 
  dane.x.new<-dane.x 
  xxx<-houses.out[i, ] 
  dane.x.new[a1[i],]<-xxx 
  rownames(dane.x.new)<-1:dim(dane.x.new)[1] 
  
  # prediction for out-of-sample calibrated SDM model
  pred<-predict(model.sem, newdata=dane.x.new, listw=pkt.k.sym.listw, legacy.mixed=TRUE) 
  
  pred[a1[i]] # prediction for a new point 
  xxx[,c("Median_House_Value")] # empirical y value of the new point 
  forecasts1[i,1]<- pred[a1[i]] # predicted value of y 
  forecasts1[i,2]<- xxx[, c("Median_House_Value")] # empirical value of y 
  forecasts1[i,3]<-xxx[, c("Latitude")] # x coordinates 
  forecasts1[i,4]<-xxx[, c("Longitude")] # y coordinates 
} 
forecasts1[,5]<-(forecasts1[,1]-forecasts1[,2])^2 
RAMSE.sem<-(mean(forecasts1[,5]))^0.5 
head(forecasts1)
RAMSE.sem

#############################
# SAR
# tessellation of observations from the selected best model
points<-data.frame(x=houses.in[selector[,c1.sar$id.med],c("Longitude")], y=houses.in[selector[,c1.sar$id.med],c("Latitude")]) 
points.sp<-SpatialPoints(points) # new points in sp class - spherical 
proj4string(points.sp)<-CRS("+proj=longlat +datum=NAD83") # spherical 
points.sp<-spTransform(points.sp, CRS("+proj=merc +datum=NAD83")) # planar 
region.ppp<-ppp(x=points.sp@coords[,1], y=points.sp@coords[,2], window=region.owin) # points of ppp class 
region.tes<-dirichlet(region.ppp) # Dirichlet tesselation 

tes.poly<-as(region.tes, "SpatialPolygons") 
proj4string(tes.poly)<-CRS("+proj=merc +datum=NAD83") 
tes.poly<-spTransform(tes.poly, CRS("+proj=merc +datum=NAD83")) #planar 

plot(region.tes, main=" ") # tessellation plot
plot(region.ppp, add=TRUE, pch=".", col="darkblue", cex=2)

nnew<-100 # number of new points in the forecast 
forecasts1<-matrix(0, nrow=nnew, ncol=5) 
colnames(forecasts1)<-c("predicted y","real y","crds x","crds y", "diff")

# restoring the medoid model
dane.x<-houses.in[selector[,c1.sar$id.med],] # we choose the data which were used in estimation of the best model
crds<-as.matrix(dane.x[,c("Longitude", "Latitude")]) # we create W
RAMSE.med.sar<-(sum((houses.in[selector[,c1.sar$id.med], c("Median_House_Value")]-fitted.sar[,c1.sar$id.med])^2)/n.row)^(0.5) 
RAMSE.med.sar

pkt.knn<-knearneigh(crds, k=5) 
pkt.k.nb<-knn2nb(pkt.knn) 
pkt.k.sym.nb<-make.sym.nb(pkt.k.nb)
pkt.k.sym.listw<-nb2listw(pkt.k.sym.nb)

model.sar<-lagsarlm(eq, data=dane.x, pkt.k.sym.listw, method="LU")
moran.test(model.sar$residuals, pkt.k.sym.listw)

# selection of out-of-sample points for predictions
points.pred<-SpatialPoints(houses.out[1:nnew, c("Longitude", "Latitude")]) # selection of a new point
proj4string(points.pred)<-CRS("+proj=longlat +datum=NAD83") # spherical projection
points.pred<-spTransform(points.pred, CRS("+proj=merc +datum=NAD83")) 
a1<-over(points.pred, tes.poly) # assigning points to tessellation tiles

# completing the draw when NA occurs # determining the number of new points to be drawn (as from-to) 
a2<-nnew+1 # from … 
a3<-which(is.na(a1)) 
a4<-a2+length(a3)-1 # to … 
points.pred2<-SpatialPoints(houses.out[a2:a4, c("Longitude", "Latitude")]) # new points 
proj4string(points.pred2)<-CRS("+proj=longlat +datum=NAD83") # spherical 
points.pred2<-spTransform(points.pred2, CRS("+proj=merc +datum=NAD83")) 
points.pred2 
a5<-over(points.pred2, tes.poly) # putting new points on the tile 
a5 
a1[which(is.na(a1))]<-a5 # overwriting with new points

# loop for forecasts for new points # there is a separate match for each point 
for(i in 1:nnew){ 
  # point by point - assigning new data to the old data set 
  dane.x.new<-dane.x 
  xxx<-houses.out[i, ] 
  dane.x.new[a1[i],]<-xxx 
  rownames(dane.x.new)<-1:dim(dane.x.new)[1] 
  
  # prediction for out-of-sample calibrated SDM model
  pred<-predict(model.sar, newdata=dane.x.new, listw=pkt.k.sym.listw, legacy.mixed=TRUE) 
  
  pred[a1[i]] # prediction for a new point 
  xxx[,c("Median_House_Value")] # empirical y value of the new point 
  forecasts1[i,1]<- pred[a1[i]] # predicted value of y 
  forecasts1[i,2]<- xxx[, c("Median_House_Value")] # empirical value of y 
  forecasts1[i,3]<-xxx[, c("Longitude")] # x coordinates 
  forecasts1[i,4]<-xxx[, c("Latitude")] # y coordinates 
} 
forecasts1[,5]<-(forecasts1[,1]-forecasts1[,2])^2 
RAMSE.sar<-(mean(forecasts1[,5]))^0.5 
head(forecasts1)
RAMSE.sar

#########################################
# SDM
# tessellation of observations from the selected best model
points<-data.frame(x=houses.in[selector[,c1.sdm$id.med],c("Longitude")], y=houses.in[selector[,c1.sdm$id.med],c("Latitude")]) 
points.sp<-SpatialPoints(points) # new points in sp class - spherical 
proj4string(points.sp)<-CRS("+proj=longlat +datum=NAD83") # spherical 
points.sp<-spTransform(points.sp, CRS("+proj=merc +datum=NAD83")) # planar 
region.ppp<-ppp(x=points.sp@coords[,1], y=points.sp@coords[,2], window=region.owin) # points of ppp class 
region.tes<-dirichlet(region.ppp) # Dirichlet tesselation 

tes.poly<-as(region.tes, "SpatialPolygons") 
proj4string(tes.poly)<-CRS("+proj=merc +datum=NAD83") 
tes.poly<-spTransform(tes.poly, CRS("+proj=merc +datum=NAD83")) #planar 

plot(region.tes, main=" ") # tessellation plot
plot(region.ppp, add=TRUE, pch=".", col="darkblue", cex=2)

nnew<-100 # number of new points in the forecast 
forecasts1<-matrix(0, nrow=nnew, ncol=5) 
colnames(forecasts1)<-c("predicted y","real y","crds x","crds y", "diff")

# restoring the medoid model
dane.x<-houses.in[selector[,c1.sdm$id.med],] # we choose the data which were used in estimation of the best model
crds<-as.matrix(dane.x[,c("Longitude", "Latitude")]) # we create W
RAMSE.med.sdm<-(sum((houses.in[selector[,c1.sdm$id.med], c("Median_House_Value")]-fitted.sdm[,c1.sdm$id.med])^2)/n.row)^(0.5) 
RAMSE.med.sdm

pkt.knn<-knearneigh(crds, k=5) 
pkt.k.nb<-knn2nb(pkt.knn) 
pkt.k.sym.nb<-make.sym.nb(pkt.k.nb)
pkt.k.sym.listw<-nb2listw(pkt.k.sym.nb)

eq<-Median_House_Value ~ Median_Income + Median_Age +  Tot_Rooms + Tot_Bedrooms + Population + Households + Distance_to_coast + Distance_to_LA + Distance_to_SanDiego + Distance_to_SanJose + Distance_to_SanFrancisco


model.sdm<-lagsarlm(eq, data=dane.x, pkt.k.sym.listw, method="LU", type="mixed")
moran.test(model.sdm$residuals, pkt.k.sym.listw)

# selection of out-of-sample points for predictions
points.pred<-SpatialPoints(houses.out[1:nnew, c("Longitude", "Latitude")]) # selection of a new point
proj4string(points.pred)<-CRS("+proj=longlat +datum=NAD83") # spherical projection
points.pred<-spTransform(points.pred, CRS("+proj=merc +datum=NAD83")) 
a1<-over(points.pred, tes.poly) # assigning points to tessellation tiles

# completing the draw when NA occurs # determining the number of new points to be drawn (as from-to) 
a2<-nnew+1 # from … 
a3<-which(is.na(a1)) 
a4<-a2+length(a3)-1 # to … 
points.pred2<-SpatialPoints(houses.out[a2:a4,  c("Longitude", "Latitude")]) # new points 
proj4string(points.pred2)<-CRS("+proj=longlat +datum=NAD83") # spherical 
points.pred2<-spTransform(points.pred2, CRS("+proj=merc +datum=NAD83")) 
points.pred2 
a5<-over(points.pred2, tes.poly) # putting new points on the tile 
a5 
a1[which(is.na(a1))]<-a5 # overwriting with new points

# loop for forecasts for new points # there is a separate match for each point 
for(i in 1:nnew){ 
  # point by point - assigning new data to the old data set 
  dane.x.new<-dane.x 
  xxx<-houses.out[i, ] 
  dane.x.new[a1[i],]<-xxx 
  rownames(dane.x.new)<-1:dim(dane.x.new)[1] 
  
  # prediction for out-of-sample calibrated SDM model
  pred<-predict(model.sdm, newdata=dane.x.new, listw=pkt.k.sym.listw, legacy.mixed=TRUE) 
  
  pred[a1[i]] # prediction for a new point 
  xxx[,c("Median_House_Value")] # empirical y value of the new point 
  forecasts1[i,1]<- pred[a1[i]] # predicted value of y 
  forecasts1[i,2]<- xxx[, c("Median_House_Value")] # empirical value of y 
  forecasts1[i,3]<-xxx[, c("Longitude")] # x coordinates 
  forecasts1[i,4]<-xxx[, c("Latitude")] # y coordinates 
} 
forecasts1[,5]<-(forecasts1[,1]-forecasts1[,2])^2 
RAMSE.sdm<-(mean(forecasts1[,5]))^0.5 
head(forecasts1)
RAMSE.sdm

#########################################################
# all RAMSE together – for models and for forecasts
RAMSE.ols
RAMSE.sem
RAMSE.sar
RAMSE.sdm

RAMSE.med.ols
RAMSE.med.sem
RAMSE.med.sar
RAMSE.med.sdm

# full sample SDM ??? to RAMSE of forecast
dane.x <- houses.in
crds<-as.matrix(houses.in[ ,c("Longitude", "Latitude")])
# tessellation of observations from the selected best model
points<-data.frame(x=houses.in[ ,c("Longitude")], y=houses.in[ ,c("Latitude")]) 
points.sp<-SpatialPoints(points) # new points in sp class - spherical 
proj4string(points.sp)<-CRS("+proj=longlat +datum=NAD83") # spherical 
points.sp<-spTransform(points.sp, CRS("+proj=merc +datum=NAD83")) # planar 
region.ppp<-ppp(x=points.sp@coords[,1], y=points.sp@coords[,2], window=region.owin) # points of ppp class 
region.tes<-dirichlet(region.ppp) # Dirichlet tesselation 

tes.poly<-as(region.tes, "SpatialPolygons") 
proj4string(tes.poly)<-CRS("+proj=merc +datum=NAD83") 
tes.poly<-spTransform(tes.poly, CRS("+proj=merc +datum=NAD83")) #planar 

plot(region.tes, main=" ") # tessellation plot
plot(region.ppp, add=TRUE, pch=".", col="darkblue", cex=2)

nnew<-100 # number of new points in the forecast 
forecasts1<-matrix(0, nrow=nnew, ncol=5) 
colnames(forecasts1)<-c("predicted y","real y","crds x","crds y", "diff")

pkt.knn<-knearneigh(crds, k=5) 
pkt.k.nb<-knn2nb(pkt.knn) 
pkt.k.sym.nb<-make.sym.nb(pkt.k.nb)
pkt.k.sym.listw<-nb2listw(pkt.k.sym.nb)

eq
model.sdm.full<-lagsarlm(eq, data=dane.x, pkt.k.sym.listw, method="LU", type="mixed")
moran.test(model.sdm.full$residuals, pkt.k.sym.listw)

# selection of out-of-sample points for predictions
points.pred<-SpatialPoints(houses.out[1:nnew, c("Longitude", "Latitude")]) # selection of a new point
proj4string(points.pred)<-CRS("+proj=longlat +datum=NAD83") # spherical projection
points.pred<-spTransform(points.pred, CRS("+proj=merc +datum=NAD83")) 
a1<-over(points.pred, tes.poly) # assigning points to tessellation tiles

# completing the draw when NA occurs # determining the number of new points to be drawn (as from-to) 
a2<-nnew+1 # from … 
a3<-which(is.na(a1)) 
a4<-a2+length(a3)-1 # to … 
points.pred2<-SpatialPoints(houses.out[a2:a4, c("Longitude", "Latitude")]) # new points 
proj4string(points.pred2)<-CRS("+proj=longlat +datum=NAD83") # spherical 
points.pred2<-spTransform(points.pred2, CRS("+proj=merc +datum=NAD83")) 
points.pred2 
a5<-over(points.pred2, tes.poly) # putting new points on the tile 
a5 
a1[which(is.na(a1))]<-a5 # overwriting with new points

# loop for forecasts for new points # there is a separate match for each point 
for(i in 1:nnew){ 
  # point by point - assigning new data to the old data set 
  dane.x.new<-dane.x 
  xxx<-houses.out[i, ] 
  dane.x.new[a1[i],]<-xxx 
  rownames(dane.x.new)<-1:dim(dane.x.new)[1] 
  
  # prediction for out-of-sample calibrated SDM model
  pred<-predict(model.sdm.full, newdata=dane.x.new, listw=pkt.k.sym.listw, legacy.mixed=TRUE) 
  
  pred[a1[i]] # prediction for a new point 
  xxx[,c("Median_House_Value")] # empirical y value of the new point 
  forecasts1[i,1]<- pred[a1[i]] # predicted value of y 
  forecasts1[i,2]<- xxx[, c("Median_House_Value")] # empirical value of y 
  forecasts1[i,3]<-xxx[, c("Longitude")] # x coordinates 
  forecasts1[i,4]<-xxx[, c("Latitude")] # y coordinates 
} 
forecasts1[,5]<-(forecasts1[,1]-forecasts1[,2])^2 
RAMSE.sdm.full<-(mean(forecasts1[,5]))^0.5 
head(forecasts1)
RAMSE.sdm.full


