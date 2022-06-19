library(spdep)
library(sp)
library(spatialreg)
library(rgdal)
library(maptools)
library(spatstat)
library(GADMTools)
library(GISTools)
library(lmtest)
library(texreg)
library(broom)

# Read California shapefile
setwd('/Users/ewelinaplachimowicz/Documents/GitHub/housing-prices-spatial-econometrics/model/')

county<-readOGR(".", "CA_Counties_TIGER2016")
county<-spTransform(county, CRS("+proj=longlat +datum=NAD83"))
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

# preparing spatial weights matrix
cont.nb.county<-poly2nb(as(county, "SpatialPolygons"))
cont.listw.county<-nb2listw(cont.nb.county, style="W")

# add poly names to dataframe
houses_poly_names <- point.in.poly(houses.sp,county)
houses <- as.data.frame(houses_poly_names)

# group data to polygons names
houses <- houses %>% group_by(NAME) %>% summarise(Median_House_Value = median(Median_House_Value),
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
                                                  coords.x1 = min(coords.x1),
                                                  coords.x2 = min(coords.x2))

# drop NAs
houses <- houses %>% drop_na()

houses.sp <- houses

# class change from data.frame to sp
coordinates(houses.sp)<-c("coords.x1","coords.x2") # conversion to sp

# adding a projection to point data
proj4string(houses.sp)<-CRS("+proj=longlat +datum=NAD83") 
houses.sp<-spTransform(houses.sp, CRS("+proj=longlat +datum=NAD83")) 

# map of dependent variable
shades<-auto.shading(houses$Median_House_Value)
choropleth(county, houses$Median_House_Value, shades)
choro.legend(-130, 36, sh = shades, bty="n", title = 'Median Value in $')
title(main="Median House Value in California")

##################### First model - OLS ##################### 

# pick vars - OLS model

eq_houses <- Median_House_Value ~ Median_Income + Median_Age + Tot_Rooms + Tot_Bedrooms + Population + Households + 
  Distance_to_coast + Distance_to_LA + Distance_to_SanDiego + Distance_to_SanJose + Distance_to_SanFrancisco

model.lm.houses<-lm(eq_houses, data=houses)
summary(model.lm.houses)

eq_houses <- Median_House_Value ~ Median_Income + Median_Age + Tot_Rooms + Tot_Bedrooms + Population + Households + 
  Distance_to_coast + Distance_to_SanDiego + Distance_to_SanJose + Distance_to_SanFrancisco

model.lm.houses<-lm(eq_houses, data=houses)
summary(model.lm.houses)

eq_houses <- Median_House_Value ~ Median_Income + Median_Age + Tot_Rooms + Tot_Bedrooms + Population + 
  Distance_to_coast + Distance_to_SanDiego + Distance_to_SanJose + Distance_to_SanFrancisco

model.lm.houses<-lm(eq_houses, data=houses)
summary(model.lm.houses)

eq_houses <- Median_House_Value ~ Median_Income + Median_Age + Tot_Rooms + Tot_Bedrooms + 
  Distance_to_coast + Distance_to_SanDiego + Distance_to_SanJose + Distance_to_SanFrancisco

model.lm.houses<-lm(eq_houses, data=houses)
summary(model.lm.houses)

eq_houses <- Median_House_Value ~ Median_Income + Median_Age + Tot_Rooms + Tot_Bedrooms + 
  Distance_to_coast + Distance_to_SanDiego + Distance_to_SanJose

model.lm.houses<-lm(eq_houses, data=houses)
summary(model.lm.houses)

eq_houses <- Median_House_Value ~ Median_Income + Median_Age + Tot_Rooms + Tot_Bedrooms + 
  Distance_to_coast + Distance_to_SanJose

model.lm.houses<-lm(eq_houses, data=houses)
summary(model.lm.houses)

eq_houses <- Median_House_Value ~ Median_Income + Median_Age + Tot_Rooms + Tot_Bedrooms + Distance_to_SanJose

model.lm.houses<-lm(eq_houses, data=houses)
summary(model.lm.houses)

eq_houses <- Median_House_Value ~ Median_Income + Median_Age + Tot_Rooms + Tot_Bedrooms

model.lm.houses<-lm(eq_houses, data=houses)
summary(model.lm.houses)

# diagnostics of linear model
bptest(model.lm.houses)

# Ramsey’s test for functional form
resettest(model.lm.houses, power=2, type="regressor") 	

# spatial distribution of OLS residuals
summary(model.lm.houses$residuals)
res_houses <- model.lm.houses$residuals
brks_houses <- c(min(res_houses), mean(res_houses)-sd(res_houses), mean(res_houses), 
                 mean(res_houses)+sd(res_houses), max(res_houses))
cols_houses<-c("steelblue4", "lightskyblue", "thistle1", "plum3")
plot(county, col=cols_houses[findInterval(res_houses,brks_houses)])
title(main="Residuals in OLS model")
legend("bottomleft", legend=c("<mean-sd", "(mean-sd, mean)", "(mean, mean+sd)", ">mean+sd"), 
       leglabs(brks1), fill=cols_houses, bty="n")

# Moran test
lm.morantest(model.lm.houses, cont.listw.county) # are residuals spatially random?
moran.test(res_houses, cont.listw.county)

# For the Global Moran's I statistic, the null hypothesis states that the attribute being analyzed is randomly 
# distributed among the features in your study area; said another way, the spatial processes promoting the observed 
# pattern of values is random chance.

# test join.count for residuals (positive vs. negative)
resid_houses <- factor(cut(res_houses, breaks=c(-51000, 0, 60000), labels=c("negative","positive")))
joincount.test(resid_houses, cont.listw.county)

##################### Estimation of spatial models with 3,2, or 1 spatial components ##################### 

# Manski model
GNS_1_houses<-sacsarlm(eq_houses, data=houses, listw=cont.listw.county, type="sacmixed", method="LU")
summary(GNS_1_houses)

# SAC (or SARAR)
SAC_1_houses<-sacsarlm(eq_houses, data=houses, listw=cont.listw.county)
summary(SAC_1_houses)

# SEM / SDEM - Spatial Error Model + Spatial Durbin Error Model
SDEM_1_houses<-errorsarlm(eq_houses, data=houses, listw=cont.listw.county, etype="emixed")
summary(SDEM_1_houses)
SEM_1_houses<-errorsarlm(eq_houses, data=houses, listw=cont.listw.county)
summary(SEM_1_houses)

# SAR / SDM - Spatial Lag Model + Spatial Durbin Model

SDM_1_houses<-lagsarlm(eq_houses, data=houses, listw=cont.listw.county, type="mixed")
summary(SDM_1_houses)
SAR_1_houses<-lagsarlm(eq_houses, data=houses, listw=cont.listw.county)
summary(SAR_1_houses)

# SLX – OLS with thetas only
SLX_1_houses<-lmSLX(eq_houses, data=houses, listw=cont.listw.county)
summary(SLX_1_houses)

# OLS
OLS_1_houses<-lm(eq_houses, data=houses)
summary(OLS_1_houses)

# summary
screenreg(list(GNS_1_houses, SAC_1_houses, SDEM_1_houses, SEM_1_houses, SDM_1_houses, SAR_1_houses, SLX_1_houses, OLS_1_houses), 
          custom.model.names=c("GNS", "SAC", "SDEM", "SEM", "SDM", "SAR", "SLX", "OLS"))

################################### Direct and indirect impacts ###################################  

# estimation of spatial lag model
SAR_1_houses<-lagsarlm(eq_houses, data=houses, listw=cont.listw.county)

W.c<-as(as_dgRMatrix_listw(cont.listw.county), "CsparseMatrix") 
trMat<-trW(W.c, type="mult") 

SAR_1_imp<-impacts(SAR_1_houses, tr=trMat)
SAR_1_imp

# estimation of Spatial Durbin Model (SDM)
SDM_1<-lagsarlm(eq_houses, data=houses, listw=cont.listw.county, type="mixed") 
summary(SDM_1)

W.c<-as(as_dgRMatrix_listw(cont.listw.county), "CsparseMatrix") 
trMat<-trW(W.c, type="mult") 

SDM_1_imp<-impacts(SDM_1, tr=trMat)
SDM_1_imp

