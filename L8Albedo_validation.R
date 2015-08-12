library(raster)
library(sp)
library(rgdal)
library(proj4)

rasterOptions(timer=TRUE)

setwd("~/Documentos/Doctorado/MODIS/val.albedo/")
load("~/Documentos/Doctorado/L8METRICforR/MfR_functions.RData")

aoi <- create.aoi(topleft = c(493350, -3592750), bottomright = c(557200, -3700000))
plot(aoi)

### Albedo from MOD09, using Tasumi coefficients #####
setwd("~/Documentos/Doctorado/MODIS/albedo.valset/")



val.sites <- read.csv("~/Documentos/Doctorado/L8METRICforR/test.modisgrid.ON.csv")
coordinates(val.sites) <- ~ X + Y
points(val.sites)

### MCD43####

MCD43BSASW313 <- raster("~/Documentos/Doctorado/MODIS/albedo.valset/CROPMCD43_2013313.Albedo_BSA_shortwave.tif")*0.001
MCD43BSASW313 <- crop(MCD43BSASW313, aoi)
plot(MCD43BSASW313)

MCD43BSASW337 <- raster("~/Documentos/Doctorado/MODIS/albedo.valset/CROPMCD43_2013337.Albedo_BSA_shortwave.tif")*0.001
MCD43BSASW337 <- crop(MCD43BSASW337, aoi)
plot(MCD43BSASW337)


### MOD09 ####
setwd("~/Documentos/Doctorado/MODIS/albedo.valset/")
albedo.coeff <- c(0.215, 0.215, 0.242, 0.129, 0.101, 0.062, 0.036)
files <- list.files(pattern="CROPMOD09_2013313.sur_refl_b0") 
stack1 <- list()
for(i in 1:length(files)){
  stack1[[i]] <- raster(files[i])/10000*albedo.coeff[i]
  }
temp <- do.call(stack, stack1) 
MOD09.2013313 <- stackApply(temp, indices = c(1,1,1,1,1,1,1), 
                            fun=sum, filename="albedo_MOD09.2013313.tif", overwrite=TRUE)
MOD09.2013313 <- crop(MOD09.2013313, aoi)
plot(MOD09.2013313)
MOD09.2013313 <- projectRaster(MOD09.2013313, destino, method = "ngb")


files <- list.files(pattern="CROPMOD09_2013337.sur_refl_b0") 
stack1 <- list()
for(i in 1:length(files)){
  stack1[[i]] <- raster(files[i])/10000*albedo.coeff[i]
}
temp <- do.call(stack, stack1) 
MOD09.2013337 <- stackApply(temp, indices = c(1,1,1,1,1,1,1), 
                            fun=sum, filename="albedo_MOD09.2013337.tif", overwrite=TRUE)
MOD09.2013337 <- crop(MOD09.2013337, aoi)
plot(MOD09.2013337)
MOD09.2013337 <- projectRaster(MOD09.2013337, destino, method = "ngb")

### Albedo from L8 using Tasumio coefficients ####
setwd("~/Documentos/Doctorado/")

destino  <-  projectExtent(MCD43BSASW313, MCD43BSASW313@crs)

wbTasumi <- c(0.254, 0.149, 0.147, 0.311, 0.103, 0.036)

removeTmpFiles(h=0)
path = "ESPA/LC82320832013319-SC20150618080812_nov15/"
wb <- wbTasumi # Calculated using SMARTS for Kimberly and Direct normal irradiance
srb2 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band2.tif"), sep="")[1]), fun=function(x){x /10000*wb[1]})
srb3 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band3.tif"), sep="")[1]), fun=function(x){x /10000*wb[2]})
srb4 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band4.tif"), sep="")[1]), fun=function(x){x /10000*wb[3]})
srb5 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band5.tif"), sep="")[1]), fun=function(x){x /10000*wb[4]})
srb6 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band6.tif"), sep="")[1]), fun=function(x){x /10000*wb[5]})
srb7 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band7.tif"), sep="")[1]), fun=function(x){x /10000*wb[6]})
l8.albedo <-  stack(srb2, srb3, srb4, srb5, srb6, srb7)
l8.albedo <- stackApply(l8.albedo, indices = c(1,1,1,1,1,1), fun=sum)
l8.albedo <- crop(l8.albedo, aoi)
L82013319.Tasumi <- projectRaster(l8.albedo, destino, method = "bilinear", 
                           filename = "AlbedoL8/Albedo_L82013319_MODGRID_Tasumi.tif", overwrite=T)
plot(L82013319.Tasumi)

removeTmpFiles(h=0)
path = "ESPA/LC82320832013335-SC20150618080841_dec01/"
wb <- wbTasumi # Calculated using SMARTS for Kimberly and Direct normal irradiance
srb2 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band2.tif"), sep="")[1]), fun=function(x){x /10000*wb[1]})
srb3 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band3.tif"), sep="")[1]), fun=function(x){x /10000*wb[2]})
srb4 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band4.tif"), sep="")[1]), fun=function(x){x /10000*wb[3]})
srb5 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band5.tif"), sep="")[1]), fun=function(x){x /10000*wb[4]})
srb6 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band6.tif"), sep="")[1]), fun=function(x){x /10000*wb[5]})
srb7 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band7.tif"), sep="")[1]), fun=function(x){x /10000*wb[6]})
l8.albedo <-  stack(srb2, srb3, srb4, srb5, srb6, srb7)
l8.albedo <- stackApply(l8.albedo, indices = c(1,1,1,1,1,1), fun=sum)
l8.albedo <- crop(l8.albedo, aoi)
L82013335.Tasumi <- projectRaster(l8.albedo, destino, method = "bilinear", 
                           filename = "AlbedoL8/Albedo_L82013335_MODGRID_Tasumi.tif", overwrite=T)
plot(L82013335.Tasumi)

### Albedo from L8 using my coefficients ####
setwd("~/Documentos/Doctorado/")

destino  <-  projectExtent(MCD43BSASW313, MCD43BSASW313@crs)

wbOlmedo <- c(0.246, 0.146, 0.191, 0.304, 0.105, 0.008)

removeTmpFiles(h=0)
path = "ESPA/LC82320832013319-SC20150618080812_nov15/"
wb <- wbOlmedo # Calculated using SMARTS for Kimberly and Direct normal irradiance
srb2 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band2.tif"), sep="")[1]), fun=function(x){x /10000*wb[1]})
srb3 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band3.tif"), sep="")[1]), fun=function(x){x /10000*wb[2]})
srb4 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band4.tif"), sep="")[1]), fun=function(x){x /10000*wb[3]})
srb5 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band5.tif"), sep="")[1]), fun=function(x){x /10000*wb[4]})
srb6 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band6.tif"), sep="")[1]), fun=function(x){x /10000*wb[5]})
srb7 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band7.tif"), sep="")[1]), fun=function(x){x /10000*wb[6]})
l8.albedo <-  stack(srb2, srb3, srb4, srb5, srb6, srb7)
l8.albedo <- stackApply(l8.albedo, indices = c(1,1,1,1,1,1), fun=sum)
l8.albedo <- crop(l8.albedo, aoi)
L82013319 <- projectRaster(l8.albedo, destino, method = "bilinear", 
                           filename = "AlbedoL8/Albedo_L82013319_MODGRID.tif", overwrite=T)
plot(L82013319)

removeTmpFiles(h=0)
path = "ESPA/LC82320832013335-SC20150618080841_dec01/"
wb <- wbOlmedo # Calculated using SMARTS for Kimberly and Direct normal irradiance
srb2 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band2.tif"), sep="")[1]), fun=function(x){x /10000*wb[1]})
srb3 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band3.tif"), sep="")[1]), fun=function(x){x /10000*wb[2]})
srb4 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band4.tif"), sep="")[1]), fun=function(x){x /10000*wb[3]})
srb5 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band5.tif"), sep="")[1]), fun=function(x){x /10000*wb[4]})
srb6 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band6.tif"), sep="")[1]), fun=function(x){x /10000*wb[5]})
srb7 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band7.tif"), sep="")[1]), fun=function(x){x /10000*wb[6]})
l8.albedo <-  stack(srb2, srb3, srb4, srb5, srb6, srb7)
l8.albedo <- stackApply(l8.albedo, indices = c(1,1,1,1,1,1), fun=sum)
l8.albedo <- crop(l8.albedo, aoi)
L82013335 <- projectRaster(l8.albedo, destino, method = "bilinear", 
                           filename = "AlbedoL8/Albedo_L82013335_MODGRID.tif", overwrite=T)
plot(L82013335)


### Validation ########################

Albedo <- stack(MCD43BSASW313, MCD43BSASW337, MOD09.2013313, MOD09.2013337,
                L82013319.Tasumi, L82013335.Tasumi, L82013319, L82013335)


library(calibrate)
val.sites <- read.csv("~/Documentos/Doctorado/AlbedoL8/test.modisgrid.ON2.csv")
val.sites <- val.sites[!val.sites$landcover=="city",]
coordinates(val.sites) <- ~ X + Y

val.sites <- cbind(val.sites@data, extract(Albedo, val.sites))

landcover <- c(as.character(val.sites$landcover),as.character(val.sites$landcover))
MCD43 <- c(val.sites$CROPMCD43_2013313.Albedo_BSA_shortwave, val.sites$CROPMCD43_2013337.Albedo_BSA_shortwave)
MOD09 <- c(val.sites$albedo_MOD09.2013313, val.sites$albedo_MOD09.2013337)
L8TAS <- c(val.sites$Albedo_L82013319_MODGRID_Tasumi, val.sites$Albedo_L82013335_MODGRID_Tasumi)
L8OLM <- c(val.sites$Albedo_L82013319_MODGRID, val.sites$Albedo_L82013335_MODGRID)



plot( L8OLM,MCD43, cex=0.5)
textxy( L8OLM, MCD43, landcover)
textxy( L8OLM,MCD43, c(1:58))
abline(0,1)
model.MCD43L8OLM <-  lm(MCD43~L8OLM)
summary(model.MCD43L8OLM)
plot(model.MCD43L8OLM)

abline(model.MCD43L8OLM,col="red")
newx<-seq(0.05,0.3, length.out = 59)
prd<-predict(model.MCD43L8OLM,newdata=data.frame(L8OLM=newx),interval = c("confidence"), 
             level = 0.90,type="response")
lines(newx,prd[,2],col="red",lty=2)
lines(newx,prd[,3],col="red",lty=2)

plot(MCD43, L8TAS, cex=0.5)
abline(0,1)
summary(lm(MCD43~L8TAS))

plot(MOD09, L8OLM, cex=0.5)
abline(0,1)
summary(lm(MOD09 ~ L8OLM))
