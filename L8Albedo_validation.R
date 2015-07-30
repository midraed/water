library(raster)
library(sp)
library(rgdal)
library(proj4)


setwd("~/Documentos/Doctorado/MODIS/val.albedo/")
load("~/Documentos/Doctorado/L8METRICforR/MfR_functions.RData")

aoi <- create.aoi(topleft = c(370470, -3537625), bottomright = c(591337, -3799132))
plot(aoi)

# files <- list.files(path = getwd(), pattern = "tif") #2.8 Gb
# for(i in 1:length(files)){
#   temp <- raster(files[i])
#   temp <- crop(temp, aoi)
#   writeRaster(temp, filename = paste("CROP",files[i], sep=""), overwrite=TRUE )
# } # 70Mb a bit better :P
removeTmpFiles(h=0)

val.sites <- read.csv("~/Documentos/Doctorado/L8METRICforR/test.points.csv")
coordinates(val.sites) <- ~ X + Y

#####

setwd("~/Documentos/Doctorado/MODIS/albedo.valset/")

### Albedo from MOD09, using Tasumi coefficients #####
albedo.coeff <- c(0.215, 0.215, 0.242, 0.129, 0.101, 0.062, 0.036)

files <- list.files(pattern="CROPMOD09_2013313.sur_refl_b0") 
stack1 <- list()
for(i in 1:length(files)){
  stack1[[i]] <- raster(files[i])/10000*albedo.coeff[i]
  }
temp <- do.call(stack, stack1) 
MOD09.2013313 <- stackApply(temp, indices = c(1,1,1,1,1,1,1), 
                            fun=sum, filename="albedo_MOD09.2013313.tif")


files <- list.files(pattern="CROPMOD09_2013337.sur_refl_b0") 
stack1 <- list()
for(i in 1:length(files)){
  stack1[[i]] <- raster(files[i])/10000*albedo.coeff[i]
}
temp <- do.call(stack, stack1) 
MOD09.2013337 <- stackApply(temp, indices = c(1,1,1,1,1,1,1), 
                            fun=sum, filename="albedo_MOD09.2013337.tif")


files <- list.files(pattern="CROPMOD09_2013353.sur_refl_b0") 
stack1 <- list()
for(i in 1:length(files)){
  stack1[[i]] <- raster(files[i])/10000*albedo.coeff[i]
}
temp <- do.call(stack, stack1) 
MOD09.2013353 <- stackApply(temp, indices = c(1,1,1,1,1,1,1), 
                            fun=sum, filename="albedo_MOD09.2013353.tif")

files <- list.files(pattern="CROPMOD09_2014001.sur_refl_b0") 
stack1 <- list()
for(i in 1:length(files)){
  stack1[[i]] <- raster(files[i])/10000*albedo.coeff[i]
}
temp <- do.call(stack, stack1) 
MOD09.2014001 <- stackApply(temp, indices = c(1,1,1,1,1,1,1), 
                            fun=sum, filename="albedo_MOD09.2014001.tif")
  
files <- list.files(pattern="CROPMOD09_2014025.sur_refl_b0") 
stack1 <- list()
for(i in 1:length(files)){
  stack1[[i]] <- raster(files[i])/10000*albedo.coeff[i]
}
temp <- do.call(stack, stack1) 
MOD09.2014025 <- stackApply(temp, indices = c(1,1,1,1,1,1,1), 
                            fun=sum, filename="albedo_MOD09.2014025.tif")

### Albedo from L8 using Olmedo coefficients ####

L8

