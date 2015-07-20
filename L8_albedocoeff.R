library(raster)
setwd("~/Documentos/Doctorado/")

b2 <- raster("ESPA/LC82320832013319-SC20150618080812_nov15/LC82320832013319LGN00_sr_band2.tif")
b3 <- raster("ESPA/LC82320832013319-SC20150618080812_nov15/LC82320832013319LGN00_sr_band3.tif")
b4 <- raster("ESPA/LC82320832013319-SC20150618080812_nov15/LC82320832013319LGN00_sr_band4.tif")
b5 <- raster("ESPA/LC82320832013319-SC20150618080812_nov15/LC82320832013319LGN00_sr_band5.tif")
b6 <- raster("ESPA/LC82320832013319-SC20150618080812_nov15/LC82320832013319LGN00_sr_band6.tif")
b7 <- raster("ESPA/LC82320832013319-SC20150618080812_nov15/LC82320832013319LGN00_sr_band7.tif")

l8image <- stack(b2,b3,b4,b5,b6,b7)


### DEM
files <- list.files(path= "L8METRICforR/",  pattern=paste("^[sn]\\d{2}_[we]\\d{3}_1arc_v3.", "tif", "$", sep="")) 
stack1 <- list()
for(i in 1:length(files)){
  stack1[[i]] <- raster(paste("L8METRICforR/", files[i], sep=""))}
stack1$fun <- mean
SRTMmosaic <- do.call(mosaic, stack1)
destino  <-  projectExtent(l8image, l8image@crs)
mosaicp <- projectRaster(SRTMmosaic, destino)

samples <- read.csv("L8METRICforR/l8albedo_groundsamples.csv")
coordinates(samples) <- ~X+Y
plot(samples)
plot(b4)
points(samples)
samples$ALT <- extract(mosaicp, samples)

samples@data
names(l8image) <- c("b2", "b3", "b4", "b5", "b6", "b7")
samples@data <- cbind(samples@data, extract(l8image, samples))
samples@data[,4:9] <- samples@data[,4:9] * 0.0001
ts.plot(t(samples@data[,4:9]), col=rainbow(20))
samples$SUM <- rowSums(samples@data[,4:9])

### As in Tasumi we are going to arbitraliry divide 
### midway between band edges to cover regions between 
### bands and calculate LOb - UOb
# b2 | 0.45 - 0.51 | 0.300 - 0.520 | 0.220
# b3 | 0.53 - 0.59 | 0.520 - 0.615 | 0.095
# b4 | 0.64 - 0.67 | 0.615 - 0.760 | 0.145
# b5 | 0.85 - 0.88 | 0.760 - 1.225 | 0.465
# b6 | 1.57 - 1.65 | 1.225 - 1.880 | 0.655
# b7 | 2.11 - 2.29 | 1.880 - 4.000 | 2.120

samples$w2 <- (samples$b2 * 0.220) / (samples$SUM * 3.7)
samples$w3 <- (samples$b3 * 0.095) / (samples$SUM * 3.7)
samples$w4 <- (samples$b4 * 0.145) / (samples$SUM * 3.7)
samples$w5 <- (samples$b5 * 0.465) / (samples$SUM * 3.7)
samples$w6 <- (samples$b6 * 0.655) / (samples$SUM * 3.7)
samples$w7 <- (samples$b7 * 2.120) / (samples$SUM * 3.7)

colSums(samples@data[-c(2,6,16),11:16])/20

# > colSums(samples@data[-c(2,6,16),11:16])/20
# w2          w3          w4          w5          w6          w7 
# 0.003871170 0.002449617 0.004805069 0.024487127 0.035110190 0.099525157

