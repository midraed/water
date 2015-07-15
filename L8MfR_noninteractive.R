### This script uses a lot of temp raster files
### This could fill you temp folder making this script to fail.
### You should watch the free space of your temp folder
### Your temp folder it's located on: rasterOptions()$tmpdir

library(raster)
library(sp)
library(rgdal)
library(proj4)


setwd("~/Documentos/Doctorado/")
load("L8METRICforR/MfR_functions.RData")

lujan <- create.aoi(topleft = c(507022, -3645139), bottomright = c(518821, -3658440))
plot(lujan)

# load L8 data, any L8-band it's ok.
raw.image <-  load_L8data("ESPA/LC82320832013319-SC20150618080812_nov15/LC82320832013319LGN00_B1.TIF", aoi= lujan)
plot(raw.image[[1]])

# check (but not get) SRTM 1 arc-second data
checkSRTMgrids(raw.image)

# This may takes some (many) minutes with a full scene. 
surface.model <- prepareSRTMdata(path = "L8METRICforR/", format = "tif", extent = raw.image)
plot(surface.model)


### Incoming Solar Radiation
solar.angles.r <- solar.angles("ESPA/LC82320832013319-SC20150618080812_nov15/LC82320832013319LGN00_MTL.txt",
                           raw.image = raw.image, slope = surface.model$Slope, aspect = surface.model$Aspect)
plot(solar.angles.r)

tau.sw <- sw.trasmisivity(ea= 0.38, dem = surface.model$DEM, incidence.hor = solar.angles.r$incidence.hor)
plot(tau.sw)

Rs.inc <- incoming.solar.radiation(solar.angles.r$incidence.rel, tau.sw, DOY = 319)
plot(Rs.inc)

### Surface albedo
l8.albedo <- albedo(path = "ESPA/LC82320832013319-SC20150618080812_nov15/", aoi = lujan)
plot(l8.albedo)

### Outgoing Long-Wave Radiation
LAI <- LAI.metric(path="ESPA/LC82320832013319-SC20150618080812_nov15/", L=0.5, aoi=lujan)
plot(LAI)

Rl.out <- outgoing.lw.radiation(path = "ESPA/LC82320832013319-SC20150618080812_nov15/", LAI, lujan)
plot(Rl.out)

### Incoming Long-Wave Radiation 
Rl.inc <- incoming.lw.radiation(298.15, DEM = surface.model$DEM, tau.sw, lujan)
plot(Rl.inc)

### Net Radiation
surf.emissivity <- 0.95 + 0.01 * LAI 
Rn <- Rs.inc - l8.albedo*Rs.inc + Rl.inc - Rl.out - (1-surf.emissivity)*Rl.inc

plot(Rn, main="Net Radiation")

### Soil Heat Flux
G <- soil.heat.flux1(path="ESPA/LC82320832013319-SC20150618080812_nov15/", l8.albedo, Rn, lujan)
plot(G, main="Soil Heat Flux")

###  
