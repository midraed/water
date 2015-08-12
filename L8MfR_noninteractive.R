### This script uses a lot of temp raster files
### This could fill you temp folder making this script to fail.
### You should watch the free space of your temp folder
### Your temp folder it's located on: rasterOptions()$tmpdir
start.all <- Sys.time ()
library(raster)
library(sp)
library(rgdal)
library(proj4)

WeatherStation  <- c(wind=2.3, ETr= 8.94, ea= 0.5322, Ta=294.05)

setwd("~/Documentos/Doctorado/")
load("L8METRICforR/MfR_functions.RData")

lujan <- create.aoi(topleft = c(507022, -3645139), bottomright = c(518821, -3658440))
lujan <- create.aoi(topleft = c(493350, -3592750), bottomright = c(557200, -3700000))
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

removeTmpFiles(h=0)
### Surface albedo
l8.albedo <- albedo(path = "ESPA/LC82320832013319-SC20150618080812_nov15/", aoi = lujan)
plot(l8.albedo)

Ts <- surface.temperature(path="ESPA/LC82320832013319-SC20150618080812_nov15/", lujan)
plot(Ts-273)

### Outgoing Long-Wave Radiation
LAI <- LAI.from.L8(path="ESPA/LC82320832013319-SC20150618080812_nov15/", L=0.5, aoi=lujan)
plot(LAI)

Rl.out <- outgoing.lw.radiation(path = "ESPA/LC82320832013319-SC20150618080812_nov15/", LAI, lujan)
plot(Rl.out)



### Incoming Long-Wave Radiation 
Rl.inc <- incoming.lw.radiation(WeatherStation["Ta"], DEM = surface.model$DEM, tau.sw, lujan)
plot(Rl.inc)

### Net Radiation
surf.emissivity <- 0.95 + 0.01 * LAI 
Rn <- Rs.inc - l8.albedo*Rs.inc + Rl.inc - Rl.out - (1-surf.emissivity)*Rl.inc

plot(Rn, main="Net Radiation")

### Soil Heat Flux
G <- soil.heat.flux1(path="ESPA/LC82320832013319-SC20150618080812_nov15/", l8.albedo, Rn, lujan)
plot(G, main="Soil Heat Flux")



###  Sensible Heat Flux
Z.om <- momentum.roughness.length(method="short.crops", LAI=LAI, mountainous = TRUE, surf.model = surface.model)
plot(Z.om)

print (Sys.time () - start.all)

delta <- aerodynamic.transport.i(Ts = Ts, LAI = LAI, n = 1, wind = WeatherStation["wind"], 
                                 ETr = WeatherStation["ETr"], ETp.coef = 1.2, Z.om.ws = 0.03, 
                                 anchors.method = "random", height.ws = 2, mountainous = TRUE, 
                                 elev.ws = 920, DEM = surface.model$DEM, Rn = Rn, G = G, 
                                 plots = TRUE)

H <- sensible.heat.flux(Ts, dem=surface.model$DEM, dT, r.ah)
plot(H, main="sensible.heat.flux")



###################### Evapotranspiration

LE <- Rn - G - H
plot(LE, main="latent energy")

ET.inst <- 3600*LE/(latent.heat.vaporization*1000)
plot(ET.inst)

ETr.F <- ET.inst/ETr
plot(ETr.F)

ET_daily <- ETrF * (ETr*12)