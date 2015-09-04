start.all <- Sys.time ()

library(water)

setwd("~/Documentos/Doctorado")

results = "METRIC_results/"

WeatherStation  <- data.frame(hours=11.5,
                     wind=2.3,
                     RH=24.8, 
                     Ta=21.4,
                     Gr.Rad=1083, 
                     height=2, 
                     lat.ws=-32.987, 
                     long.ws=68.738, 
                     elev.ws=764,
                     hours=10.5)




lujan <- create.aoi(topleft = c(493350, -3592750), bottomright = c(557200, -3700000))
raw.image <-  load.image.DN(path="ESPA/LC82320832013319-SC20150618080812_nov15/", result.folder = results, aoi = lujan)
surface.model <-METRIC.topo("old_L8METRICforR/DEM.ON.tif")
solar.angles.r <- solar.angles(L8MTL = "ESPA/LC82320832013319-SC20150618080812_nov15/LC82320832013319LGN00_MTL.txt",
                               raw.image = raw.image, slope = surface.model$Slope, aspect = surface.model$Aspect)
tau.sw <- sw.trasmisivity(ea= 0.38, dem = surface.model$DEM, incidence.hor = solar.angles.r$incidence.hor)
Rs.inc <- incoming.solar.radiation(solar.angles.r$incidence.rel, tau.sw, DOY = 319)
l8.albedo <- albedo(path = "ESPA/LC82320832013319-SC20150618080812_nov15/", coeff = "Tasumi", aoi = lujan)
Ts <- surface.temperature(path="ESPA/LC82320832013319-SC20150618080812_nov15/", lujan)
LAI <- LAI.from.L8(path="ESPA/LC82320832013319-SC20150618080812_nov15/", L=0.5, method="metric", aoi=lujan)
Rl.out <- outgoing.lw.radiation(path = "ESPA/LC82320832013319-SC20150618080812_nov15/", LAI, lujan)
Rl.inc <- incoming.lw.radiation(WeatherStation$Ta, DEM = surface.model$DEM, tau.sw, lujan)
surf.emissivity <- 0.95 + 0.01 * LAI 
Rn <- Rs.inc - l8.albedo*Rs.inc + Rl.inc - Rl.out - (1-surf.emissivity)*Rl.inc
G <- soil.heat.flux1(path="ESPA/LC82320832013319-SC20150618080812_nov15/", l8.albedo, Rn, lujan)
Z.om <- momentum.roughness.length(method="short.crops", LAI=LAI, mountainous = TRUE, surf.model = surface.model)
ETo.hourly <- ETo.PM.hourly(WeatherStation, hours=11.5, DOY=319)
H.and.dT <- H(Ts = Ts, LAI = LAI, n = 1, WeatherStation = WeatherStation, ETo.hourly = ETo.hourly, ETp.coef = 1.2, Z.om.ws = 0.03, anchors.method = "random",  mountainous = TRUE, DEM = surface.model$DEM, Rn = Rn, G = G, plots = TRUE)

LE <- Rn - G - H.and.dT$H
ET.inst <- 3600*LE/(latent.heat.vaporization*1000)
###################### Evapotranspiration
LE <- Rn - G - H.and.dT$H
ET.inst <- 3600*LE/((2.501 - 0.00236 * (Ts - 273.15)) * (1e6))
ET.daily <- 7
ETo.hourly <- ETo.PM.hourly(WeatherStation, WeatherStation$hours, WeatherStation$DOY)
ETr.F <- ETo.hourly/ET.daily

print (Sys.time () - start.all)


