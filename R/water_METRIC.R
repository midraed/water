#' Estimates Net Radiation as in METRIC Model
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
METRIC.Rn <- function(image.DN, DEM, WeatherStation, aoi){
  path=getwd()
  surface.model <-METRICtopo(DEM)
  solar.angles.r <- solarAngles(surface.model = surface.model)
  Rs.inc <- incSWradiation(surface.model = surface.model, solar.angles = solar.angles.r, WeatherStation = WeatherStation)
  image.TOAr <- calcTOAr(image.DN = image.DN, incidence.rel = solar.angles.r$incidence.rel)
  image.SR <- calcSR(path=path, image.TOAr=image.TOAr, 
                      surface.model=surface.model, 
                      incidence.hor = solar.angles.r$incidence.hor, 
                      WeatherStation=WeatherStation, sat="auto", ESPA = F)
  albedo <- albedo(image.SR = image.SR, sat="auto")
  LAI <- LAI(method = "metric2010", image = image.TOAr, L=0)
  Ts <- surfaceTemperature(sat = "auto", LAI=LAI)
  Rl.out <- outLWradiation(LAI = LAI, Ts=Ts)
  Rl.inc <- incLWradiation(WeatherStation,DEM = surface.model$DEM, solar.angles = solar.angles.r)
  Rn <- netRadiation(LAI, albedo, Rs.inc, Rl.inc, Rl.out)
  plot(Rn, main="METRIC Net Radiation")
  return(Rn)
}

#' Estimates Net Radiation as in METRIC Model
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
METRIC.G <- function(Rn, DEM, image.DN, WeatherStation=WeatherStation){
  path=getwd()
  surface.model <-METRICtopo(DEM)
  solar.angles.r <- solarAngles(surface.model = surface.model)
  image.TOAr <- calcTOAr(image.DN = image.DN, incidence.rel = solar.angles.r$incidence.rel)
  image.SR <- calcSR(path=path, image.TOAr=image.TOAr, 
                      surface.model=surface.model, 
                      incidence.hor = solar.angles.r$incidence.hor, 
                      WeatherStation=WeatherStation, sat="auto", ESPA = F)
  albedo <- albedo(image.SR = image.SR, sat="auto")
  Ts <- surfaceTemperature(sat = "auto", image.TOAr = image.TOAr)
  G <- soilHeatFlux(image = image.SR, Ts=Ts,albedo=albedo, Rn)
}


#' Estimates Energy Balance using METRIC2010 Model
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
METRIC.EB <- function(image.DN, DEM, WeatherStation, aoi, MTL, sat, plain=TRUE,
                      thermalband,
                      alb.coeff = "Tasumi", LAI.method = "metric2010", 
                      Zom.method = "short.crops", anchors.method = "CITRA-MCB",
                      ETp.coef= 1.05, Z.om.ws=0.0018, ESPA = FALSE){
  path=getwd()
  pb <- txtProgressBar(min = 0, max = 100, style = 3)
  if(plain==TRUE){
    DEM <- raster(image.DN[[1]])
    values(DEM) <- WeatherStation$location$elev
  } else {DEM <- prepareSRTMdata(extent = image.DN)}
  setTxtProgressBar(pb, 3)
  surface.model <-METRICtopo(DEM)
  solar.angles.r <- solarAngles(surface.model = surface.model, 
                                WeatherStation = WeatherStation, MTL = MTL)
  Rs.inc <- incSWradiation(surface.model = surface.model, 
                           solar.angles = solar.angles.r, 
                           WeatherStation = WeatherStation)
  image.TOAr <- calcTOAr(image.DN = image.DN, sat=sat,
                         incidence.rel = solar.angles.r$incidence.rel)
  image.SR <- calcSR(image.TOAr=image.TOAr, sat = sat, 
                     surface.model=surface.model, 
                     incidence.hor = solar.angles.r$incidence.hor, 
                     WeatherStation=WeatherStation, ESPA = F)
  albedo <- albedo(image.SR = image.SR,  coeff=alb.coeff, sat = sat)
  setTxtProgressBar(pb, 6)
  LAI <- LAI(method = LAI.method, image = image.TOAr, L=0.1, sat = sat)
  Ts <- surfaceTemperature(LAI=LAI, sat = sat, thermalband = thermalband,
                           WeatherStation = WeatherStation)
  setTxtProgressBar(pb, 35)
  Rl.out <- outLWradiation(LAI = LAI, Ts=Ts)
  Rl.inc <- incLWradiation(WeatherStation,DEM = surface.model$DEM, 
                           solar.angles = solar.angles.r, Ts= Ts)
  Rn <- netRadiation(LAI, albedo, Rs.inc, Rl.inc, Rl.out)
  setTxtProgressBar(pb, 40)
  G <- soilHeatFlux(image = image.SR, Ts=Ts,albedo=albedo, 
                    Rn=Rn, image.SR, LAI=LAI, sat = sat)
  Z.om <- momentumRoughnessLength(LAI=LAI, mountainous = TRUE, 
                                  method = Zom.method, 
                                  surface.model = surface.model)
  hot.and.cold <- calcAnchors(image = image.TOAr, Ts = Ts, LAI = LAI, plots = F,
                              albedo = albedo, Z.om = Z.om, n = 1, 
                              anchors.method = anchors.method, sat = sat, 
                              deltaTemp = 5, verbose = FALSE)
  setTxtProgressBar(pb, 45)
  on.meta <-  TRUE
  H <- calcH(anchors = hot.and.cold, Ts = Ts, Z.om = Z.om, 
             WeatherStation = WeatherStation, ETp.coef = ETp.coef, sat = sat, 
             Z.om.ws = Z.om.ws, DEM = DEM, Rn = Rn, G = G, verbose = FALSE)
  setTxtProgressBar(pb, 99)
  H <-  H$H
  LE <- Rn - G - H
  EB <- stack(Rn, G, H, LE, Ts)
  EB <- saveLoadClean(imagestack = EB,
                stack.names = c("NetRadiation", "SoilHeat", "SensibleHeat", 
                                "LatentHeat", "surfaceTemperature"), 
                file = "EB", overwrite=TRUE)
  return(EB)
}
  
