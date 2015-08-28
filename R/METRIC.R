#' Estimates Net Radiation as in METRIC Model
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
METRIC.Rn <- function(path=getwd(), image.DN, DEM, WeatherStation, aoi){
  surface.model <-METRIC.topo(DEM)
  solar.angles.r <- solar.angles(surface.model = surface.model)
  Rs.inc <- incoming.sw.radiation(surface.model, solar.angles.r, WeatherStation)
  image.TOAr <- calc.TOAr(image.DN = image.DN, incidence.rel = solar.angles.r$incidence.rel)
  image.SR <- calc.SR(path=path, image.TOAr=image.TOAr, 
                      surface.model=surface.model, 
                      incidence.hor = solar.angles.r$incidence.hor, 
                      WeatherStation=WeatherStation, sat="auto", ESPA = F)
  albedo <- albedo(image.SR = image.SR, sat="auto")
  LAI <- LAI(method = "metric2010", image = image.TOAr, L=0)
  Ts <- surface.temperature(sat = "auto", LAI=LAI)
  Rl.out <- outgoing.lw.radiation(LAI = LAI, Ts=Ts)
  Rl.inc <- incoming.lw.radiation(WeatherStation,DEM = surface.model$DEM, solar.angles = solar.angles.r)
  surf.emissivity <- 0.95 + 0.01 * LAI 
  Rn <- Rs.inc - albedo*Rs.inc + Rl.inc - Rl.out - (1-surf.emissivity)*Rl.inc
  plot(Rn, main="METRIC Net Radiation")
  return(Rn)
}

#' Estimates Net Radiation as in METRIC Model
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
METRIC.G <- function(Rn, DEM, image.DN){
  surface.model <-METRIC.topo(DEM)
  solar.angles.r <- solar.angles(surface.model = surface.model)
  image.TOAr <- calc.TOAr(image.DN = image.DN, incidence.rel = solar.angles.r$incidence.rel)
  image.SR <- calc.SR(path=path, image.TOAr=image.TOAr, 
                      surface.model=surface.model, 
                      incidence.hor = solar.angles.r$incidence.hor, 
                      WeatherStation=WeatherStation, sat="auto", ESPA = F)
  albedo <- albedo(image.SR = image.SR, sat="auto")
  Ts <- surface.temperature(sat = "auto", image.TOAr = image.TOAr)
  G <- soil.heat.flux(image = image.SR, Ts=Ts,albedo=albedo, Rn)
}