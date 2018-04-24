#' Estimates Net Radiation as in METRIC Model
#' @param image.DN         raw imagen in digital counts to evaluate
#' @param WeatherStation   Weather Station data, can be a waterWeatherStation 
#' object
#' @param MTL              Landsat metadata file
#' @param sat              Landsat satellite version. "L7" or "L8"
#' @param thermalband      Landsat low gain thermalband
#' @param plain            Logical. If TRUE surface is assumed plain
#' @param DEM              Digital Elevation Model of the study area. Not needed
#' if plain = TRUE
#' @param aoi              SpatialPolygon object with limits of Area of interest
#' @param alb.coeff        coefficient to transform narrow to broad band albedo.
#' See Details.
#' @param LAI.method       Method used to estimate LAI from spectral data. 
#' See Details.
#' @family METRIC model functions
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
METRIC.Rn <- function(image.DN, WeatherStation, MTL, sat = "auto", thermalband, 
                      alb.coeff = "Tasumi", LAI.method = "metric2010", 
                      plain = TRUE, DEM, aoi){
  path=getwd()
  #pb <- txtProgressBar(min = 0, max = 100, style = 3)
  if(plain==TRUE){
    DEM <- raster(image.DN[[1]])
    values(DEM) <- WeatherStation$location$elev
  }
  surface.model <-METRICtopo(DEM)
  #setTxtProgressBar(pb, 3)
  solar.angles.r <- solarAngles(surface.model = surface.model, 
                                WeatherStation = WeatherStation, MTL = MTL)
  Rs.inc <- incSWradiation(surface.model = surface.model, 
                           solar.angles = solar.angles.r, 
                           WeatherStation = WeatherStation)
  if(sat=="L7" | sat=="L8"){
  image.TOAr <- calcTOAr(image.DN = image.DN, sat=sat, MTL = MTL, 
                         incidence.rel = solar.angles.r$incidence.rel)
  image.SR <- calcSR(image.TOAr=image.TOAr, sat = sat, 
                     surface.model=surface.model, 
                     incidence.hor = solar.angles.r$incidence.hor, 
                     WeatherStation=WeatherStation)}
  if(sat=="MODIS"){image.SR <- image.DN}
  albedo <- albedo(image.SR = image.SR,  coeff=alb.coeff, sat=sat)
  #setTxtProgressBar(pb, 6)
  if(sat=="MODIS"){image.TOAr <- image.DN} # Only used for LAI estimation,
                                           # and some LAI models, use SR
  LAI <- LAI(method = LAI.method, image = image.TOAr, L=0.1)
  if(sat=="L7" | sat=="L8"){
  Ts <- surfaceTemperature(LAI=LAI, sat = sat, thermalband = thermalband,
                           WeatherStation = WeatherStation)}
  if(sat=="MODIS"){Ts <- image.DN$LST}
  #setTxtProgressBar(pb, 35)
  Rl.out <- outLWradiation(LAI = LAI, Ts=Ts)
  Rl.inc <- incLWradiation(WeatherStation,DEM = surface.model$DEM, 
                           solar.angles = solar.angles.r, Ts= Ts)
  Rn <- netRadiation(LAI, albedo, Rs.inc, Rl.inc, Rl.out)
  plot(Rn, main="METRIC Net Radiation")
  return(Rn)
}

#' Estimates Net Radiation as in METRIC Model
#' @param image.DN         raw imagen in digital counts to evaluate
#' @param WeatherStation   Weather Station data, can be a waterWeatherStation 
#' object
#' @param Rn               RasterLayer with Net Radiation data in W/m2
#' @param plain            Logical. If TRUE surface is assumed plain
#' @param DEM              Digital Elevation Model of the study area. Not needed
#' if plain = TRUE
#' @param aoi              SpatialPolygon object with limits of Area of interest
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @family METRIC model functions
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
METRIC.G <- function(image.DN, WeatherStation=WeatherStation, Rn,  
                     plain = TRUE, DEM, aoi){
  path=getwd()
  if(plain==TRUE){
    DEM <- raster(image.DN[[1]])
    values(DEM) <- WeatherStation$location$elev
  } 
  surface.model <-METRICtopo(DEM)
  solar.angles.r <- solarAngles(surface.model = surface.model)
  image.TOAr <- calcTOAr(image.DN = image.DN, 
                         incidence.rel = solar.angles.r$incidence.rel)
  image.SR <- calcSR(image.TOAr=image.TOAr, 
                      surface.model=surface.model, 
                      incidence.hor = solar.angles.r$incidence.hor, 
                      WeatherStation=WeatherStation, sat="auto")
  albedo <- albedo(image.SR = image.SR)
  Ts <- surfaceTemperature(sat = "auto" )
  G <- soilHeatFlux(image = image.SR, Ts=Ts,albedo=albedo, Rn)
}


#' Estimates Energy Balance using METRIC2010 Model
#' @param image.DN         raw imagen in digital counts to evaluate
#' @param image.SR         L8 ONLY. Surface reflectance imagen. water package does not 
#' include a model to calculate surface reflectance for Landsat 8 images. Landsat 8 users 
#' should download precalculated surface reflectances from espa website 
#' (espa.cr.usgs.gov). 
#' @param WeatherStation   Weather Station data, can be a waterWeatherStation 
#' object
#' @param MTL              Landsat metadata file
#' @param sat              Landsat satellite version. "L7" or "L8"
#' @param thermalband      Landsat low gain thermalband
#' @param plain            Logical. If TRUE surface is assumed plain
#' @param DEM              Digital Elevation Model of the study area. Not needed
#' if plain = TRUE
#' @param aoi              SpatialPolygon object with limits of Area of interest
#' @param G.method         method used for the G estimation. Currently implemeted are 
#'                         "Tasumi" for Tasumi,2003 or "Bastiaanssen" for Bastiaanssen, 2000

#' @param alb.coeff        coefficient to transform narrow to broad band albedo.
#' See Details.
#' @param LST.method       Method for land surface temperature estimation. "SC" 
#' for single channel or "SW" for split window algorithm. 
#' "SW" is only available for L8. See \code{water::surfaceTemperature}
#' @param LAI.method       Method used to estimate LAI from spectral data. 
#' See Details.
#' @param L                L value for SAVI calculation
#' @param Zom.method       method selected to calculate momentum roughness 
#' length. Use "short.crops" for short crops methods from Allen et al (2007); 
#' "custom" for custom method also in Allen et al (2007); Or "Perrier" to use 
#' Perrier equation as in Santos et al (2012) and Pocas et al (2014).
#' @param anchors.method   method for the automatic selection of the anchor pixels. 
#' @param anchors          data.frame or SpatialPointsDataFrame with the anchor
#'                         pixels. The data frame must include a "type" column 
#'                         with "hot" and "cold" values.
#' @param n                number of pair of anchors pixels to calculate
#' @param ETp.coef         ETp coefficient usually 1.05 or 1.2 for alfalfa
#' @param Z.om.ws          momentum roughness lenght for WeatherStation. Usually
#' 0.0018 or 0.03 for long grass
#' @param verbose          Logical. If TRUE will print aditional data to console
#' @param extraParameters  Extra parameters for the non default methods. i.e. 
#' Zom.method = "Perrier", needs two extra parameters: fLAI, h. See 
#' help(momentumRoughnessLength).
#' @details
#' There are differents models to convert narrowband data to broadband albedo. 
#' You can choose alb.coeff ="Tasumi" to use Tasumi et al (2008) coefficients, 
#' calculated for Landsat 7; alb.coeff ="Liang" to use Liang Landsat 7 
#' coefficients or "Olmedo" to use Olmedo coefficients for Landsat 8.
#' @section Extra Parameters:
#' Extra Paramenters for functions inside METRIC.EB() include:
#' * for momentumRoughness when Zom.method = "Perrier": fLAI, h.
#' * for calcAnchors(): minDist, WSbuffer, deltaTemp
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @family METRIC model functions
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' González, Arturo & Hay, Christopher & Kjaersgaard, Jeppe & Neale, Christopher. (2015). Use of Remote Sensing to Generate Crop Coefficient and Estimate Actual Crop Evapotranspiration. 10.13031/aim.20152190105. 
#' @examples 
#' ### Data preparation
#' library(water)
#' aoi <- createAoi(topleft = c(500000, -3644000), bottomright = c(526000, -3660000))
#' raw_data_folder <- system.file("extdata", package="water")
#' image <- loadImage(path=raw_data_folder, aoi=aoi, sat="L8")
#' image.SR <- loadImageSR(path=raw_data_folder, aoi=aoi)
#' csvfile <- system.file("extdata", "INTA.csv", package="water")
#' MTLfile <- system.file("extdata", "LC82320832016040LGN00_MTL.txt", package="water")
#' \dontrun{
#' WeatherStation <- read.WSdata(WSdata = csvfile, 
#'                               datetime.format =  "%Y/%m/%d %H:%M", 
#'                               columns = c("datetime", "temp",
#'                                           "RH", "pp", "radiation", "wind"), 
#'                               lat=-33.00513, long= -68.86469, elev=927, height= 2,
#'                               MTL=MTLfile)
#'                               
#' ### LSEB with default methods and no extra parameters                               
#' Energy.Balance <- METRIC.EB(image.DN = image, image.SR = image.SR,
#'                             plain=TRUE, aoi=aoi, n = 5, WeatherStation = WeatherStation, 
#'                             ETp.coef = 1.2, sat="L8", alb.coeff = "Olmedo", LST.method = "SW", 
#'                             LAI.method = "metric2010", Z.om.ws = 0.03, MTL = MTLfile)
#'                             
#'                             
#' ### LSEB with "Perrier" method for Zom and extra parameters                               
#' Energy.Balance <- METRIC.EB(image.DN = image, image.SR = image.SR,
#'                             plain=TRUE, aoi=aoi, n = 5, WeatherStation = WeatherStation, 
#'                             ETp.coef = 1.2, sat="L8", alb.coeff = "Olmedo", LST.method = "SW", 
#'                             LAI.method = "metric2010", Zom.method = "Perrier", Z.om.ws = 0.03, 
#'                             MTL = MTLfile, extraParameters = c(fLAI = 0.5, h = 1.8) ) 
#' }
#' @export
METRIC.EB <- function(image.DN, image.SR, WeatherStation, MTL, sat = "auto",
                      thermalband, plain=TRUE, DEM, aoi, G.method = "Tasumi",
                      alb.coeff = "Tasumi", LST.method = "SC",
                      LAI.method = "metric2010", L = 0.1,
                      Zom.method = "short.crops", anchors.method = "CITRA-MCB",
                      anchors, n = 1, ETp.coef= 1.05, Z.om.ws=0.0018,
                      verbose = FALSE, extraParameters = vector()){
  path=getwd()
  #pb <- txtProgressBar(min = 0, max = 100, style = 3)
  if(plain==TRUE){
    DEM <- raster(image.DN[[1]])
    values(DEM) <- WeatherStation$location$elev
  }
  surface.model <-METRICtopo(DEM)
  #setTxtProgressBar(pb, 3)
  if(missing(MTL)){MTL <- list.files(path = path, pattern = "MTL.txt", full.names = T)}
  solar.angles.r <- solarAngles(surface.model = surface.model, 
                                WeatherStation = WeatherStation, MTL = MTL)
  Rs.inc <- incSWradiation(surface.model = surface.model, 
                           solar.angles = solar.angles.r, 
                           WeatherStation = WeatherStation)
  if(sat=="L7" | sat=="L8"){
    image.TOAr <- calcTOAr(image.DN = image.DN, sat=sat, MTL = MTL,
                           incidence.rel = solar.angles.r$incidence.rel)
    if(sat=="L7"){
      image.SR <- calcSR(image.TOAr=image.TOAr, sat = sat, 
                         surface.model=surface.model, 
                         incidence.hor = solar.angles.r$incidence.hor, 
                         WeatherStation=WeatherStation)}
    }
  if(sat=="MODIS"){image.SR <- image.DN}
  albedo <- albedo(image.SR = image.SR,  coeff=alb.coeff, sat=sat)
  #setTxtProgressBar(pb, 6)
  if(sat=="MODIS"){image.TOAr <- image.DN} # Only used for LAI estimation,
  # and some LAI models, use SR
  if(LAI.method == "vineyard"){LAI <- LAI(method = LAI.method, image = image.DN, L=L)}
  if(LAI.method == "turner"){LAI <- LAI(method = LAI.method, image = image.SR, L=L)}
  if (LAI.method != "vineyard" & LAI.method != "turner"){LAI <- LAI(method = LAI.method, image = image.TOAr, L=L)}
  if(sat=="L7" | sat=="L8"){
    Ts <- surfaceTemperature(LAI=LAI, sat = sat, image.DN=image.DN,
                             WeatherStation = WeatherStation, method = LST.method)}
  if(sat=="MODIS"){Ts <- image.DN$LST}
  #setTxtProgressBar(pb, 35)
  Rl.out <- outLWradiation(LAI = LAI, Ts=Ts)
  Rl.inc <- incLWradiation(WeatherStation,DEM = surface.model$DEM, 
                           solar.angles = solar.angles.r, Ts= Ts)
  Rn <- netRadiation(LAI, albedo, Rs.inc, Rl.inc, Rl.out)
  # Rn[Rn < 0]  <-  0 # see H
  #setTxtProgressBar(pb, 40)
  G <- soilHeatFlux(image = image.SR, Ts=Ts,albedo=albedo, 
                    Rn=Rn, LAI=LAI, method = G.method)
  # G[G < 0]  <-  0 # see H
  if(Zom.method == "short.crops"){
    Z.om <- momentumRoughnessLength(LAI=LAI, mountainous = !plain, 
                                    method = Zom.method, 
                                    surface.model = surface.model)
  }
  if(Zom.method == "custom"){
    Z.om <- momentumRoughnessLength(LAI=LAI, mountainous = !plain, 
                                    method = Zom.method, a = extraParameters["a"],
                                    b = extraParameters["b"],
                                    surface.model = surface.model)
  }
  if(Zom.method == "Perrier"){
    Z.om <- momentumRoughnessLength(LAI=LAI, mountainous = !plain, 
                                    method = Zom.method, fLAI = extraParameters["fLAI"],
                                    h = extraParameters["h"],
                                    surface.model = surface.model)
  }
  if(Zom.method != "short.crops"){Z.om.sc <- momentumRoughnessLength(LAI=LAI, mountainous = !plain, 
                                                             method = "short.crops", 
                                                             surface.model = surface.model)
  } else {Z.om.sc <- Z.om}
  par(mfrow=c(1,2))
  if(missing(anchors)){
    if(is.na(extraParameters["deltaTemp"])){extraParameters['deltaTemp'] = 5}
    if(is.na(extraParameters["minDist"])){extraParameters['minDist'] = 500}
    if(is.na(extraParameters["WSbuffer"])){extraParameters['WSbuffer'] = 30000}
    anchors <- calcAnchors(image = image.TOAr, Ts = Ts, LAI = LAI, plots = T,
                           albedo = albedo, Z.om = Z.om.sc, n = n, 
                           anchors.method = anchors.method, WeatherStation = WeatherStation,
                           deltaTemp = extraParameters['deltaTemp'],
                           minDist = extraParameters['minDist'],
                           WSbuffer = extraParameters['WSbuffer'],
                           verbose = verbose)
    print(anchors)
  }
  #setTxtProgressBar(pb, 45)
  if(anchors.method == "flexible"){
     flex.cold <- anchors$flex
     NDVIcold <- mean(as.numeric(anchors$NDVI[anchors$type == "cold"]))
     print(paste0("mean NDVI for cold pixel: ", NDVIcold))
     # NDVI and Kc relationship from:
     #González, Arturo & Hay, Christopher & Kjaersgaard, Jeppe & Neale, 
     #Christopher. (2015). Use of Remote Sensing to Generate Crop Coefficient 
     #and Estimate Actual Crop Evapotranspiration. 10.13031/aim.20152190105. 
     ETp.coef <-  2.11 * NDVIcold - 0.4989 
     print(paste0("ETp.coef updated to: ", ETp.coef, " (Gonzalez et al, 2015)"))

  }
  H <- calcH(anchors = anchors, Ts = Ts, Z.om = Z.om, mountainous = !plain,
             WeatherStation = WeatherStation, ETp.coef = ETp.coef,
             Z.om.ws = Z.om.ws, DEM = DEM, Rn = Rn, G = G, verbose = verbose)
  par(mfrow=c(1,1))
  #setTxtProgressBar(pb, 99)
  H <-  H$H
  # H[H < 0]  <-  0  not a good solution...! not for any of the components
  LE <- Rn - G - H
  LE[LE < 0]  <-  0   # see H
  EB <- stack(Rn, G, H, LE, Ts)
  EB <- saveLoadClean(imagestack = EB,
                stack.names = c("NetRadiation", "SoilHeat", "SensibleHeat", 
                                "LatentHeat", "surfaceTemperature"), 
                file = "EB", overwrite=TRUE)
  #setTxtProgressBar(pb, 100)
  result <- list()
  result$EB <- EB
  result$WeatherStation <- WeatherStation
  coordinates(anchors) <- ~ X + Y
  result$anchors <- anchors
  result$methods <- c(sat = sat, alb.coeff = alb.coeff, LST.method = LST.method,
                      LAI.method = LAI.method, Zom.method = Zom.method, 
                      anchors.method = anchors.method, plain = plain,
                      ETp.coef= ETp.coef, Z.om.ws=Z.om.ws,
                      extraParameters = extraParameters)
  class(result) <- "waterLSEB"
  return(result)
}
  
