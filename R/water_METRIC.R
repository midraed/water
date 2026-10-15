#' Estimates Net Radiation as in METRIC Model
#' @description
#' Estimates the net radiation (Rn) as in the METRIC model, with the same
#' inputs and steps as METRIC.EB().
#' For Landsat 8 and 9, use Collection 2 Level-2 products: surface reflectance
#' as image.SR and surface temperature as surf.temp (both from loadImageSR()).
#' For Landsat 7, image.SR and Ts can also be computed from image.DN.
#' @param image.DN         Raw image in digital counts. Optional if image.SR and
#' surf.temp are supplied. If supplied (L7/L8/L9) it is used for TOA
#' reflectance (LAI) and, without surf.temp, for surface temperature.
#' @param WeatherStation   waterWeatherStation object with the conditions at
#' the satellite overpass (from read.WSdata() with the MTL file)
#' @param MTL              path to the Landsat metadata (MTL) file. If missing,
#' the only *MTL.txt file in the working directory is used
#' @param sat              Satellite: "L7", "L8", "L9", "MODIS" or "auto" to
#' detect it from the MTL file
#' @param thermalband      Thermal band in digital numbers (L7: band 6 low gain;
#' L8/L9: band 10). If missing, image.DN$Thermal1 is used. Not needed if
#' surf.temp is supplied
#' @param plain            Logical. If TRUE (default) the surface is assumed
#' flat: DEM is replaced by a constant raster at the weather station elevation
#' @param DEM              SpatRaster with a digital elevation model (m) of the
#' study area, on the same grid as the image. Not needed if plain = TRUE
#' @param aoi              Optional sf polygon or SpatVector. If given, the
#' inputs are cropped to it
#' @param alb.coeff        coefficient to transform narrow to broad band albedo.
#' See Details.
#' @param LAI.method       Method used to estimate LAI from spectral data.
#' See Details.
#' @param image.SR         Surface reflectance image (B, G, R, NIR, SWIR1,
#' SWIR2), e.g. from loadImageSR(). Required for L8 and L9; for L7 it is
#' computed from image.DN if missing.
#' @param surf.temp        Optional. Surface temperature raster in Kelvin, on the
#' same grid as image.SR (e.g. image.SR$Band_10 from loadImageSR()). If
#' supplied it is used as Ts and surfaceTemperature() is not called.
#' @param L                L value for SAVI calculation. See LAI()
#' @param LST.method       Method for land surface temperature estimation: "SC"
#' for single channel or "SW" for split window (only for L8 and L9). Ignored
#' if surf.temp is supplied. See surfaceTemperature()
#' @details
#' There are different models to convert narrowband data to broadband albedo.
#' You can choose alb.coeff = "Tasumi" to use Tasumi et al. (2008)
#' coefficients, calculated for Landsat 7; alb.coeff = "Liang" to use Liang
#' (2001) Landsat 7 coefficients or alb.coeff = "Olmedo" to use Olmedo
#' coefficients for Landsat 8 and 9. See \code{\link{albedo}}.
#'
#' LAI.method can be "metric2010" (default), "metric", "MCB", "turner" or
#' "vineyard". See \code{\link{LAI}}. Without image.DN, the leaf area index is
#' computed from the surface reflectance (image.SR), although the "metric",
#' "metric2010" and "MCB" equations were described for top-of-atmosphere
#' reflectance.
#' @return SpatRaster with net radiation (W/m2). The result is also plotted.
#' @family METRIC model functions
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
METRIC.Rn <- function(image.DN, WeatherStation, MTL, sat = "auto", thermalband,
                      alb.coeff = "Tasumi", LAI.method = "metric2010",
                      plain = TRUE, DEM, aoi, image.SR, surf.temp, L = 0.1,
                      LST.method = "SC"){
  path=getwd()
  if(sat=="auto"){sat <- getSat(path, MTL = if(!missing(MTL)) MTL)}
  MTL <- findMTL(if(!missing(MTL)) MTL, path = path)
  surface <- metricSurface(image.DN = image.DN, image.SR = image.SR,
                           WeatherStation = WeatherStation, MTL = MTL, sat = sat,
                           thermalband = thermalband, surf.temp = surf.temp,
                           alb.coeff = alb.coeff, LAI.method = LAI.method, L = L,
                           LST.method = LST.method, plain = plain, DEM = DEM,
                           aoi = aoi)
  Rn <- metricRn(surface, WeatherStation)
  plot(Rn, main="METRIC Net Radiation")
  return(Rn)
}

#' Estimates Soil Heat Flux as in METRIC Model
#' @description
#' Estimates the soil heat flux (G) as in the METRIC model.
#' Inputs are the same as in METRIC.Rn(). For Landsat 8 and 9, use Collection 2
#' Level-2 surface reflectance as image.SR and surface temperature as surf.temp.
#' @inheritParams METRIC.Rn
#' @inherit METRIC.Rn details
#' @param Rn               Net radiation in W/m2, e.g. from METRIC.Rn(). If
#' missing it is computed as in METRIC.Rn()
#' @param G.method         method used for the G estimation: "Tasumi" (Tasumi,
#' 2003; default) or "Bastiaanssen" (Bastiaanssen, 2000). See soilHeatFlux()
#' @return SpatRaster with soil heat flux G (W/m2). The result is also plotted.
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @family METRIC model functions
#' @references
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
METRIC.G <- function(image.DN, WeatherStation, MTL, sat = "auto", thermalband,
                     Rn, alb.coeff = "Tasumi", LAI.method = "metric2010",
                     G.method = "Tasumi", plain = TRUE, DEM, aoi, image.SR,
                     surf.temp, L = 0.1, LST.method = "SC"){
  path=getwd()
  if(sat=="auto"){sat <- getSat(path, MTL = if(!missing(MTL)) MTL)}
  MTL <- findMTL(if(!missing(MTL)) MTL, path = path)
  surface <- metricSurface(image.DN = image.DN, image.SR = image.SR,
                           WeatherStation = WeatherStation, MTL = MTL, sat = sat,
                           thermalband = thermalband, surf.temp = surf.temp,
                           alb.coeff = alb.coeff, LAI.method = LAI.method, L = L,
                           LST.method = LST.method, plain = plain, DEM = DEM,
                           aoi = aoi)
  if(missing(Rn)){Rn <- metricRn(surface, WeatherStation)}
  G <- soilHeatFlux(image = surface$image.SR, Ts = surface$Ts,
                    albedo = surface$albedo, Rn = Rn, LAI = surface$LAI,
                    method = G.method)
  plot(G, main="METRIC Soil Heat Flux")
  return(G)
}

# Surface inputs shared by METRIC.EB(), METRIC.Rn() and METRIC.G(): supplied
# image.SR > SR derived from image.DN (L7, MODIS), and supplied surf.temp >
# surfaceTemperature()
metricSurface <- function(image.DN, image.SR, WeatherStation, MTL, sat,
                          thermalband, surf.temp, alb.coeff, LAI.method,
                          L = 0.1, LST.method = "SC", plain, DEM, aoi){
  if(missing(image.SR) && (missing(image.DN) || !sat %in% c("L7", "MODIS"))){
    stop("image.SR is missing. It can only be computed from image.DN for L7 or ",
         "MODIS; for L8/L9 provide it (e.g., loadImageSR()).")
  }
  if(missing(surf.temp) && missing(image.DN)){
    stop("Provide surf.temp (e.g., image.SR$Band_10 from loadImageSR()) or ",
         "image.DN with a thermal band.")
  }
  if(!missing(aoi)){
    if(!missing(image.DN)){image.DN <- crop(image.DN, aoi)}
    if(!missing(image.SR)){image.SR <- crop(image.SR, aoi)}
    if(!missing(surf.temp)){surf.temp <- crop(surf.temp, aoi)}
    if(!missing(thermalband)){thermalband <- crop(thermalband, aoi)}
    if(!plain){DEM <- crop(DEM, aoi)}
  }
  if(plain==TRUE){
    DEM <- rast(if(missing(image.SR)) image.DN[[1]] else image.SR[[1]])
    values(DEM) <- WeatherStation$location$elev
  }
  surface.model <- METRICtopo(DEM)
  solar.angles <- solarAngles(surface.model = surface.model,
                              WeatherStation = WeatherStation, MTL = MTL)
  if(!missing(image.DN) && sat %in% c("L7", "L8", "L9")){
    image.TOAr <- calcTOAr(image.DN = image.DN, sat = sat, MTL = MTL,
                           incidence.rel = solar.angles$incidence.rel)
  }
  if(missing(image.SR)){
    if(sat=="L7"){
      image.SR <- calcSR(image.TOAr = image.TOAr, sat = sat,
                         surface.model = surface.model,
                         incidence.hor = solar.angles$incidence.hor,
                         WeatherStation = WeatherStation)
    } else {
      image.SR <- image.DN  # MODIS
    }
  }
  # Without DN (or for MODIS) LAI uses SR instead of TOAr
  if(missing(image.DN) || sat=="MODIS"){image.TOAr <- image.SR}
  # aoi = extent of image.SR, so albedo() does not crop to a global aoi object
  albedo <- albedo(image.SR = image.SR, aoi = ext(image.SR), coeff = alb.coeff,
                   sat = sat)
  if(LAI.method == "vineyard"){
    if(missing(image.DN)){stop("LAI.method = 'vineyard' needs image.DN.")}
    LAI <- LAI(method = LAI.method, L = L,
               image = calcRadiance(image.DN, sat = sat, MTL = MTL))
  } else if(LAI.method == "turner"){
    LAI <- LAI(method = LAI.method, image = image.SR, L = L)
  } else {
    LAI <- LAI(method = LAI.method, image = image.TOAr, L = L)
  }
  if(!missing(surf.temp)){
    Ts <- surf.temp
  } else if(sat %in% c("L7", "L8", "L9")){
    Ts <- surfaceTemperature(LAI = LAI, sat = sat, image.DN = image.DN,
                             method = LST.method, thermalband = thermalband,
                             MTL = MTL, WeatherStation = WeatherStation)
  } else {
    Ts <- image.DN$LST  # MODIS
  }
  list(DEM = DEM, surface.model = surface.model, solar.angles = solar.angles,
       image.TOAr = image.TOAr, image.SR = image.SR, albedo = albedo,
       LAI = LAI, Ts = Ts)
}

# Net radiation from the output of metricSurface()
metricRn <- function(surface, WeatherStation){
  Rs.inc <- incSWradiation(surface.model = surface$surface.model,
                           solar.angles = surface$solar.angles,
                           WeatherStation = WeatherStation)
  Rl.out <- outLWradiation(LAI = surface$LAI, Ts = surface$Ts)
  Rl.inc <- incLWradiation(WeatherStation, DEM = surface$surface.model$DEM,
                           solar.angles = surface$solar.angles, Ts = surface$Ts)
  netRadiation(surface$LAI, surface$albedo, Rs.inc, Rl.inc, Rl.out)
}


#' Estimates Energy Balance using METRIC2010 Model
#' @description
#' Estimates the instantaneous land surface energy balance (Rn, G, H and LE)
#' with the METRIC model (Allen et al., 2007). For Landsat 8 and 9 use
#' Collection 2 Level-2 surface reflectance (image.SR) and, optionally, surface
#' temperature (surf.temp), both from loadImageSR().
#' @param image.DN         Raw image in digital counts. Optional if image.SR is
#' supplied. If supplied (L7/L8/L9) it is used for TOA reflectance (LAI and
#' anchors) and for surface temperature; if missing, image.SR is used instead.
#' Needed for LAI.method = "vineyard" and, without surf.temp, for MODIS (LST layer).
#' @param image.SR         Surface reflectance image (B, G, R, NIR, SWIR1,
#' SWIR2), e.g. from loadImageSR(). If supplied it is used directly; if
#' missing it is computed from image.DN (L7 via calcSR(), MODIS uses image.DN).
#' Required for L8 and L9.
#' @param WeatherStation   waterWeatherStation object with the conditions at
#' the satellite overpass (from read.WSdata() with the MTL file)
#' @param MTL              path to the Landsat metadata (MTL) file. If missing,
#' the only *MTL.txt file in the working directory is used
#' @param sat              Satellite: "L7", "L8", "L9", "MODIS" or "auto" to
#' detect it from the MTL file
#' @param thermalband      Thermal band in digital numbers (L7: band 6 low gain;
#' L8/L9: band 10). If missing, image.DN$Thermal1 is used. Not needed if
#' surf.temp is supplied
#' @param surf.temp        Optional. Surface temperature raster in Kelvin, on the
#' same grid as image.SR (e.g. image.SR$Band_10 from loadImageSR()). If
#' supplied it is used as Ts and surfaceTemperature() is not called.
#' @param plain            Logical. If TRUE (default) the surface is assumed
#' flat: DEM is replaced by a constant raster at the weather station elevation
#' @param DEM              SpatRaster with a digital elevation model (m) of the
#' study area, on the same grid as the image. Not needed if plain = TRUE
#' @param aoi              Optional sf polygon or SpatVector. If given, the
#' inputs are cropped to it
#' @param G.method         method used for the G estimation: "Tasumi" (Tasumi,
#' 2003; default) or "Bastiaanssen" (Bastiaanssen, 2000). See soilHeatFlux()
#' @param alb.coeff        coefficient to transform narrow to broad band albedo.
#' See Details.
#' @param LST.method       Method for land surface temperature estimation: "SC"
#' for single channel or "SW" for split window (only for L8 and L9). Ignored
#' if surf.temp is supplied. See \code{\link{surfaceTemperature}}
#' @param LAI.method       Method used to estimate LAI from spectral data.
#' See Details.
#' @param L                L value for SAVI calculation. See LAI()
#' @param Zom.method       method selected to calculate momentum roughness
#' length. Use "short.crops" for short crops methods from Allen et al (2007);
#' "custom" for custom method also in Allen et al (2007); Or "Perrier" to use
#' Perrier equation as in Santos et al (2012) and Pocas et al (2014).
#' @param anchors.method   method for the automatic selection of the anchor
#' pixels: "flexible" (default), "best" or "random". See calcAnchors()
#' @param anchors          data.frame with X, Y (coordinates in the CRS of the
#' image) and type ("hot" or "cold") columns, e.g. from calcAnchors(), or an sf
#' POINT object with a type column. If missing, the anchors are selected
#' automatically with calcAnchors()
#' @param n                number of pairs of anchor pixels to select
#' @param ETp.coef         coefficient for the ET of the cold pixels: calcH()
#' multiplies it by the hourly grass reference ET (ETo) at the satellite
#' overpass. Usually 1.05, or 1.2 for alfalfa. It can be updated with
#' anchors.method = "flexible" (see Details)
#' @param Z.om.ws          momentum roughness length (m) of the weather station
#' surface. The default is 0.0018 here and 0.03 in calcH(); see calcH() for
#' typical values
#' @param verbose          Logical. If TRUE will print additional data to console
#' @param extraParameters  Extra parameters for the non default methods, e.g.
#' Zom.method = "Perrier" needs two extra parameters: fLAI, h. See
#' \code{\link{momentumRoughnessLength}} and the Extra Parameters section.
#' @details
#' There are different models to convert narrowband data to broadband albedo.
#' You can choose alb.coeff = "Tasumi" to use Tasumi et al. (2008)
#' coefficients, calculated for Landsat 7; alb.coeff = "Liang" to use Liang
#' (2001) Landsat 7 coefficients or alb.coeff = "Olmedo" to use Olmedo
#' coefficients for Landsat 8 and 9. See \code{\link{albedo}}.
#'
#' LAI.method can be "metric2010" (default), "metric", "MCB", "turner" or
#' "vineyard". See \code{\link{LAI}}. Without image.DN, the leaf area index is
#' computed and the anchor pixels are selected from the surface reflectance
#' (image.SR), although the "metric", "metric2010" and "MCB" LAI equations were
#' described for top-of-atmosphere reflectance.
#'
#' The ET of the cold pixels is ETp.coef times the hourly grass reference ET
#' (ETo) at the satellite overpass (see calcH()), while ET24h() and dailyET()
#' use the alfalfa reference ET (ETr) by default. With anchors.method =
#' "flexible", if the cold-pixel criteria had to be relaxed and the mean NDVI of
#' the cold pixels is below 0.8, ETp.coef is replaced by 2.11 * NDVI - 0.4989
#' (Gonzalez et al., 2015).
#' @section Extra Parameters:
#' Extra parameters for functions inside METRIC.EB() include:
#' \itemize{
#'   \item for momentumRoughnessLength() when Zom.method = "Perrier": fLAI, h;
#'   \item for momentumRoughnessLength() when Zom.method = "custom": a, b;
#'   \item for calcAnchors(): minDist (default 500), WSbuffer (default 30000)
#'     and deltaTemp (default 5).
#' }
#' @return An object of class \code{waterLSEB}, a list with:
#' \describe{
#'   \item{EB}{SpatRaster with layers NetRadiation, SoilHeat, SensibleHeat and
#'     LatentHeat (W/m2), and surfaceTemperature (K).}
#'   \item{WeatherStation}{the weather station data used.}
#'   \item{anchors}{sf POINT object with the hot and cold anchor pixels.}
#'   \item{methods}{named character vector with the methods and parameters used.}
#' }
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @family METRIC model functions
#' @references
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#'
#' González, Arturo & Hay, Christopher & Kjaersgaard, Jeppe & Neale, Christopher. (2015). Use of Remote Sensing to Generate Crop Coefficient and Estimate Actual Crop Evapotranspiration. 10.13031/aim.20152190105. 
#' @examples
#' ### Data preparation: Landsat 8 Collection 2 Level-2 example crop
#' library(water)
#' scene_folder <- system.file("extdata", "L8_C2", package = "water")
#' image.SR <- loadImageSR(path = scene_folder)
#' image.SR <- cfmask(path = scene_folder, image = image.SR)
#' MTLfile <- list.files(scene_folder, pattern = "MTL.txt$", full.names = TRUE)
#' csvfile <- system.file("extdata", "INTA.csv", package = "water")
#' WeatherStation <- read.WSdata(WSdata = csvfile, datetime.format = "%Y/%m/%d %H:%M",
#'                               columns = c("date" = 1, "time" = 1, "radiation" = 5,
#'                                           "wind" = 6, "RH" = 3, "temp" = 2, "rain" = 4),
#'                               lat = -33.00513, long = -68.86469, elev = 927, height = 2,
#'                               tz = "America/Argentina/Mendoza", MTL = MTLfile)
#' \donttest{
#' ### LSEB with surface reflectance and surface temperature
#' Energy.Balance <- METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10,
#'                             WeatherStation = WeatherStation, MTL = MTLfile,
#'                             alb.coeff = "Olmedo", n = 5, ETp.coef = 1.2,
#'                             Z.om.ws = 0.03)
#'
#' ### LSEB with "Perrier" method for Zom and extra parameters
#' Energy.Balance <- METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10,
#'                             WeatherStation = WeatherStation, MTL = MTLfile,
#'                             alb.coeff = "Olmedo", n = 5, ETp.coef = 1.2,
#'                             Z.om.ws = 0.03, Zom.method = "Perrier",
#'                             extraParameters = c(fLAI = 0.5, h = 1.8))
#' }
#' @export
METRIC.EB <- function(image.DN, image.SR, WeatherStation, MTL, sat = "auto",
                      thermalband, plain=TRUE, DEM, aoi, G.method = "Tasumi",
                      alb.coeff = "Tasumi", LST.method = "SC",
                      LAI.method = "metric2010", L = 0.1,
                      Zom.method = "short.crops", anchors.method = "flexible",
                      anchors, n = 1, ETp.coef= 1.05, Z.om.ws=0.0018,
                      verbose = FALSE, extraParameters = vector(), surf.temp){
  path=getwd()
  if(sat=="auto"){sat <- getSat(path, MTL = if(!missing(MTL)) MTL)}
  if(missing(image.SR) && (missing(image.DN) || !sat %in% c("L7", "MODIS"))){
    stop("image.SR is missing. It can only be computed from image.DN for L7 or ",
         "MODIS; provide it (e.g., loadImageSR()).")
  }

  if(missing(surf.temp) && missing(thermalband) && sat %in% c("L7", "L8", "L9")){
    if(missing(image.DN)){
      stop("Provide surf.temp, thermalband or image.DN with a thermal layer.")
    }
    thermal.layer <- intersect(c("Thermal1", "thermal.low"), names(image.DN))
    if(length(thermal.layer) == 0){
      stop("thermalband is missing. Provide a thermal layer (e.g., image.DN$Thermal1).")
    }
    thermalband <- image.DN[[thermal.layer[1]]]
  }

  if(missing(surf.temp) && sat %in% c("L7", "L8", "L9") && is.null(thermalband)){
    stop("thermalband is NULL. Check the band name in image.DN (e.g., names(image.DN)).")
  }

  MTL <- findMTL(if(!missing(MTL)) MTL, path = path)
  surface <- metricSurface(image.DN = image.DN, image.SR = image.SR,
                           WeatherStation = WeatherStation, MTL = MTL, sat = sat,
                           thermalband = thermalband, surf.temp = surf.temp,
                           alb.coeff = alb.coeff, LAI.method = LAI.method, L = L,
                           LST.method = LST.method, plain = plain, DEM = DEM,
                           aoi = aoi)
  DEM <- surface$DEM
  surface.model <- surface$surface.model
  image.TOAr <- surface$image.TOAr
  image.SR <- surface$image.SR
  albedo <- surface$albedo
  LAI <- surface$LAI
  Ts <- surface$Ts
  Rn <- metricRn(surface, WeatherStation)
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
    NDVI <- (image.TOAr$NIR - image.TOAr$R) / (image.TOAr$NIR + image.TOAr$R)
    Z.om <- momentumRoughnessLength(LAI=LAI, NDVI = NDVI, albedo = albedo,
                                    mountainous = !plain,
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
  if(anchors.method == "flexible" && !is.null(anchors$flex)){
     flex.cold <- anchors$flex[anchors$type == "cold"][1]
     NDVIcold <- mean(as.numeric(anchors$NDVI[anchors$type == "cold"]))
     if(flex.cold > 0 & NDVIcold < 0.8){
       print(paste0("mean NDVI for cold pixel: ", NDVIcold))
       # NDVI and Kc relationship from:
       #González, Arturo & Hay, Christopher & Kjaersgaard, Jeppe & Neale, 
       #Christopher. (2015). Use of Remote Sensing to Generate Crop Coefficient 
       #and Estimate Actual Crop Evapotranspiration. 10.13031/aim.20152190105. 
       ETp.coef <-  2.11 * NDVIcold - 0.4989 
       print(paste0("ETp.coef updated to: ", ETp.coef, " (Gonzalez et al, 2015)"))
     }
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
  EB <- c(Rn, G, H, LE, Ts)
  EB <- saveLoadClean(imagestack = EB,
                stack.names = c("NetRadiation", "SoilHeat", "SensibleHeat", 
                                "LatentHeat", "surfaceTemperature"), 
                file = "EB", overwrite=TRUE)
  #setTxtProgressBar(pb, 100)
  result <- list()
  result$EB <- EB
  result$WeatherStation <- WeatherStation
  if(!inherits(anchors, "sf")){
    anchors <- sf::st_as_sf(anchors, coords = c("X", "Y"), crs = terra::crs(Ts))
  }
  result$anchors <- anchors
  result$methods <- c(sat = sat, alb.coeff = alb.coeff, LST.method = LST.method,
                      LAI.method = LAI.method, Zom.method = Zom.method, 
                      anchors.method = anchors.method, plain = plain,
                      ETp.coef= ETp.coef, Z.om.ws=Z.om.ws,
                      extraParameters = extraParameters)
  class(result) <- "waterLSEB"
  return(result)
}
  
