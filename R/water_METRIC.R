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
#' @param MTL              path to the Landsat metadata (*_MTL.txt) file.
#' Required
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
#' @param image.SR         Surface reflectance image with layers B, G, R, NIR,
#' SWIR1 and SWIR2, e.g. from loadImageSR(). Required for L8 and L9; for L7
#' it is computed from image.DN if missing.
#' Layers are selected by name; a stack without the standard names is used by
#' position, with a warning, and a stack with these names that lacks a needed
#' layer is an error.
#' @param surf.temp        Optional. Surface temperature raster in Kelvin, on the
#' same grid as image.SR (e.g. image.SR$Band_10 from loadImageSR()). If
#' supplied it is used as Ts and surfaceTemperature() is not called.
#' Must be in Kelvin: values that look like Celsius, or more than 1\% of the
#' cells below 200 K (unmasked fill or clouds), stop; fewer cells below 200 K
#' give a warning.
#' @param L                L value for SAVI calculation. See \code{\link{LAI}}
#' @param LST.method       Land surface temperature when surf.temp is not
#' supplied: "auto" (default) uses the Level-2 surface temperature
#' (image.SR$Band_10 from loadImageSR()) for L8 and L9 when present, and the
#' METRIC equations (Allen et al., 2007, Eqs. 20-21) with image.DN otherwise;
#' "L2" requires Band_10; "METRIC" uses the thermal band of image.DN; "SW"
#' (split window, L8 and L9) is deprecated. See \code{\link{surfaceTemperature}}
#' @param Kt               clearness (turbidity) coefficient of the atmosphere,
#' 0 < Kt <= 1: 1 (default) for clean air, 0.5 for extremely turbid, dusty or
#' polluted air (Allen et al., 2007). See \code{\link{SWtrasmisivity}}
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
#'
#' Method arguments (sat, alb.coeff, LAI.method, LST.method and, where they
#' apply, G.method, Zom.method and anchors.method) must be one of the listed
#' names, exactly (case included): an unknown name stops before any
#' computation, with a hint for typos.
#' @return SpatRaster with net radiation (W/m2). The result is also plotted.
#' @family METRIC model functions
#' @author Guillermo Federico Olmedo
#' @references
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
#' @examples
#' scene_folder <- system.file("extdata", "L8_C2", package = "water")
#' image.SR <- loadImageSR(path = scene_folder)
#' MTLfile <- list.files(scene_folder, pattern = "MTL.txt$", full.names = TRUE)
#' csvfile <- system.file("extdata", "INTA.csv", package = "water")
#' WeatherStation <- read.WSdata(WSdata = csvfile, datetime.format = "%Y/%m/%d %H:%M",
#'                               columns = c("date" = 1, "time" = 1, "radiation" = 5,
#'                                           "wind" = 6, "RH" = 3, "temp" = 2, "rain" = 4),
#'                               lat = -33.00513, long = -68.86469, elev = 927, height = 2,
#'                               tz = "America/Argentina/Mendoza", MTL = MTLfile)
#' Rn <- METRIC.Rn(image.SR = image.SR, surf.temp = image.SR$Band_10,
#'                 WeatherStation = WeatherStation, MTL = MTLfile,
#'                 alb.coeff = "Olmedo")
METRIC.Rn <- function(image.DN, WeatherStation, MTL, sat = "auto", thermalband, 
                      alb.coeff = "Tasumi", LAI.method = "metric2010",
                      plain = TRUE, DEM, aoi, image.SR, surf.temp, L = 0.1,
                      LST.method = "auto", Kt = 1){
  checkMetricArgs(sat = sat, alb.coeff = alb.coeff, LAI.method = LAI.method,
                  LST.method = LST.method)
  checkKt(Kt)
  if(sat=="auto"){sat <- getSat(MTL = if(!missing(MTL)) MTL)}
  checkLSTmethod(LST.method, sat, missing(surf.temp))
  MTL <- findMTL(if(!missing(MTL)) MTL)
  surface <- metricSurface(image.DN = image.DN, image.SR = image.SR,
                           WeatherStation = WeatherStation, MTL = MTL, sat = sat,
                           thermalband = thermalband, surf.temp = surf.temp,
                           alb.coeff = alb.coeff, LAI.method = LAI.method, L = L,
                           LST.method = LST.method, plain = plain, DEM = DEM,
                           aoi = aoi, Kt = Kt)
  Rn <- metricRn(surface, WeatherStation, Kt = Kt)
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
#' 2003; default) or "Bastiaanssen" (Bastiaanssen, 2000). See \code{\link{soilHeatFlux}}
#' @return SpatRaster with soil heat flux G (W/m2). The result is also plotted.
#' @author Guillermo Federico Olmedo
#' @family METRIC model functions
#' @references
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
#' @examples
#' scene_folder <- system.file("extdata", "L8_C2", package = "water")
#' image.SR <- loadImageSR(path = scene_folder)
#' MTLfile <- list.files(scene_folder, pattern = "MTL.txt$", full.names = TRUE)
#' csvfile <- system.file("extdata", "INTA.csv", package = "water")
#' WeatherStation <- read.WSdata(WSdata = csvfile, datetime.format = "%Y/%m/%d %H:%M",
#'                               columns = c("date" = 1, "time" = 1, "radiation" = 5,
#'                                           "wind" = 6, "RH" = 3, "temp" = 2, "rain" = 4),
#'                               lat = -33.00513, long = -68.86469, elev = 927, height = 2,
#'                               tz = "America/Argentina/Mendoza", MTL = MTLfile)
#' G <- METRIC.G(image.SR = image.SR, surf.temp = image.SR$Band_10,
#'               WeatherStation = WeatherStation, MTL = MTLfile,
#'               alb.coeff = "Olmedo")
METRIC.G <- function(image.DN, WeatherStation, MTL, sat = "auto", thermalband,
                     Rn, alb.coeff = "Tasumi", LAI.method = "metric2010",
                     G.method = "Tasumi", plain = TRUE, DEM, aoi, image.SR,
                     surf.temp, L = 0.1, LST.method = "auto", Kt = 1){
  checkMetricArgs(sat = sat, alb.coeff = alb.coeff, LAI.method = LAI.method,
                  LST.method = LST.method, G.method = G.method)
  checkKt(Kt)
  if(sat=="auto"){sat <- getSat(MTL = if(!missing(MTL)) MTL)}
  checkLSTmethod(LST.method, sat, missing(surf.temp))
  MTL <- findMTL(if(!missing(MTL)) MTL)
  surface <- metricSurface(image.DN = image.DN, image.SR = image.SR,
                           WeatherStation = WeatherStation, MTL = MTL, sat = sat,
                           thermalband = thermalband, surf.temp = surf.temp,
                           alb.coeff = alb.coeff, LAI.method = LAI.method, L = L,
                           LST.method = LST.method, plain = plain, DEM = DEM,
                           aoi = aoi, Kt = Kt)
  if(missing(Rn)){Rn <- metricRn(surface, WeatherStation, Kt = Kt)}
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
                          L = 0.1, LST.method = "auto", plain, DEM, aoi, Kt = 1){
  if(missing(image.SR) && (missing(image.DN) || !sat %in% c("L7", "MODIS"))){
    stop("image.SR is missing. It can only be computed from image.DN for L7 or ",
         "MODIS; for L8/L9 provide it (e.g., loadImageSR()).")
  }
  lst <- resolveLST(LST.method, sat, if(!missing(image.SR)) image.SR,
                    !missing(surf.temp))
  if(lst %in% c("METRIC", "SW") && missing(image.DN) && missing(thermalband)){
    stop("Provide surf.temp (e.g., image.SR$Band_10 from loadImageSR()) or ",
         "image.DN with a thermal band.")
  }
  if(!missing(aoi)){
    if(!missing(image.DN)){image.DN <- cropToAoi(image.DN, aoi)}
    if(!missing(image.SR)){image.SR <- cropToAoi(image.SR, aoi)}
    if(!missing(surf.temp)){surf.temp <- cropToAoi(surf.temp, aoi)}
    if(!missing(thermalband)){thermalband <- cropToAoi(thermalband, aoi)}
    if(!plain){DEM <- cropToAoi(DEM, aoi)}
  }
  # thermal band from image.DN, as in METRIC.EB()
  if(lst %in% c("METRIC", "SW") && missing(thermalband) && !missing(image.DN) &&
     sat %in% c("L7", "L8", "L9")){
    thermal.layer <- intersect(c("Thermal1", "thermal.low"), names(image.DN))
    if(length(thermal.layer) == 0){
      stop("thermalband is missing. Provide a thermal layer (e.g., image.DN$Thermal1).")
    }
    thermalband <- image.DN[[thermal.layer[1]]]
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
                         WeatherStation = WeatherStation, Kt = Kt)
    } else {
      image.SR <- image.DN  # MODIS
    }
  }
  # Without DN (or for MODIS) LAI uses SR instead of TOAr
  if(missing(image.DN) || sat=="MODIS"){image.TOAr <- image.SR}
  albedo <- albedo(image.SR = image.SR, coeff = alb.coeff, sat = sat)
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
  } else if(lst == "L2"){
    message("Using the Level-2 surface temperature (image.SR$Band_10)")
    Ts <- image.SR$Band_10
  } else if(sat %in% c("L7", "L8", "L9")){
    Ts <- surfaceTemperature(LAI = LAI, sat = sat,
                             image.DN = if(!missing(image.DN)) image.DN,
                             method = lst, thermalband = thermalband,
                             MTL = MTL, WeatherStation = WeatherStation)
  } else {
    Ts <- image.DN$LST  # MODIS
  }
  checkKelvin(Ts, if(!missing(surf.temp)) "surf.temp" else "Ts")
  list(DEM = DEM, surface.model = surface.model, solar.angles = solar.angles,
       image.TOAr = image.TOAr, image.SR = image.SR, albedo = albedo,
       LAI = LAI, Ts = Ts)
}

# Net radiation from the output of metricSurface()
metricRn <- function(surface, WeatherStation, Kt = 1){
  Rs.inc <- incSWradiation(surface.model = surface$surface.model,
                           solar.angles = surface$solar.angles,
                           WeatherStation = WeatherStation, Kt = Kt)
  Rl.out <- outLWradiation(LAI = surface$LAI, Ts = surface$Ts)
  Rl.inc <- incLWradiation(WeatherStation, DEM = surface$surface.model$DEM,
                           solar.angles = surface$solar.angles, Ts = surface$Ts,
                           Kt = Kt)
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
#' @param image.SR         Surface reflectance image with layers B, G, R, NIR,
#' SWIR1 and SWIR2, e.g. from loadImageSR(). If supplied it is used directly;
#' if missing it is computed from image.DN (L7 via calcSR(), MODIS uses
#' image.DN). Required for L8 and L9.
#' Layers are selected by name; a stack without the standard names is used by
#' position, with a warning, and a stack with these names that lacks a needed
#' layer is an error.
#' @param WeatherStation   waterWeatherStation object with the conditions at
#' the satellite overpass (from read.WSdata() with the MTL file)
#' @param MTL              path to the Landsat metadata (*_MTL.txt) file.
#' Required
#' @param sat              Satellite: "L7", "L8", "L9", "MODIS" or "auto" to
#' detect it from the MTL file
#' @param thermalband      Thermal band in digital numbers (L7: band 6 low gain;
#' L8/L9: band 10). If missing, image.DN$Thermal1 is used. Not needed if
#' surf.temp is supplied
#' @param surf.temp        Optional. Surface temperature raster in Kelvin, on the
#' same grid as image.SR (e.g. image.SR$Band_10 from loadImageSR()). If
#' supplied it is used as Ts and surfaceTemperature() is not called.
#' Must be in Kelvin: values that look like Celsius, or more than 1\% of the
#' cells below 200 K (unmasked fill or clouds), stop; fewer cells below 200 K
#' give a warning.
#' @param plain            Logical. If TRUE (default) the surface is assumed
#' flat: DEM is replaced by a constant raster at the weather station elevation
#' @param DEM              SpatRaster with a digital elevation model (m) of the
#' study area, on the same grid as the image. Not needed if plain = TRUE
#' @param aoi              Optional sf polygon or SpatVector. If given, the
#' inputs are cropped to it
#' @param G.method         method used for the G estimation: "Tasumi" (Tasumi,
#' 2003; default) or "Bastiaanssen" (Bastiaanssen, 2000). See \code{\link{soilHeatFlux}}
#' @param alb.coeff        coefficient to transform narrow to broad band albedo.
#' See Details.
#' @param LST.method       Land surface temperature when surf.temp is not
#' supplied: "auto" (default) uses the Level-2 surface temperature
#' (image.SR$Band_10 from loadImageSR()) for L8 and L9 when present, and the
#' METRIC equations (Allen et al., 2007, Eqs. 20-21) with image.DN otherwise;
#' "L2" requires Band_10; "METRIC" uses the thermal band of image.DN; "SW"
#' (split window, L8 and L9) is deprecated. See \code{\link{surfaceTemperature}}
#' @param LAI.method       Method used to estimate LAI from spectral data.
#' See Details.
#' @param L                L value for SAVI calculation. See \code{\link{LAI}}
#' @param Zom.method       method selected to calculate momentum roughness
#' length. Use "short.crops" for short crops methods from Allen et al (2007);
#' "custom" for custom method also in Allen et al (2007); Or "Perrier" to use
#' Perrier equation as in Santos et al (2012) and Pocas et al (2014).
#' @param anchors.method   method for the automatic selection of the anchor
#' pixels: "percentile" (default), "flexible", "best" or "random". See
#' \code{\link{calcAnchors}}
#' @param anchors          data.frame with X, Y (coordinates in the CRS of the
#' image) and type ("hot" or "cold") columns, e.g. from calcAnchors(), or an sf
#' POINT object with a type column. If missing, the anchors are selected
#' automatically with calcAnchors(). They are checked by calcH(): type "hot" or
#' "cold" (case ignored), inside the image, inputs not NA there, and the hot
#' anchors at least 5 K warmer than the cold ones
#' @param n                number of pairs of anchor pixels to select (not
#' used by anchors.method = "percentile")
#' @param ETp.coef         ratio of the ET of the cold pixels to the
#' instantaneous alfalfa reference ET (ETr) at the satellite overpass. Allen et al. (2007)
#' use 1.05, the default. It can be updated with anchors.method = "flexible"
#' (see Details)
#' @param Z.om.ws          momentum roughness length (m) of the weather station
#' surface, used to extrapolate the wind speed to 200 m (Allen et al., 2007,
#' Eq. 32). It is about 0.123 times the height of the vegetation around the
#' station (Allen et al., 1998): 0.015 m (the default) for a standard station
#' over clipped grass 0.12 m tall, 0.03 m for vegetation about 0.25 m tall
#' @param verbose          Logical. If TRUE will print additional data to console
#' @param Kt               clearness (turbidity) coefficient of the atmosphere,
#' 0 < Kt <= 1: 1 (default) for clean air, 0.5 for extremely turbid, dusty or
#' polluted air (Allen et al., 2007). See \code{\link{SWtrasmisivity}}
#' @param extraParameters  named numeric vector with extra parameters for the
#' non default methods, e.g. c(fLAI = 0.5, h = 1.8) for Zom.method =
#' "Perrier". See \code{\link{momentumRoughnessLength}} and the Extra
#' Parameters section.
#' @param anchors.control  parameters of the automatic selection of anchor
#' pixels: a list from \code{\link{anchors.control}} (or a plain list with
#' some of its elements)
#' @param calcH.control    parameters of the iteration of H: a list from
#' \code{\link{calcH.control}} (or a plain list with some of its elements)
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
#' The ET of the cold pixels is ETp.coef times the instantaneous alfalfa reference ET
#' (ETr, ASCE standardized equation) at the satellite overpass (see \code{\link{calcH}}),
#' the same reference used by ET24h() and dailyET(), as in Allen et al. (2007).
#' With anchors.method = "flexible", if the cold-pixel criteria had to be
#' relaxed and the mean NDVI of the cold pixels is below 0.8, ETp.coef is
#' replaced by Kc * ETo / ETr, where Kc = 2.11 * NDVI - 0.4989 is a crop
#' coefficient on a grass reference basis (Reyes-González et al., 2015) and ETo / ETr
#' is the ratio of the grass and alfalfa reference ET at the overpass.
#'
#' Method arguments (sat, alb.coeff, LAI.method, LST.method and, where they
#' apply, G.method, Zom.method and anchors.method) must be one of the listed
#' names, exactly (case included): an unknown name stops before any
#' computation, with a hint for typos.
#' @section Extra Parameters:
#' extraParameters must be a named numeric vector with these names (unknown
#' names stop). a and b are required for Zom.method = "custom", and fLAI and h
#' for "Perrier":
#' \itemize{
#'   \item for momentumRoughnessLength() when Zom.method = "Perrier": fLAI, h;
#'   \item for momentumRoughnessLength() when Zom.method = "custom": a, b;
#'   \item deprecated: minDist, WSbuffer and deltaTemp for calcAnchors(); use
#'     anchors.control = anchors.control(...) instead.
#' }
#' @return An object of class \code{waterLSEB}, a list with the elements
#' below. A message summarises the calibration, and a warning is given when
#' ETrF is above 1.1 on more than 10\% of the pixels (cold anchors too warm):
#' \describe{
#'   \item{EB}{SpatRaster with layers NetRadiation, SoilHeat, SensibleHeat and
#'     LatentHeat (W/m2), and surfaceTemperature (K).}
#'   \item{WeatherStation}{the weather station data used.}
#'   \item{anchors}{sf POINT object with the hot and cold anchor pixels.}
#'   \item{calibration}{list with a and b (dT = a * Ts_datum + b), converged,
#'     iterations and u200 (wind speed at 200 m used, m/s), from calcH(); for
#'     automatic anchors, the number of candidates and the relaxation of the
#'     criteria, from calcAnchors(); and quality checks: Ts.hot and Ts.cold
#'     (mean Ts of the anchors, K, referred to the station elevation with
#'     plain = FALSE), dT.anchors (their difference), ETrF (5\%, 50\% and
#'     95\% quantiles and maximum of the reference ET fraction at the
#'     overpass, before LE is set to 0), ETrF.above.1.1, ETrF.below.0 and
#'     H.above.Rn.G (shares of the pixels). See \code{\link{summary.waterLSEB}}.}
#'   \item{terrain}{with plain = FALSE, SpatRaster with the surface model and
#'     the solar angles, for the slope correction of daily ET (see
#'     \code{\link{Crad}} and \code{\link{ET24h}}).}
#'   \item{methods}{named character vector with the methods and parameters used.}
#' }
#' @author Guillermo Federico Olmedo
#' @family METRIC model functions
#' @references
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#'
#' Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration: guidelines for computing crop water requirements. FAO Irrigation and Drainage Paper 56. FAO, Rome.
#'
#' Reyes-González, A., Hay, C., Kjaersgaard, J., & Neale, C. M. U. (2015). Use of remote sensing to generate crop coefficient and estimate actual crop evapotranspiration. 2015 ASABE Annual International Meeting, New Orleans, Louisiana. https://doi.org/10.13031/aim.20152190105
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
#'                             alb.coeff = "Olmedo", n = 5, Z.om.ws = 0.03)
#'
#' ### LSEB with "Perrier" method for Zom and extra parameters
#' Energy.Balance <- METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10,
#'                             WeatherStation = WeatherStation, MTL = MTLfile,
#'                             alb.coeff = "Olmedo", n = 5, Z.om.ws = 0.03,
#'                             Zom.method = "Perrier",
#'                             extraParameters = c(fLAI = 0.5, h = 1.8))
#' }
#' @export
METRIC.EB <- function(image.DN, image.SR, WeatherStation, MTL, sat = "auto",
                      thermalband, plain=TRUE, DEM, aoi, G.method = "Tasumi",
                      alb.coeff = "Tasumi", LST.method = "auto",
                      LAI.method = "metric2010", L = 0.1,
                      Zom.method = "short.crops", anchors.method = "percentile",
                      anchors, n = 1, ETp.coef= 1.05, Z.om.ws=0.015,
                      verbose = FALSE, extraParameters = vector(), surf.temp,
                      Kt = 1, anchors.control = NULL, calcH.control = NULL){
  checkMetricArgs(sat = sat, alb.coeff = alb.coeff, LAI.method = LAI.method,
                  LST.method = LST.method, G.method = G.method,
                  Zom.method = Zom.method)
  checkKt(Kt)
  anchors.method <- checkAnchorsMethod(anchors.method)
  checkExtraParameters(extraParameters, Zom.method)
  anchors.control <- controlList(anchors.control, "anchors.control")
  calcH.control <- controlList(calcH.control, "calcH.control")
  old <- intersect(c("deltaTemp", "minDist", "WSbuffer"), names(extraParameters))
  if(length(old)){
    message("extraParameters ", paste(old, collapse = ", "), " are deprecated: ",
            "use anchors.control = anchors.control(", old[1], " = ...)")
    anchors.control[old] <- as.list(extraParameters[old])
  }
  if(sat=="auto"){sat <- getSat(MTL = if(!missing(MTL)) MTL)}
  checkLSTmethod(LST.method, sat, missing(surf.temp))
  if(missing(anchors) && !missing(MTL)){
    checkCloudMask(MTL, if(!missing(image.SR)) image.SR, if(!missing(image.DN)) image.DN)
  }
  if(missing(image.SR) && (missing(image.DN) || !sat %in% c("L7", "MODIS"))){
    stop("image.SR is missing. It can only be computed from image.DN for L7 or ",
         "MODIS; provide it (e.g., loadImageSR()).")
  }


  MTL <- findMTL(if(!missing(MTL)) MTL)
  surface <- metricSurface(image.DN = image.DN, image.SR = image.SR,
                           WeatherStation = WeatherStation, MTL = MTL, sat = sat,
                           thermalband = thermalband, surf.temp = surf.temp,
                           alb.coeff = alb.coeff, LAI.method = LAI.method, L = L,
                           LST.method = LST.method, plain = plain, DEM = DEM,
                           aoi = aoi, Kt = Kt)
  DEM <- surface$DEM
  surface.model <- surface$surface.model
  image.TOAr <- surface$image.TOAr
  image.SR <- surface$image.SR
  albedo <- surface$albedo
  LAI <- surface$LAI
  Ts <- surface$Ts
  Rn <- metricRn(surface, WeatherStation, Kt = Kt)
  # Rn[Rn < 0]  <-  0 # see H
  #setTxtProgressBar(pb, 40)
  G <- soilHeatFlux(image = image.SR, Ts=Ts,albedo=albedo, 
                    Rn=Rn, LAI=LAI, method = G.method)
  # G[G < 0]  <-  0 # see H
  if(Zom.method == "short.crops"){
    Z.om <- momentumRoughnessLength(LAI=LAI, mountainous = !plain, slope.min1 = calcH.control$slope.min1, 
                                    method = Zom.method, 
                                    surface.model = surface.model)
  }
  if(Zom.method == "custom"){
    NDVI <- (image.TOAr$NIR - image.TOAr$R) / (image.TOAr$NIR + image.TOAr$R)
    Z.om <- momentumRoughnessLength(LAI=LAI, NDVI = NDVI, albedo = albedo,
                                    mountainous = !plain, slope.min1 = calcH.control$slope.min1,
                                    method = Zom.method, a = extraParameters["a"],
                                    b = extraParameters["b"],
                                    surface.model = surface.model)
  }
  if(Zom.method == "Perrier"){
    Z.om <- momentumRoughnessLength(LAI=LAI, mountainous = !plain, slope.min1 = calcH.control$slope.min1, 
                                    method = Zom.method, fLAI = extraParameters["fLAI"],
                                    h = extraParameters["h"],
                                    surface.model = surface.model)
  }
  if(Zom.method != "short.crops"){Z.om.sc <- momentumRoughnessLength(LAI=LAI, mountainous = !plain, slope.min1 = calcH.control$slope.min1, 
                                                             method = "short.crops", 
                                                             surface.model = surface.model)
  } else {Z.om.sc <- Z.om}
  par(mfrow=c(1,2))
  auto.anchors <- missing(anchors)
  # in mountain mode the anchors are chosen on Ts referred to the elevation of
  # the weather station, as the dT calibration of calcH()
  Ts.anchors <- if(plain) Ts else {
    Ts + (DEM - WeatherStation$location$elev) * calcH.control$lapse / 1000
  }
  if(auto.anchors){
    anchors <- muffleKelvin(
      calcAnchors(image = image.TOAr, Ts = Ts.anchors, LAI = LAI, plots = T,
                  albedo = albedo, Z.om = Z.om.sc, n = n,
                  anchors.method = anchors.method, WeatherStation = WeatherStation,
                  deltaTemp = anchors.control$deltaTemp,
                  minDist = anchors.control$minDist,
                  WSbuffer = anchors.control$WSbuffer,
                  min.candidates = anchors.control$min.candidates,
                  control = anchors.control,
                  verbose = verbose))
    if(nrow(anchors) > 20){
      # many anchors ("percentile"): their number and mean values
      print(do.call(rbind, lapply(split(anchors, anchors$type), function(a){
        data.frame(type = a$type[1], n = nrow(a), Ts = mean(a$Ts),
                   LAI = mean(a$LAI), NDVI = mean(a$NDVI))
      })), row.names = FALSE)
    } else print(anchors)
  }
  #setTxtProgressBar(pb, 45)
  if(!isTRUE(calcH.control$cold.NDVI.rule) && anchors.method == "flexible" &&
     !is.null(anchors$flex)){
     is.cold <- tolower(trimws(as.character(anchors$type))) == "cold"
     flex.cold <- anchors$flex[is.cold][1]
     NDVIcold <- mean(as.numeric(anchors$NDVI[is.cold]))
     if(flex.cold > 0 & NDVIcold < 0.8){
       print(paste0("mean NDVI for cold pixel: ", NDVIcold))
       # NDVI and Kc relationship from:
       # Reyes-González, Hay, Kjaersgaard and Neale (2015), alfalfa in Mexico.
       # doi:10.13031/aim.20152190105 
       # Kc on a grass reference (ETo) basis, converted to the alfalfa
       # reference (ETr) used for the cold pixels
       Kc <- 2.11 * NDVIcold - 0.4989
       ETp.coef <- Kc * hourlyET(WeatherStation, ET = "ETo") /
         hourlyET(WeatherStation, ET = "ETr")
       print(paste0("ETp.coef updated to: ", round(ETp.coef, 3), " (Kc = ",
                    round(Kc, 3), " from Reyes-Gonzalez et al, 2015, times ETo/ETr)"))
     }
  } else if(auto.anchors && isTRUE(calcH.control$cold.NDVI.rule)){
    # cold anchors that are not full cover: ET fraction from their NDVI, for
    # every method (calcH.control(cold.NDVI.rule = TRUE))
    is.cold <- tolower(trimws(as.character(anchors$type))) == "cold"
    NDVIcold <- mean(as.numeric(anchors$NDVI[is.cold]))
    if(NDVIcold < 0.8){
      Kc <- 2.11 * NDVIcold - 0.4989
      ETp.coef <- min(ETp.coef, Kc * hourlyET(WeatherStation, ET = "ETo") /
                        hourlyET(WeatherStation, ET = "ETr"))
      message("Mean NDVI of the cold anchors ", round(NDVIcold, 3), ": ETp.coef = ",
              round(ETp.coef, 3), " (cold.NDVI.rule)")
    }
  }
  if(calcH.control$ETrF.hot == 0){checkRecentRain(WeatherStation)}
  H <- muffleKelvin(
    calcH(anchors = anchors, Ts = Ts, Z.om = Z.om, mountainous = !plain,
          WeatherStation = WeatherStation, ETp.coef = ETp.coef,
          Z.om.ws = Z.om.ws, DEM = DEM, Rn = Rn, G = G, verbose = verbose,
          maxit = calcH.control$maxit, u200.min = calcH.control$u200.min,
          ETrF.hot = calcH.control$ETrF.hot, lapse = calcH.control$lapse,
          damping = calcH.control$damping))
  par(mfrow=c(1,1))
  #setTxtProgressBar(pb, 99)
  calibration <- H[c("a", "b", "converged", "iterations", "u200")]
  calibration$candidates <- attr(anchors, "candidates")
  calibration$relaxation <- attr(anchors, "relaxation")
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
  result$calibration <- c(calibration,
                          calibrationQC(Rn = Rn, G = G, H = H, Ts = Ts,
                                        Ts.anchors = Ts.anchors, anchors = anchors,
                                        WeatherStation = WeatherStation))
  reportCalibration(result$calibration)
  # for the slope correction of daily ET (Crad(), used by ET24h())
  if(!plain){result$terrain <- c(surface.model, surface$solar.angles)}
  result$methods <- c(sat = sat, alb.coeff = alb.coeff, LST.method = LST.method,
                      LAI.method = LAI.method, Zom.method = Zom.method, 
                      anchors.method = anchors.method, plain = plain,
                      ETp.coef= ETp.coef, Z.om.ws=Z.om.ws, Kt = Kt,
                      extraParameters = extraParameters,
                      unlist(anchors.control), unlist(calcH.control))
  class(result) <- "waterLSEB"
  return(result)
}
  
