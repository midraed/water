#' Create aoi polygon from topleft and bottomright coordinates
#' @description
#' An AOI (Area of Interest) is created based on two points (topleft and bottomright) using a coordinate reference system.
#' @param topleft     a vector with topleft x,y coordinates
#' @param bottomright a vector with bottomright x,y coordinates
#' @param EPSG        EPSG code of the coordinate reference system of the
#' coordinates (e.g. 32619 for WGS 84 / UTM zone 19N). Optional: if missing,
#' the CRS is not set
#' @return An sf data frame with a single POLYGON geometry (the rectangle
#' defined by topleft and bottomright). Its CRS is set only when EPSG is
#' supplied.
#' @family support functions
#' @author Guillermo Federico Olmedo
#' @author David Fonseca-Luengo
#' @examples 
#' tl <- c(493300, -3592700)
#' br <- c(557200, -3700000) 
#' aoi <- createAoi(topleft = tl, bottomright=br, EPSG=32619)
#' plot(aoi)
#' @import terra sf
#' @export
# Maybe i can provide some points and use CHull like in QGIS-Geostat
createAoi <- function(topleft, bottomright, EPSG){
  coords <- matrix(c(topleft[1],bottomright[1], bottomright[1],topleft[1],topleft[1],
                     topleft[2], topleft[2], bottomright[2], bottomright[2],topleft[2]),
                   ncol=2, nrow=5)
  aoi <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(coords))))
  if(!missing(EPSG)){sf::st_crs(aoi) <- paste0("EPSG:", EPSG)}
  return(aoi)
}

#' Check needed SRTM grids from image extent
#' @description
#' Lists the 1 x 1 degree SRTM tiles needed to cover the extent of an image,
#' with the file names that prepareSRTMdata() reads. The tiles can be
#' downloaded from EarthExplorer (dataset "SRTM 1 Arc-Second Global"). To
#' download a DEM automatically, use getDEM().
#' @param raw.image  SpatRaster whose extent (and CRS) is used to find the SRTM
#' tiles needed
#' @return Character vector with the file names of the SRTM tiles needed (e.g.
#' "s34_w069_1arc_v3.tif"). The number of tiles is also reported.
#' @author Guillermo Federico Olmedo
#' @export
#' @family support functions
# Get links or optionally open web pages... 
# Check if the files are present on path o in a specific SRTM local repo
checkSRTMgrids <-function(raw.image){
  tiles <- demTiles(raw.image)
  files <- sprintf("%s%02d_%s%03d_1arc_v3.tif", ifelse(tiles$lat < 0, "s", "n"),
                   abs(tiles$lat), ifelse(tiles$lon < 0, "w", "e"), abs(tiles$lon))
  message("You need ", length(files), " SRTM tile(s) of 1 x 1 degree (EarthExplorer, ",
          "dataset \"SRTM 1 Arc-Second Global\"), or use getDEM() to download a ",
          "DEM automatically")
  files
}

#' Download a digital elevation model for an image
#' @description
#' Downloads the Copernicus DEM GLO-30 (30 m) tiles that cover the image and
#' returns the DEM on the image grid, ready for METRIC.EB(plain = FALSE, DEM =
#' ). The tiles are Cloud Optimized GeoTIFFs in the AWS Open Data registry,
#' read over HTTPS: no account is needed, and only the area of the image is
#' downloaded. It needs an internet connection. To use SRTM tiles already
#' downloaded, see \code{\link{checkSRTMgrids}} and prepareSRTMdata().
#' @param image   SpatRaster (e.g. the Landsat image) that defines the extent,
#' CRS and grid of the DEM
#' @param source  DEM source. Currently only "copernicus" (Copernicus DEM
#' GLO-30)
#' @return SpatRaster with one layer named DEM: elevation (m, relative to the
#' EGM2008 geoid) on the grid of image.
#' @references
#' Copernicus DEM GLO-30: (c) DLR e.V. 2010-2014 and (c) Airbus Defence and
#' Space GmbH 2014-2018, provided under COPERNICUS by the European Union and
#' ESA; all rights reserved. Distributed through the Registry of Open Data on
#' AWS, \url{https://registry.opendata.aws/copernicus-dem/}
#' @author Guillermo Federico Olmedo
#' @family support functions
#' @examples
#' \dontrun{
#' image.DN <- rast(L7_Talca)
#' DEM <- getDEM(image.DN)
#' plot(DEM)
#' }
#' @export
getDEM <- function(image, source = "copernicus"){
  source <- match.arg(source)
  if(!inherits(image, "SpatRaster")){
    stop("image must be a SpatRaster, e.g. the Landsat image")
  }
  if(terra::crs(image) == ""){stop("image has no CRS")}
  tiles <- demTiles(image, margin = 0.01)
  lat <- sprintf("%s%02d", ifelse(tiles$lat < 0, "S", "N"), abs(tiles$lat))
  lon <- sprintf("%s%03d", ifelse(tiles$lon < 0, "W", "E"), abs(tiles$lon))
  name <- paste0("Copernicus_DSM_COG_10_", lat, "_00_", lon, "_00_DEM")
  url <- paste0("/vsicurl/https://copernicus-dem-30m.s3.amazonaws.com/", name, "/",
                name, ".tif")
  area <- ext(terra::project(terra::as.polygons(ext(image), crs = terra::crs(image)),
                             "EPSG:4326"))
  area <- ext(area[1] - 0.01, area[2] + 0.01, area[3] - 0.01, area[4] + 0.01)
  parts <- list()
  for(u in url){
    # tiles over the ocean do not exist
    tile <- tryCatch(crop(rast(u), area), error = function(e) NULL)
    if(!is.null(tile)){parts[[length(parts) + 1]] <- tile}
  }
  if(!length(parts)){
    stop("No Copernicus DEM tiles could be read for the image (only ocean, or ",
         "no internet connection)")
  }
  DEM <- if(length(parts) == 1) parts[[1]] else terra::merge(terra::sprc(parts))
  DEM <- project(DEM, image[[1]], method = "bilinear")
  DEM <- saveLoadClean(imagestack = DEM, stack.names = "DEM", file = "DEM",
                       overwrite = TRUE)
  return(DEM)
}

#' Create a DEM mosaic from SRTM grids for the image extent
#' @description
#' Mosaics the SRTM 1 arc-second tiles found in path (see \code{\link{checkSRTMgrids}}) and
#' projects the result to the grid of the image.
#' @param format  file extension of the SRTM grid files (default "tif"); files
#' must be named like s34_w069_1arc_v3.tif
#' @param extent  SpatRaster (e.g. the Landsat image) whose grid is used for the
#' output DEM, or an sf / SpatVector polygon with a CRS, in which case a 30 m
#' grid in that CRS is used
#' @param path  folder where SRTM files are stored
#' @return SpatRaster with one layer named DEM (m), on the grid defined by
#' extent.
#' @author Guillermo Federico Olmedo
#' @export
#' @family support functions
# Should use checkSRTMgrids to get the files list and not use all from the folder...!
# Also look for files on path and local repo
prepareSRTMdata <- function(path=getwd(), format="tif", extent){
  files <- list.files(path= path,  
                      pattern=paste("^[sn]\\d{2}_[we]\\d{3}_1arc_v3.", 
                              format, "$", sep=""), full.names = T) 
  if(length(files) < 1){stop(paste("You need to download SRTM grids and save them to working directory \n try: checkSRTMgrids()"))}
  stack1 <- list()
  for(i in 1:length(files)){
    stack1[[i]] <- rast(files[i])}
  stack1$fun <- "mean"
  if(length(files)>1){SRTMmosaic <- do.call(mosaic, stack1)}
  if(length(files)==1){SRTMmosaic <- stack1[[1]]}
  # Target grid: the image grid, or a 30 m grid in the CRS of a polygon
  if(inherits(extent, "SpatRaster")){
    destino <- rast(extent[[1]])
  } else if(inherits(extent, c("sf", "sfc", "SpatVector"))){
    extent <- terra::vect(extent)
    if(terra::crs(extent) == ""){
      stop("extent has no CRS; use createAoi(..., EPSG = ) or pass the image")
    }
    if(terra::is.lonlat(extent)){
      stop("extent is in longitude/latitude; pass the image or a polygon in ",
           "the projected CRS of the image (e.g. UTM)")
    }
    e <- as.vector(ext(extent))
    snap <- function(v, f){f((v - 15) / 30) * 30 + 15}
    destino <- rast(ext(snap(e[1], floor), snap(e[2], ceiling),
                        snap(e[3], floor), snap(e[4], ceiling)),
                    crs = terra::crs(extent), resolution = 30)
  } else {
    stop("extent must be a SpatRaster (e.g. the Landsat image) or an sf / ",
         "SpatVector polygon with a CRS")
  }
  mosaicp <- project(SRTMmosaic, destino, method = "bilinear")
  mosaicp <- saveLoadClean(imagestack = mosaicp, stack.names = "DEM", 
                           file = "DEM", overwrite=TRUE)
  return(mosaicp)
}

#' Calculates surface model used in METRIC
#' @description
#' The DEM is used to generate the surface representation of the image through
#' aspect and slope maps. This procedure helps to avoid differences in the
#' surface temperature (and finally evapotranspiration) caused by different
#' incidence angles and/or elevations in mountainous areas. A constant DEM
#' (plain = TRUE in METRIC.EB()) gives a flat surface everywhere, without an NA
#' border.
#' @param DEM  SpatRaster with a digital elevation model (m)
#' @return SpatRaster with layers DEM (m), Slope and Aspect (radians; Aspect is
#' shifted by -pi, as METRIC expects).
#' @author Guillermo Federico Olmedo
#' @author David Fonseca-Luengo
#' @family net radiation related functions
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
METRICtopo <- function(DEM){
  aspect <- terrain(DEM, v="aspect", unit="radians")
  slope <- terrain(DEM, v="slope", unit="radians")
  # terrain() leaves the outer cells NA; a constant DEM (plain = TRUE) is flat
  # there too, so use the values terrain() gives inside a flat surface
  if(isTRUE(diff(as.numeric(global(DEM, "range", na.rm = TRUE))) == 0)){
    slope <- ifel(is.na(DEM), NA, 0)
    aspect <- ifel(is.na(DEM), NA, 3*pi/2)
  }
  aspect_metric <- aspect-pi  #METRIC expects aspect - 1 pi
  surface.model <- c(DEM, slope, aspect_metric)
  surface.model <- saveLoadClean(imagestack = surface.model, 
                                 stack.names = c("DEM", "Slope", "Aspect"), 
                                 file = "surface.model", 
                                 overwrite=TRUE)
  return(surface.model)
}

#' Calculates solar angles
#' @description
#' Metadata, aspect and slope maps are combined to estimate solar angles for the entire image. 
#' @param surface.model   SpatRaster with layers DEM, Slope and Aspect, from
#' METRICtopo()
#' @param MTL             path to the Landsat metadata (*_MTL.txt) file.
#' Required
#' @param WeatherStation  waterWeatherStation object (see \code{\link{read.WSdata}}) or the
#' data.frame returned by getDataWS()
#' @details Latitude is computed for each row of the grid; declination and hour
#' angle are computed from the acquisition date and scene centre time in the
#' MTL file and the weather station longitude. Solar incidence angles are
#' computed for a horizontal surface (incidence.hor) and for the sloping
#' surface given by Slope and Aspect (incidence.rel), following Allen et al.
#' (2007).
#' @return SpatRaster with layers (radians) latitude, declination, hour.angle
#' (negative before solar noon), incidence.hor (solar incidence angle on a
#' horizontal surface) and incidence.rel (incidence angle on the sloping
#' surface).
#' @author Guillermo Federico Olmedo
#' @author David Fonseca-Luengo
#' @author Fernando Fuentes-Peñailillo
#' @family net radiation related functions
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
### Change to look in metadata for keyword instead of using line #
solarAngles <- function(surface.model, MTL, WeatherStation){
  if(inherits(WeatherStation, "waterWeatherStation")){
    WeatherStation <- getDataWS(WeatherStation)
  }
  MTL <- readLines(findMTL(if(!missing(MTL)) MTL), warn=FALSE)
  Elev.line <- grep("SUN_ELEVATION",MTL,value=TRUE)
  sun.elevation <- (90 - as.numeric(regmatches(Elev.line, 
                                               regexec(text=Elev.line ,
                                                       pattern="([0-9]{1,5})([.]+)([0-9]+)"))[[1]][1]))*pi/180
  Azim.line <- grep("SUN_AZIMUTH",MTL,value=TRUE)
  sun.azimuth <- as.numeric(regmatches(Azim.line, 
                                       regexec(text=Azim.line ,
                                               pattern="([0-9]{1,5})([.]+)([0-9]+)"))[[1]][1])
  # latitude
  latitude <- surface.model[[1]]
  xy <- terra::xyFromCell(latitude, terra::cellFromRowCol(latitude, 1:nrow(latitude), 1))
  xy <- sf::st_as_sf(data.frame(x = xy[,1], y = xy[,2]), coords = c("x", "y"), crs = terra::crs(latitude))
  lat <- sf::st_coordinates(sf::st_transform(xy, 4326))[,2]
  values(latitude) <- rep(lat*pi/180,each=ncol(latitude))
  # declination
  time.line <- grep("SCENE_CENTER_TIME",MTL,value=TRUE)
  date.line <- grep("DATE_ACQUIRED",MTL,value=TRUE)
  sat.time <-regmatches(time.line,regexec(text=time.line,
                                          pattern="([0-9]{2})(:)([0-9]{2})(:)([0-9]{2})(.)([0-9]{2})"))[[1]][1]
  sat.date <-regmatches(date.line,regexec(text=date.line,
                                          pattern="([0-9]{4})(-)([0-9]{2})(-)([0-9]{2})"))[[1]][1]
  sat.datetime <- strptime(paste(sat.date, sat.time), 
                           format = "%Y-%m-%d %H:%M:%S", tz="GMT")
  DOY <-  sat.datetime$yday +1
  declination <- surface.model[[1]]
  values(declination) <- 0.409*sin((2*pi/365*DOY)-1.39)
  # hour angle
  hour.angle <- surface.model[[1]]
  time.line <- grep("SCENE_CENTER_TIME",MTL,value=TRUE)
  time.hours <-regmatches(time.line,regexec(text=time.line,
                                            pattern="([0-9]{2})(:)([0-9]{2})(:)([0-9]{2})(.)([0-9]{2})"))[[1]][1]
  time.hours <- strptime(time.hours, format = "%H:%M:%S", tz="GMT")
  Sc <- (0.1645*sin(4*pi*(DOY-81)/364)-0.1255*cos(2*pi*(DOY-81)/364)-0.025*sin(2*pi*(DOY-81)/364))*60
  time.hours <- (time.hours$hour*3600+time.hours$min*60+time.hours$sec+WeatherStation$long/15*3600)/60 + Sc
  values(hour.angle) <- (as.numeric(time.hours)/60-12)*15/180*pi
  #hour.angle <- asin(-1*(cos(sun.elevation)*sin(sun.azimuth)/cos(declination))) first
  
  ## solar incidence angle, for horizontal surface
  incidence.hor <- acos(sin(declination) * sin(latitude) + cos(declination)*
                          cos(latitude)*cos(hour.angle))
  slope <- surface.model$Slope
  aspect <- surface.model$Aspect
  ##solar incidence angle, for sloping surface
  incidence.rel <- acos(sin(declination)*sin(latitude)*cos(slope) 
                        - sin(declination)*cos(latitude)*sin(slope)*cos(aspect)
                        + cos(declination)*cos(latitude)*cos(slope)*cos(hour.angle)
                        + cos(declination)*sin(latitude)*sin(slope)*cos(aspect)*cos(hour.angle)
                        + cos(declination)*sin(aspect)*sin(slope)*sin(hour.angle))
  ## End
  solarAngles <- c(latitude, declination, hour.angle, incidence.hor, incidence.rel)
  solarAngles <- saveLoadClean(imagestack = solarAngles, 
                               stack.names = c("latitude", "declination", 
                                               "hour.angle", "incidence.hor", "incidence.rel"), 
                               file = "solarAngles", overwrite=TRUE)
  return(solarAngles)
}

#' Calculates Incoming Solar Radiation
#' @description 
#' This function calculates incoming solar radiation from the surface model,
#' the solar angles and the weather station data.
#' @param surface.model   SpatRaster with layers DEM, Slope and Aspect, from
#' METRICtopo()
#' @param solar.angles    SpatRaster with layers latitude, declination,
#' hour.angle, incidence.hor and incidence.rel (radians), from solarAngles()
#' @param WeatherStation  waterWeatherStation object (see \code{\link{read.WSdata}}) or the
#' data.frame returned by getDataWS()
#' @param Kt              clearness (turbidity) coefficient of the atmosphere,
#' 0 < Kt <= 1: 1 (default) for clean air, 0.5 for extremely turbid, dusty or
#' polluted air (Allen et al., 2007). See \code{\link{SWtrasmisivity}}
#' @return SpatRaster with incoming short-wave radiation (W/m2) at the time of
#' the satellite overpass.
#' @author Guillermo Federico Olmedo
#' @author Daniel de la Fuente-Sáiz
#' @author David Fonseca-Luengo
#' @family net radiation related functions
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
incSWradiation <- function(surface.model, solar.angles, WeatherStation, Kt = 1){
  checkKt(Kt)
  if(inherits(WeatherStation, "waterWeatherStation")){
    WeatherStation <- getDataWS(WeatherStation)
  }
  ea <- (WeatherStation$RH/100)*0.6108*exp((17.27*WeatherStation$temp)/(WeatherStation$temp+237.3))
  tau.sw <- SWtrasmisivity(Kt = Kt, ea, surface.model$DEM, solar.angles$incidence.hor)
  DOY <- WeatherStation$DOY
  d <- sqrt(1/(1+0.033*cos(DOY * (2 * pi/365))))
  Rs.inc <- 1367 * cos(solar.angles$incidence.rel) * tau.sw / d^2
  Rs.inc <- saveLoadClean(imagestack = Rs.inc, 
                          file = "Rs.inc", overwrite=TRUE)
  return(Rs.inc)
}

#' Calculates Broadband Albedo from Landsat data
#' @description 
#' Broadband surface albedo is estimated considering the integration of all
#' narrowband at-surface reflectances following a weighting function with
#' empirical coefficients (Tasumi et al., 2008).
#' @param image.SR   surface reflectance SpatRaster with layers named B, G, R,
#' NIR, SWIR1 and SWIR2 (MODIS: R, NIR, B, G, SWIR1, SWIR2 and SWIR3), selected
#' by name. Layers without these names are used by position, with a warning
#' @param aoi        Optional area of interest (e.g. from createAoi()) to crop
#' the result. If missing, the result is not cropped
#' @param coeff      coefficient to transform narrow to broad band albedo.
#' See Details.
#' @param sat        "L7", "L8", "L9", "MODIS" or "auto": MODIS when there is a
#' layer named SWIR3 (as from loadImage()), Landsat otherwise (L7, L8 and L9
#' use the same coefficients)
#' @details
#' There are different models to convert narrowband data to broadband albedo.
#' You can choose coeff = "Tasumi" to use Tasumi et al. (2008) coefficients,
#' calculated for Landsat 7; coeff = "Liang" to use Liang (2001) Landsat 7
#' coefficients or coeff = "Olmedo" to use Olmedo coefficients for Landsat 8
#' and 9.
#' @return SpatRaster with broadband surface albedo (unitless).
#' @author Guillermo Federico Olmedo
#' @author David Fonseca-Luengo
#' @family net radiation related functions
#' @references
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' Tasumi, M., Allen, R. G., & Trezza, R. (2008). At-surface reflectance and albedo from satellite for operational calculation of land surface energy balance. Journal of Hydrologic Engineering, 13(2), 51-63. https://doi.org/10.1061/(ASCE)1084-0699(2008)13:2(51) \cr
#' Liang, S. (2001). Narrowband to broadband conversions of land surface albedo I: Algorithms. Remote Sensing of Environment, 76(2), 213-238. https://doi.org/10.1016/S0034-4257(00)00205-4 \cr
#' @export
albedo <- function(image.SR, aoi, coeff="Tasumi", sat="auto"){
  checkChoice(coeff, waterChoices$alb.coeff, "coeff")
  checkChoice(sat, waterChoices$sat, "sat")
  if(sat=="auto"){
    # L7, L8 and L9 use the same coefficients; only MODIS has SWIR3
    sat <- if("SWIR3" %in% names(image.SR)) "MODIS" else "L8"
  }
  if(sat %in% c("L7", "L8", "L9")){
    if(coeff=="Tasumi"){wb <- c(0.254, 0.149, 0.147, 0.311, 0.103, 0.036)}# Tasumi 2008
    if(coeff=="Olmedo") {wb <- c(0.246, 0.146, 0.191, 0.304, 0.105, 0.008)} # Calculated using SMARTS for Kimberly2-noc13 and Direct Normal Irradiance
    if(coeff=="Liang") {wb <- c(0.356, 0, 0.130, 0.373, 0.085, 0.072)} # Liang 2001
  }
  if(sat=="MODIS"){wb <- c(0.037, 0.479, -0.068, 0.0976, 0.266, 0.0757, 0.107)} # Liang 2001 direct-NIR coefficients
  # bands by name, in the order of the coefficients
  image.SR <- if(sat == "MODIS"){
    landsatBands(image.SR, c("R", "NIR", "B", "G", "SWIR1", "SWIR2", "SWIR3"))
  } else landsatBands(image.SR, c("B", "G", "R", "NIR", "SWIR1", "SWIR2"))
  albedo <- app(image.SR[[1]], fun=function(x){x *wb[1]})
    for(i in 2:length(wb)){
      albedo <- albedo + app(image.SR[[i]], fun=function(x){x *wb[i]})
  }
  albedo <- aoiCrop(albedo, aoi) 
  if(coeff=="Liang"){
    albedo <- albedo - 0.0018
  }
  albedo <- saveLoadClean(imagestack = albedo, 
                          file = "albedo", overwrite=TRUE)
  return(albedo)
}

#' Estimate LAI from Landsat Data
#' @description
#' This function implements empirical models to estimate LAI (Leaf Area Index)
#' from satellite images. Models were extracted from METRIC publications and
#' other works developed on different crops.
#' @param method   Method used to estimate LAI from spectral data: "metric2010"
#' (default), "metric", "MCB", "vineyard" or "turner". See Details.
#' @param image    SpatRaster with layers named R and NIR, selected by name
#' (layers without standard names: the third and fourth, with a warning):
#' top-of-atmosphere reflectance for "metric",
#' "metric2010" and "MCB"; surface reflectance for "turner"; radiance for
#' "vineyard" (see \code{\link{calcRadiance}})
#' @param aoi      Not used; kept for backward compatibility
#' @param L        L factor used in method = "metric" or "metric2010" to estimate SAVI, defaults to 0.1
#' @param sat      Not used; kept for backward compatibility
#' @details LAI is computed from top-of-atmosphere reflectance ("metric",
#' "metric2010" and "MCB"), surface reflectance ("turner") or radiance
#' ("vineyard"). LAI and other indices such as NDVI and SAVI are used to
#' describe the vegetation, depending on the preferences of the user.
#' Available methods are: "metric", "metric2010", "MCB", "vineyard" and "turner".
#'
#' "turner" inverts the relationship between surface reflectance NDVI and LAI
#' of Turner et al. (1999, Table 3), NDVI = 0.5724 + 0.0989 LAI - 0.0114 LAI^2 +
#' 0.0004 LAI^3, fitted for forest, grassland and shrub plots: LAI is 0 for NDVI
#' below 0.572 and 6 for NDVI above 0.842.
#'
#' When METRIC.EB(), METRIC.Rn() or METRIC.G() run without image.DN, the
#' "metric", "metric2010" and "MCB" methods receive surface reflectance,
#' although their equations were described for top-of-atmosphere reflectance.
#' @return SpatRaster with leaf area index (m2/m2) in a layer named LAI;
#' negative values are set to 0.
#' @author Guillermo Federico Olmedo
#' @author David Fonseca-Luengo
#' @author Fernando Fuentes-Peñailillo
#' @family net radiation related functions
#' @references
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' Carrasco-Benavides, M., Ortega-Farias, S., Lagos, L., Kleissl, J., Morales-Salinas, L., & Kilic, A. (2014). Parameterization of the Satellite-Based Model (METRIC) for the Estimation of Instantaneous Surface Energy Balance Components over a Drip-Irrigated Vineyard. Remote Sensing, 6(11), 11342-11371. http://doi.org/10.3390/rs61111342\cr
#' Johnson, L. F. (2003). Temporal stability of an NDVI-LAI relationship in a Napa Valley vineyard. Australian Journal of Grape and Wine Research, 9(2), 96-101. https://doi.org/10.1111/j.1755-0238.2003.tb00258.x \cr
#' Turner, D. P., Cohen, W. B., Kennedy, R. E., Fassnacht, K. S., & Briggs, J. M. (1999). Relationships between leaf area index and Landsat TM spectral vegetation indices across three temperate zone sites. Remote Sensing of Environment, 70(1), 52–68. http://doi.org/10.1016/S0034-4257(99)00057-7
#' @export
## Cite Pocas work for LAI from METRIC2010
LAI <- function(method="metric2010", image, aoi, L=0.1, sat){
  checkChoice(method, waterChoices$LAI.method, "method")
  if(!missing(sat)){warning("LAI(): the sat argument is not used")}
  # R and NIR by name (layers 3 and 4 of an unnamed stack)
  if(method=="metric" | method=="metric2010" | method=="vineyard" | method=="MCB"){
      toa.4.5 <- landsatBands(image, c("R", "NIR"), pos = 3:4)}
  if(method=="turner"){
      sr.4.5 <- landsatBands(image, c("R", "NIR"), pos = 3:4)}
  if(method=="metric2010"){
    # SAVI
    SAVI_ID <- (1 + L) * (toa.4.5[[2]] - toa.4.5[[1]]) / (L + toa.4.5[[1]] + toa.4.5[[2]])
    SAVI_ID <- saveLoadClean(imagestack = SAVI_ID, stack.names = "SAVI_ID",
                             file = "SAVI_ID", overwrite = TRUE)
    
    # LAI from SAVI
    LAI <- ifel(
      SAVI_ID <= 0, 0,
      ifel(SAVI_ID <= 0.817, 11 * (SAVI_ID^3), 6)
    )
    names(LAI) <- "LAI"
  }
  if(method=="metric"){
    SAVI_ID <- (1 + L)*(toa.4.5[[2]] - toa.4.5[[1]])/(L + toa.4.5[[1]] + toa.4.5[[2]])
    LAI <- ifel(SAVI_ID > 0.687, 6, log((0.69-SAVI_ID)/0.59)/0.91 *-1)
  }
  if(method=="vineyard"){
    # image must be the radiance image (calcRadiance())
    NDVI <- (toa.4.5[[2]] - toa.4.5[[1]])/(toa.4.5[[1]] + toa.4.5[[2]])
    LAI <- 4.9 * NDVI -0.46 # Johnson 2003
  }
  ## method carrasco
  if(method=="MCB"){
    NDVI <- (toa.4.5[[2]] - toa.4.5[[1]])/(toa.4.5[[1]] + toa.4.5[[2]])
    LAI <- 1.2 - 3.08*exp(-2013.35*NDVI^6.41) 
  }
  if(method=="turner"){
    NDVI <- (sr.4.5[[2]] - sr.4.5[[1]])/(sr.4.5[[1]] + sr.4.5[[2]])
    # Turner et al. (1999, Table 3) regressed the surface reflectance NDVI on
    # LAI: NDVI = 0.5724 + 0.0989 LAI - 0.0114 LAI^2 + 0.0004 LAI^3, which
    # increases up to LAI = 6.7. Inverted here for LAI between 0 and 6
    LAI.grid <- seq(0, 6, by = 0.001)
    NDVI.grid <- 0.5724 + 0.0989*LAI.grid - 0.0114*LAI.grid^2 + 0.0004*LAI.grid^3
    LAI <- app(NDVI, fun = function(x){
      stats::approx(NDVI.grid, LAI.grid, xout = x, rule = 2)$y
    })
  }
  LAI[LAI<0]  <- 0
  LAI <- saveLoadClean(imagestack = LAI, stack.names = "LAI", 
                       file = "LAI", overwrite=TRUE)
  return(LAI)
}

#' Calculates broadband short-wave transmissivity
#' @description
#' Short-wave transmissivity is estimated for broad-band considering an
#' extended equation developed by Allen (1996), based on Majumdar et al.
#' (1972), using coefficients developed by ASCE-EWRI (2005).
#' @param Kt            unitless turbidity coefficient 0<Kt<=1.0, where Kt=1.0 for clean air and Kt=0.5 for extremely turbid, dusty, or polluted air
#' @param ea            near-surface vapor pressure (kPa)
#' @param dem           SpatRaster with a digital elevation model (m)
#' @param incidence.hor solar incidence angle (radians) on a horizontal surface,
#' from solarAngles()
#' @return SpatRaster with broadband short-wave atmospheric transmissivity
#' (unitless).
#' @author Guillermo Federico Olmedo
#' @author David Fonseca-Luengo
#' @family net radiation related functions
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#'
#' Majumdar, N.; Mathur, B. & Kaushik, S. Prediction of direct solar radiation for low atmospheric turbidity Solar Energy, Elsevier, 1972, 13, 383-394 \cr
#'
#' ASCE-EWRI The ASCE Standardized Reference Evapotranspiration Equation Report of the ASCE-EWRI Task Committee on Standardization of Reference Evapotranspiration, 2005 \cr
#' @export
SWtrasmisivity <- function(Kt = 1, ea, dem, incidence.hor){
  P <- 101.3*((293-0.0065 * dem)/293)^5.26
  W <- 0.14 * ea * P + 2.1
  tauB <- 0.98 * exp(((-0.00149 * P )/ (Kt * 
                      cos(incidence.hor)))-0.075*((W / cos(incidence.hor))^0.4))
  tauD <- rast(tauB)
  tauD <- 0.35 - 0.36 * tauB 
  tauD[tauB < 0.15] <- 0.18 + 0.82 * tauB[tauB < 0.15]
  tau.sw <- tauB + tauD
  # Next one it's from METRIC 2007, previous from METRIC 2010
  #sw.t <- 0.35 + 0.627 * exp((-0.00149 * P / Kt * 
  #                              cos(incidence.hor))-0.075*(W / cos(incidence.hor))^0.4)
  return(tau.sw)
}

#' Estimates Land Surface Temperature from Landsat Data
#' @description
#' Surface temperature is estimated using a modified Planck equation
#' considering empirical constants for Landsat images. In addition, this model
#' implements a correction of thermal radiance following Wukelic et al.
#' (1989). For Landsat 8 and 9, a single-channel (Jimenez-Munoz et al., 2014)
#' or a split-window algorithm can be used. With Landsat 8/9 Collection 2
#' Level-2 products, the surface temperature band can be used directly instead
#' (surf.temp in METRIC.EB()).
#' @param image.DN        raw image in digital numbers, e.g. from loadImage()
#' @param thermalband     thermal band in digital numbers. Optional: if missing,
#' image.DN$Thermal1 is used (L7: band 6 low gain; L8/L9: band 10)
#' @param sat             "L7", "L8", "L9" or "auto" to detect it from the MTL
#' file (MTL is then required)
#' @param LAI             SpatRaster with leaf area index, from LAI()
#' @param aoi             Not used; kept for backward compatibility
#' @param method          "SC" for single channel or "SW" for split window
#'                        algorithm. "SW" is only available for L8 and L9
#' @param WeatherStation  waterWeatherStation object (see \code{\link{read.WSdata}}) or the
#' data.frame returned by getDataWS()
#' @param MTL             path to the Landsat metadata (MTL) file. For L8 and
#'                        L9 the thermal calibration constants are read from
#'                        it. Required for L9; without it L8 uses Landsat 8
#'                        defaults
#' @return SpatRaster with land surface temperature (K).
#' @family net radiation related functions
#' @author Guillermo Federico Olmedo
#' @author David Fonseca-Luengo
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for 
#' mapping evapotranspiration with internalized calibration (METRIC) - 
#' Model" Journal of Irrigation and Drainage Engineering, vol. 133, 
#' p. 380, 2007 \cr
#'
#' Wukelic G. E.; Gibbons D. E.; Martucci L. M. & Foote, H. P. Radiometric 
#' calibration of Landsat thematic mapper thermal band Remote Sensing of 
#' Environment, 1989, 28, (339-347) \cr
#'
#' Jimenez-Munoz, J. C., Sobrino, J. A., Skokovic, D., Mattar, C., & Cristobal, 
#' J. (2014). Land surface temperature retrieval methods from Landsat-8 thermal
#' infrared sensor data.  IEEE Geoscience and Remote Sensing Letters, 11(10), 
#' 1840–1843. http://doi.org/10.1109/LGRS.2014.2312032 \cr 
#' @export
## Add Sobrino and Qin improvements to LST in ETM+
surfaceTemperature <- function(image.DN, sat="auto", LAI, aoi, 
                               method = "SC", WeatherStation, thermalband, MTL){
  checkChoice(sat, c("auto", "L7", "L8", "L9"), "sat")
  checkChoice(method, waterChoices$LST.method, "method")
  checkLSTmethod(method, sat, TRUE, "method")
  if(inherits(WeatherStation, "waterWeatherStation")){
    WeatherStation <- getDataWS(WeatherStation)
  }
  epsilon_NB <- 0.97 + 0.0033 * LAI  
  epsilon_NB[LAI > 3] <- 0.98
  # Downwelling atm radiance
  
  R_sky <- (1.807e-10)*(WeatherStation$temp+273.15)^4*(1-0.26*
                                            exp(-7.77e-4*WeatherStation$temp^2))
  # path radiance
  Rp <- 0.91     
  # atm transmissivity
  tau_NB <- 0.75+2e-5*WeatherStation$elev
  if(sat=="auto"){sat = getSat(MTL = if(!missing(MTL)) MTL)}
  checkLSTmethod(method, sat, TRUE, "method")
  if(sat %in% c("L8", "L9")){
    # Landsat 8 defaults; L9 (TIRS-2) constants differ, so read them from MTL
    cal <- c(RADIANCE_MULT_BAND_10 = 3.3420E-04, RADIANCE_ADD_BAND_10 = 0.1,
             RADIANCE_MULT_BAND_11 = 3.3420E-04, RADIANCE_ADD_BAND_11 = 0.1,
             K1_CONSTANT_BAND_10 = 774.8853, K2_CONSTANT_BAND_10 = 1321.0789,
             K1_CONSTANT_BAND_11 = 480.8883, K2_CONSTANT_BAND_11 = 1201.1442)
    if(!missing(MTL) && length(MTL)){
      cal[] <- mtlValues(MTL, names(cal))
    } else if(sat=="L9"){
      stop("Landsat 9 needs MTL for the thermal calibration constants ",
           "(or use Level-2 surface temperature, surf.temp in METRIC.EB())")
    }
    B10 <- if(!missing(thermalband)) thermalband else image.DN$Thermal1
    L_b10 <- B10*cal[["RADIANCE_MULT_BAND_10"]]+cal[["RADIANCE_ADD_BAND_10"]]
    bT_b10 <- cal[["K2_CONSTANT_BAND_10"]] / log((cal[["K1_CONSTANT_BAND_10"]]/L_b10)+1)
    if(method == "SC"){
      #atm.coeff from Jimenez-Munoz, 2014
      atm.coeff <- matrix(data=c(0.04019, 0.02916, 1.01523,
                                 -0.38333, -1.50294, 0.20324,
                                 0.00918, 1.36072, -0.27514),
                          byrow= TRUE, nrow= 3, ncol=3)
      ea = (WeatherStation$RH/100)*0.6108*exp((17.27*WeatherStation$temp)/
                                                (WeatherStation$temp+237.3))
      atm.cond <- matrix(data=c(ea^2, ea, 1), nrow= 3, ncol=1) #not ea, vwc in g/cm2
      AF.j <- atm.coeff %*% atm.cond
      AF.m <- c(1/tau_NB, -R_sky - (Rp/tau_NB), R_sky)
      delta <- (bT_b10^2)/(1324*L_b10 )
      gama <- bT_b10 - (bT_b10^2/1324)
      Ts <- delta*((1/epsilon_NB)*(AF.m[1]*L_b10+AF.m[2])+AF.m[3])+gama
    }
    if(method=="SW"){
      L_b11 <- image.DN$Thermal2*cal[["RADIANCE_MULT_BAND_11"]]+cal[["RADIANCE_ADD_BAND_11"]]
      bT_b11 <- cal[["K2_CONSTANT_BAND_11"]] / log((cal[["K1_CONSTANT_BAND_11"]]/L_b11)+1)
      SW.c <- c(-0.268, 1.378, 0.183, 54.30, -2.238, -129.20, 16.40)
      w <- 1.2 #Atm water vapor content in g cm-2
      delta_e <- 0 # emissivity difference 
      Ts <- bT_b10 + SW.c[2] * (bT_b10 - bT_b11) +
        SW.c[3] * (bT_b10 - bT_b11)^2 + SW.c[1] +
        (SW.c[4] + SW.c[5] * w) * (1 - epsilon_NB) 
      # If I had two different emissivity, I should change epsilon_NB parameter
      # with the mean emissivity and add +(SW.c[6] + SW.c[7] * w) * delta_e
    }
  }
  if(sat=="L7"){
    if(!missing(thermalband)){
      B6 <- thermalband
    } else {
      B6 <- image.DN$Thermal1  
    }
    L_t_6 <-  0.067 * B6 - 0.06709
    L7_K1 <- 666.09 
    L7_K2 <- 1282.71 
    Rc <- ((L_t_6 - Rp) / tau_NB) - (1-epsilon_NB)*R_sky
    Ts <- L7_K2 / log((epsilon_NB*L7_K1/Rc)+1)}
  Ts <- saveLoadClean(imagestack = Ts, 
                      file = "Ts", overwrite=TRUE)
  return(Ts)
}

#' Calculates outgoing long-wave radiation
#' @description
#' This function estimates the outgoing long-wave radiation using the
#' Stefan-Boltzmann equation.
#' @param LAI  SpatRaster with leaf area index, from LAI()
#' @param Ts   land surface temperature (K), from surfaceTemperature()
#' @return SpatRaster with outgoing long-wave radiation (W/m2).
#' @author Guillermo Federico Olmedo
#' @author David Fonseca-Luengo
#' @family net radiation related functions
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for 
#' mapping evapotranspiration with internalized calibration (METRIC) - Model" 
#' Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
outLWradiation <- function(LAI, Ts){
  surf.emissivity <- 0.95 + 0.01 * LAI
  surf.emissivity[LAI>3] <- 0.98
  Rl.out <- surf.emissivity * 5.67e-8 * Ts^4
  Rl.out <- saveLoadClean(imagestack = Rl.out, 
                          file = "Rl.out", overwrite=TRUE)
  return(Rl.out)
}

#' Calculates incoming long-wave radiation
#' @description
#' This function estimates the incoming long-wave radiation using the
#' Stefan-Boltzmann equation. In addition, the empirical equation of
#' Bastiaanssen (1995) with coefficients developed by Allen (2000) is used to
#' estimate the effective atmospheric emissivity.
#' @param WeatherStation   waterWeatherStation object (see \code{\link{read.WSdata}}) or the
#' data.frame returned by getDataWS()
#' @param DEM              SpatRaster with a digital elevation model (m)
#' @param solar.angles     SpatRaster with layers latitude, declination,
#' hour.angle, incidence.hor and incidence.rel (radians), from solarAngles()
#' @param Ts               land surface temperature (K), from
#' surfaceTemperature()
#' @param Kt               clearness (turbidity) coefficient of the atmosphere,
#' 0 < Kt <= 1, used for the short-wave transmissivity: 1 (default) for clean
#' air. See \code{\link{incSWradiation}}
#' @return SpatRaster with incoming long-wave radiation (W/m2).
#' @author Guillermo Federico Olmedo
#' @author David Fonseca-Luengo
#' @family net radiation related functions
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#'
#' Bastiaanssen W. Regionalization of surface flux densities and moisture indicators in composite terrain: A remote sensing approach under clear skies in Mediterranean climates. Ph.D. dissertation, CIP Data Koninklijke Bibliotheek, Den Haag, The Netherlands, 1995, p273 \cr
#'
#' Allen R. RAPID long-wave radiation calculations and model comparisons Internal report, University of Idaho, Kimberly, Idaho, 2000 \cr
#' @export
# Add a function to get "cold" pixel temperature so in can be used in the next function
incLWradiation <- function(WeatherStation, DEM, solar.angles, Ts, Kt = 1){
  checkKt(Kt)
  if(inherits(WeatherStation, "waterWeatherStation")){
    WeatherStation <- getDataWS(WeatherStation)
  }
  ea <- (WeatherStation$RH/100)*0.6108*exp((17.27*WeatherStation$temp)/
                                             (WeatherStation$temp+237.3))
  tau.sw <- SWtrasmisivity(Kt = Kt, ea, DEM, solar.angles$incidence.hor)
  #temp <-  WeatherStation$temp+273.15 - (DEM - WeatherStation$elev) * 6.49 / 1000 
  ## Temperature in Kelvin
  #Mountain lapse effects from International Civil Aviation Organization
  ef.atm.emissivity  <- 0.85*(-1*log(tau.sw))^0.09
  Rl.in <- ef.atm.emissivity * 5.67e-8 * Ts^4
  #Rl.in <- ef.atm.emissivity * 5.67e-8 * temp^4
  Rl.in <- saveLoadClean(imagestack = Rl.in, 
                         file = "Rl.in", overwrite=TRUE)
  return(Rl.in)
}

#' Estimates net radiation
#' @description
#' This function estimates net radiation considering a surface radiation
#' balance. Equations use the information from the image, in addition to
#' measurements of actual vapor pressure and altitude.
#' @param LAI     SpatRaster with leaf area index, from LAI()
#' @param albedo  broadband surface albedo, from albedo()
#' @param Rs.inc  incoming short-wave radiation (W/m2), from incSWradiation()
#' @param Rl.inc  incoming long-wave radiation (W/m2), from incLWradiation()
#' @param Rl.out  outgoing long-wave radiation (W/m2), from outLWradiation()
#' @return SpatRaster with net radiation (W/m2).
#' @author Guillermo Federico Olmedo
#' @author David Fonseca-Luengo
#' @family net radiation related functions
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
netRadiation <- function(LAI, albedo, Rs.inc, Rl.inc, Rl.out){
  surf.emissivity <- 0.95 + 0.01 * LAI 
  surf.emissivity[LAI>3] <- 0.98
  Rn <- (1- albedo)*Rs.inc + Rl.inc - Rl.out - (1-surf.emissivity)*Rl.inc
  Rn <- saveLoadClean(imagestack = Rn, 
                      file = "Rn", overwrite=TRUE)
  return(Rn)
}


#' Estimates Soil Heat Flux
#' @description
#' This function implements models to estimate soil heat flux for different
#' surfaces and considering different inputs.
#' @param image    surface reflectance SpatRaster with layers named R and NIR,
#' used for NDVI (layers without standard names: the third and fourth, with a
#' warning)
#' @param Ts       land surface temperature (K), from surfaceTemperature()
#' @param albedo   broadband surface albedo, from albedo()
#' @param LAI      SpatRaster with leaf area index, from LAI()
#' @param Rn       net radiation (W/m2), from netRadiation()
#' @param aoi      Not used; kept for backward compatibility
#' @param method   method used for the G estimation: "Tasumi" (Tasumi, 2003;
#' default) or "Bastiaanssen" (Bastiaanssen, 2000)
#' @return SpatRaster with soil heat flux (W/m2).
#' @author Guillermo Federico Olmedo
#' @author David Fonseca-Luengo
#' @family net radiation related functions
#' @references
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#'
#' Tasumi, M. (2003). Progress in operational estimation of regional evapotranspiration using satellite imagery. Ph.D. dissertation, University of Idaho, Moscow, Idaho. \cr
#'
#' Bastiaanssen, W. G. M. (2000). SEBAL-based sensible and latent heat fluxes in the irrigated Gediz Basin, Turkey. Journal of Hydrology, 229(1-2), 87-100. https://doi.org/10.1016/S0022-1694(99)00202-4 \cr
#' @export
soilHeatFlux <- function(image, Ts, albedo, LAI, Rn, aoi, method = "Tasumi"){
  checkChoice(method, waterChoices$G.method, "method")
  sr.4.5 <- landsatBands(image, c("R", "NIR"), pos = 3:4)
  NDVI <- (sr.4.5[[2]] - sr.4.5[[1]])/(sr.4.5[[1]] + sr.4.5[[2]])
  if(method == "Bastiaanssen"){
    G <- ((Ts - 273.15)*(0.0038+0.0074*albedo)*(1-0.98*NDVI^4))*Rn
  }
  if(method == "Tasumi"){
    G <- rast(NDVI)
    e <- 2.71828
    G <- (0.05 + 0.18 * e^(-0.521*LAI)) * Rn
    G[LAI < 0.5] <- ((1.8*(Ts[LAI < 0.5] - 273.16) / Rn[LAI < 0.5]) + 0.084) * 
      Rn[LAI < 0.5]
  }
  G <- saveLoadClean(imagestack = G, file = "G", overwrite=TRUE)
  return(G)
}

