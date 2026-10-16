#' Load Landsat data from a folder
#' @description
#' Loads Landsat 7, 8 or 9 Level-1 bands (digital numbers), or MODIS bands,
#' from a folder. Landsat file names can follow the Collection 1/2 naming
#' (e.g. \code{LC08_L1TP_232083_20160209_20200907_02_T1_B2.TIF}, any tier) or
#' the pre-Collection naming (e.g. \code{LC82320832016040LGN00_band2.tif});
#' Level-2 surface reflectance files are skipped (use loadImageSR()). If the
#' folder holds more than one image, the first one is used, with a warning.
#' @param path  folder where band files are stored
#' @param sat   "L7" for Landsat 7, "L8" for Landsat 8, "L9" for Landsat 9,
#' "MODIS" for MODIS or "auto" to detect the Landsat satellite from the
#' SPACECRAFT_ID of an MTL file in path or, failing that, from the Landsat file
#' names
#' @param aoi   area of interest (e.g. from createAoi()) to crop the image. If
#' missing and waterOptions("autoAoi") is TRUE, an object called aoi in
#' .GlobalEnv is used if it exists
#' @return SpatRaster in digital numbers. Landsat 8/9: layers B, G, R, NIR,
#' SWIR1, SWIR2, Thermal1 and Thermal2 (bands 2 to 7, 10 and 11). Landsat 7: B,
#' G, R, NIR, SWIR1, SWIR2 and Thermal1 (bands 1 to 5, 7 and 6 low gain).
#' MODIS: R, NIR, B, G, SWIR1, SWIR2, SWIR3 (reflectance), LST (K) and Time
#' (hours).
#' @author Guillermo Federico Olmedo
#' @author David Fonseca-Luengo
#' @family remote sensing support functions
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
loadImage <-  function(path = getwd(), sat="auto", aoi){
  ## TODO: For L8 i should load the sirfrefl directly, like in MODIS!
  if(sat=="auto"){sat = getSat(path)} #DRY!
  if(sat %in% c("L8", "L9")){bands <- c(2:7, 10, 11)}
  if(sat=="L7"){bands <- c(1:5,7, 6)}
  if(sat=="MODIS"){bands <- c(1:7)}
  ## Check for more than 1 image on the same folder
  if(sat %in% c("L7", "L8", "L9")){
    # Level-1 band files: <scene>_B<n>.TIF (Collection 1/2, any tier) or
    # <scene>_band<n>.tif (pre-collection); Level-2 SR/ST files are skipped
    image_list <- list.files(path=path, pattern = paste0("^L[CEOT]0?[789].*_(b|B|band)",
                                                         bands[1], "\\.(TIF|tif)$"))
    image_list <- image_list[!grepl("_(SR|ST)_B", image_list)]
    if(length(image_list) == 0){
      stop("No Landsat Level-1 band files found in ", path, ". Expected names ",
           "like LC08_L1TP_232083_20160209_20200907_02_T1_B2.TIF or ",
           "LC82320832016040LGN00_band2.tif. For Landsat 8/9 Collection 2 ",
           "Level-2 surface reflectance use loadImageSR()")
    }
    # scene prefix, including the "B"/"band" before the band number
    image_pattern <- sub(paste0(bands[1], "\\.(TIF|tif)$"), "", image_list)
    if(length(image_list) > 1) {  ## Check if there are more images present on folder
      warning(paste("More than 1 image present on path. Using",
                    sub("_(b|B|band)$", "", image_pattern[1])))
    }
    image_pattern <- image_pattern[1]
    bandnames <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2", "Thermal1")
    if(sat %in% c("L8", "L9")){bandnames <- c(bandnames, "Thermal2")}
  }
  if(sat=="MODIS"){
    image_list <- list.files(path=path, pattern = paste0(".sur_refl_b0",
                                                         bands[1] ,"_1.(TIF|tif)$"))
    if(length(image_list) > 1) { ## Check if there are more images present on folder
      image_pattern <- substr(image_list[[1]], 0, nchar(image_list[[1]])-7)
      warning(paste("More than 1 image present on path. Using", 
                    substr(image_pattern, 0, nchar(image_pattern)-2)))
    } else {
      image_pattern <- substr(image_list[[1]], 0, nchar(image_list[[1]])-7)
    }
    bandnames <- c("R", "NIR", "B", "G", "SWIR1", "SWIR2", "SWIR3", "LST", "Time") # band names for MOD09GA
  }
  
  stack1 <- list()
  files <- list.files(path = path)
  for(i in 1:length(bands)){
    band.file <- files[startsWith(files, image_pattern) &
                         grepl(paste0("^", bands[i], "(_1)?(_VCID_1)?\\.(TIF|tif)$"),
                               substring(files, nchar(image_pattern) + 1))]
    if(length(band.file) == 0){
      stop("Band ", bands[i], " not found for ", image_pattern, " in ", path)
    }
    stack1[[i]] <- rast(file.path(path, band.file[1]))
  }
  if(sat == "MODIS"){
    thermal <- list.files(path=path, pattern = paste0(".LST_Day_1km",
                                                      ".(TIF|tif)$"), full.names = T)[1]
    stack1[[8]] <- rast(thermal)
    time <- list.files(path=path, pattern = paste0(".Day_view_time",
                                                      ".(TIF|tif)$"), full.names = T)[1]
    stack1[[9]] <- rast(time)
  }
  raw.image <- do.call(c, stack1)
  raw.image <- aoiCrop(raw.image, aoi)
  if(sat=="MODIS"){for(i in 1:7){
    raw.image[[i]] <- raw.image[[i]]*0.0001
  }
    raw.image[[8]] <- raw.image[[8]]*0.02
    raw.image[[9]] <- raw.image[[9]]*0.1
  }
  raw.image <- saveLoadClean(imagestack = raw.image, 
                             stack.names = bandnames, 
                             file = "imageDN", 
                             overwrite=TRUE)
  return(raw.image) 
}  



#' Load Landsat 8/9 surface reflectance (and surface temperature) from a folder
#' @description
#' Loads Landsat 8 or 9 Collection 2 Level-2 surface reflectance bands 2 to 7
#' (\code{*_SR_B2.TIF} to \code{*_SR_B7.TIF}) from a folder, scaled as
#' DN * 0.0000275 - 0.2, with the fill value 0 set to NA. All six bands are
#' required, as other functions use them by position.
#' If \code{*_ST_B10.TIF} is present it is added as the last layer, "Band_10",
#' in Kelvin (DN * 0.00341802 + 149, fill value 0 set to NA).
#' It can be passed as surf.temp to METRIC.EB(), METRIC.Rn() or METRIC.G().
#' Landsat 7 Level-2 products are not supported (different band numbers).
#' @param path  folder where band files are stored
#' @param aoi   area of interest (e.g. from createAoi()) to crop the image. If
#' missing and waterOptions("autoAoi") is TRUE, an object called aoi in
#' .GlobalEnv is used if it exists
#' @return SpatRaster with surface reflectance (unitless) in layers B, G, R,
#' NIR, SWIR1 and SWIR2 (Landsat 8/9 bands 2 to 7). If \code{*_ST_B10.TIF} is
#' present, a seventh layer Band_10 with surface temperature (K) is added.
#' @author Guillermo Federico Olmedo
#' @family remote sensing support functions
#' @export
#' @examples
#' scene_folder <- system.file("extdata", "L8_C2", package = "water")
#' image.SR <- loadImageSR(path = scene_folder)
#' names(image.SR)
#' plot(image.SR[["NIR"]])
loadImageSR <-  function(path = getwd(),  aoi){
  files <- list.files(path = path, pattern = "_SR_B[2-7]\\.(TIF|tif)$", full.names = T)
  if(!length(files)){
    stop("No Collection 2 Level-2 surface reflectance bands (*_SR_B2.TIF to ",
         "*_SR_B7.TIF) found in ", path)
  }
  band <- sub(".*_SR_B([2-7])\\.(TIF|tif)$", "\\1", basename(files))
  if(anyDuplicated(band)) stop("More than 1 image present on path: ", path)
  if(!setequal(band, as.character(2:7))){
    stop("Missing surface reflectance bands in ", path, ": ",
         paste(setdiff(as.character(2:7), band), collapse = ", "))
  }
  bandnames <- c("2" = "B", "3" = "G", "4" = "R", "5" = "NIR", "6" = "SWIR1", "7" = "SWIR2")
  ord <- order(band)
  image_SR <- rast(files[ord])
  image_SR <- aoiCrop(image_SR, aoi)
  bandnames <- unname(bandnames[band[ord]])
  image_SR[image_SR == 0] <- NA  # Collection 2 fill value
  image_SR <- image_SR * 0.0000275 - 0.2
  # L8/L9 surface temperature, appended last so SR band positions don't change
  st <- list.files(path = path, pattern = "^LC0[89]_.*_ST_B10\\.(TIF|tif)$", full.names = T)
  if(length(st)){
    st <- aoiCrop(rast(st[1]), aoi)
    st[st == 0] <- NA
    image_SR <- c(image_SR, st * 0.00341802 + 149)
    bandnames <- c(bandnames, "Band_10")
  }
  image_SR <- saveLoadClean(imagestack = image_SR,
                             stack.names = bandnames, 
                             file = "image_SR", 
                             overwrite=TRUE)
  return(image_SR)}  


#' Calculates at-sensor radiance
#' @description
#' Converts digital numbers to at-sensor spectral radiance using the
#' RADIANCE_MULT_BAND_x and RADIANCE_ADD_BAND_x values of the MTL file.
#' @param image.DN      raw image in digital numbers, e.g. from loadImage().
#' Bands are selected by name (B, G, R, NIR, SWIR1, SWIR2, Thermal1, Thermal2)
#' when present
#' @param sat           "L7", "L8", "L9" or "auto" to detect it from the MTL
#' file (or, if MTL is missing, from an MTL file or Landsat file names in the
#' working directory)
#' @param MTL           path to the Landsat metadata (MTL) file. If missing,
#' the only *MTL.txt file in the working directory is used
#' @return SpatRaster with at-sensor spectral radiance (W m-2 sr-1 um-1) in
#' layers B, G, R, NIR, SWIR1, SWIR2, Thermal1 and, for Landsat 8/9, Thermal2.
#' @author Guillermo Federico Olmedo
#' @author María Victoria Munafó
#' @family remote sensing support functions
#' @references
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#'
#' LPSO. (2004). Landsat 7 science data users handbook, Landsat Project Science Office, NASA Goddard Space Flight Center, Greenbelt, Md., (http://landsathandbook.gsfc.nasa.gov/) (Feb. 5, 2007) \cr
#'
#' U.S. Geological Survey (2019). Landsat 8 (L8) Data Users Handbook, LSDS-1574. EROS, Sioux Falls, South Dakota. https://www.usgs.gov/landsat-missions/landsat-8-data-users-handbook \cr
#' @export
calcRadiance <- function(image.DN, sat = "auto", MTL){
  path <- getwd()
  if(sat=="auto"){sat = getSat(path, MTL = if(!missing(MTL)) MTL)}
  if(sat %in% c("L8", "L9")){bands <- c(2:7, 10, 11)}
  if(sat=="L7"){bands <- c(1:5,7, 6)}
  MTL <- findMTL(if(!missing(MTL)) MTL)
  keys <- bands
  keys[keys == 6] <- "6_VCID_1"  # Landsat 7 band 6, low gain
  ADD <- mtlValues(MTL, paste0("RADIANCE_ADD_BAND_", keys))
  MULT <- mtlValues(MTL, paste0("RADIANCE_MULT_BAND_", keys))
  bandnames <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2", "Thermal1")
  if(sat %in% c("L8", "L9")){bandnames <- c(bandnames, "Thermal2")}
  image <- landsatBands(image.DN, bandnames) * MULT + ADD
  image <- saveLoadClean(imagestack = image, 
                            stack.names = bandnames, 
                            file = "image_Rad", 
                            overwrite=TRUE)
  return(image)
}

#' Calculates top-of-atmosphere reflectance
#' @description
#' Calculates the TOA (top-of-atmosphere) reflectance, divided by the cosine of
#' the solar incidence angle of each pixel. For Landsat 8 and 9 the reflectance
#' rescaling coefficients of the Level-1 products (2e-05 and -0.1) are used; for
#' Landsat 7, the gains, biases and ESUN values of the Landsat 7 handbook and
#' the acquisition date in the MTL file.
#' @param image.DN      raw image in digital numbers, e.g. from loadImage().
#' Bands are selected by name (B, G, R, NIR, SWIR1, SWIR2) when present
#' @param sat           "L7", "L8", "L9" or "auto" to detect it from the MTL
#' file (or, if MTL is missing, from an MTL file or Landsat file names in the
#' working directory)
#' @param aoi           Not used; kept for backward compatibility. Crop the
#' inputs beforehand (e.g. loadImage(..., aoi = aoi))
#' @param incidence.rel solar incidence angle (radians) considering the relief,
#' from solarAngles()
#' @param MTL           path to the Landsat metadata (MTL) file, needed for
#' Landsat 7 and for sat = "auto". If missing, the only *MTL.txt file in the
#' working directory is used
#' @return SpatRaster with top-of-atmosphere reflectance (unitless) in layers B,
#' G, R, NIR, SWIR1 and SWIR2.
#' @author Guillermo Federico Olmedo
#' @author David Fonseca-Luengo
#' @family remote sensing support functions
#' @references
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#'
#' LPSO. (2004). Landsat 7 science data users handbook, Landsat Project Science Office, NASA Goddard Space Flight Center, Greenbelt, Md., (http://landsathandbook.gsfc.nasa.gov/) (Feb. 5, 2007) \cr
#'
#' U.S. Geological Survey (2019). Landsat 8 (L8) Data Users Handbook, LSDS-1574. EROS, Sioux Falls, South Dakota. https://www.usgs.gov/landsat-missions/landsat-8-data-users-handbook \cr
#' @export
calcTOAr <- function(image.DN, sat="auto", 
                     aoi, incidence.rel, MTL){
  path = getwd()
  if(sat=="auto"){sat = getSat(path, MTL = if(!missing(MTL)) MTL)}
  if(sat %in% c("L8", "L9")){bands <- 2:7}
  if(sat=="L7"){bands <- c(1:5,7)}
  if(sat %in% c("L8", "L9")){
    refl <- landsatBands(image.DN, c("B", "G", "R", "NIR", "SWIR1", "SWIR2"))
    image_TOA <- ((2.0000E-05 * refl) + -0.100000) / cos(incidence.rel)
  }
  ### Ro TOA L7
  if(sat=="L7"){
    MTL <- readLines(findMTL(if(!missing(MTL)) MTL, path = path), warn=FALSE)
    ESUN <- c(1997, 1812, 1533, 1039, 230.8, 84.90) # Landsat 7 Handbook
    ## On sect 11.3, L7 handbook recommends using this formula for Ro TOA
    ## O using DN - QCALMIN with METRIC 2010 formula.
    Gain <- c(1.181, 1.210, 0.943, 0.969, 0.191, 0.066)
    Bias <- c(-7.38071, -7.60984, -5.94252, -6.06929, -1.19122, -0.41650)
    if(missing(image.DN)){image.DN <- loadImage(path = path)}
    time.line <- grep("SCENE_CENTER_TIME",MTL,value=TRUE)
    date.line <- grep("DATE_ACQUIRED",MTL,value=TRUE)
    sat.time <-regmatches(time.line,regexec(text=time.line,
                                            pattern="([0-9]{2})(:)([0-9]{2})(:)([0-9]{2})(.)([0-9]{2})"))[[1]][1]
    sat.date <-regmatches(date.line,regexec(text=date.line,
                                            pattern="([0-9]{4})(-)([0-9]{2})(-)([0-9]{2})"))[[1]][1]
    sat.datetime <- strptime(paste(sat.date, sat.time), 
                             format = "%Y-%m-%d %H:%M:%S", tz="GMT")
    DOY <-  sat.datetime$yday +1
    d2 <- 1/(1+0.033*cos(DOY * 2 * pi/365))
    dr <- 1 + 0.033 * cos(DOY * (2 * pi / 365))
    refl <- landsatBands(image.DN, c("B", "G", "R", "NIR", "SWIR1", "SWIR2"))
    Ro.TOAr <- list()
    for(i in 1:6){
      Ro.TOAr[[i]] <- (pi * (Gain[i] * refl[[i]] + Bias[i])) / (ESUN[i] * cos(incidence.rel) * dr)
    }
    image_TOA <- do.call(c, Ro.TOAr)
  }
  #### 
  image_TOA <- saveLoadClean(imagestack = image_TOA, 
                             stack.names = c("B", "G", "R", "NIR", "SWIR1", "SWIR2"), 
                             file = "image_TOAr", 
                             overwrite=TRUE)
  return(image_TOA)
}  

#' Calculates surface reflectance for Landsat 7
#' @description
#' Calculates surface reflectance from top-of-atmosphere reflectance using the
#' model developed by Tasumi et al. (2008) and Allen et al. (2007), applied band
#' by band. Only for Landsat 7: for Landsat 8 and 9 use the Collection 2
#' Level-2 surface reflectance (see \code{\link{loadImageSR}}).
#' @param image.TOAr      SpatRaster with top-of-atmosphere reflectance (B, G,
#' R, NIR, SWIR1, SWIR2), e.g. from calcTOAr()
#' @param sat             Only "L7" is supported; "L8" and "L9" stop with an
#' error. "auto" detects the satellite from an MTL file or Landsat file names
#' in the working directory
#' @param aoi             Not used; kept for backward compatibility
#' @param incidence.hor   solar incidence angle (radians) on a horizontal
#' surface, from solarAngles()
#' @param WeatherStation  waterWeatherStation object (see \code{\link{read.WSdata}}) or the
#' data.frame returned by getDataWS()
#' @param surface.model   SpatRaster with layers DEM, Slope and Aspect, from
#' METRICtopo()
#' @param MTL             Landsat metadata file, used to detect the satellite
#' when sat = "auto"
#' @return SpatRaster with surface reflectance (unitless) in layers B, G, R,
#' NIR, SWIR1 and SWIR2.
#' @author Guillermo Federico Olmedo
#' @author David Fonseca-Luengo
#' @family remote sensing support functions
#' @references
#' Tasumi, M., Allen, R. G., & Trezza, R. (2008). At-surface reflectance and albedo from satellite for operational calculation of land surface energy balance. Journal of Hydrologic Engineering, 13(2), 51-63. https://doi.org/10.1061/(ASCE)1084-0699(2008)13:2(51) \cr
#'
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
# incidence hor from TML?? 
calcSR <- function(image.TOAr, sat="auto", aoi, incidence.hor,
                   WeatherStation, surface.model, MTL){
  if(inherits(WeatherStation, "waterWeatherStation")){
    WeatherStation <- getDataWS(WeatherStation)
  }
  path <- getwd()
  if(sat=="auto"){sat = getSat(path, MTL = if(!missing(MTL)) MTL)}
  if(sat %in% c("L8", "L9")){stop("water package does not include a model to calculate ",
    "surface reflectance for Landsat 8/9. Download Collection 2 Level-2 products ",
    "and load them with loadImageSR().")}
  if(sat=="L7"){
    bands <- c(1:5,7)
    if(missing(image.TOAr)){image.TOAr <- calcTOAr()}
    P <- 101.3*((293-0.0065 * surface.model$DEM)/293)^5.26
    ea <- (WeatherStation$RH/100)*0.6108*exp((17.27*WeatherStation$temp)/(WeatherStation$temp+237.3))
    W <- 0.14 * ea * P + 2.1
    Kt <- 1
    Cnb <- matrix(data=c(0.987, 2.319, 0.951, 0.375, 0.234, 0.365,
                         -0.00071, -0.00016, -0.00033, -0.00048, -0.00101, -0.00097,
                         0.000036, 0.000105, 0.00028, 0.005018, 0.004336, 0.004296,
                         0.0880, 0.0437, 0.0875, 0.1355, 0.056, 0.0155,
                         0.0789, -1.2697, 0.1014, 0.6621, 0.7757, 0.639,
                         0.64, 0.31, 0.286, 0.189, 0.274, -0.186), byrow = T, nrow=6, ncol=6)
    tau_in <- list()
    tau_out <- list()
    for(i in 1:6){
      tau_in[[i]] <- Cnb[1,i] * exp((Cnb[2,i]*P/(Kt*cos(incidence.hor)))-
                                    ((Cnb[3,i]*W+Cnb[4,i])/cos(incidence.hor)))+Cnb[5,i]
    }
    eta = 0 # eta it's the satellite nadir angle
    for(i in 1:6){
      tau_out[[i]] <- Cnb[1,i] * exp((Cnb[2,i]*P/(Kt*cos(eta)))-
                                     ((Cnb[3,i]*W+Cnb[4,i])/cos(eta)))+Cnb[5,i]
    }
    path_refl <- list()
    for(i in 1:6){
      path_refl[[i]] <- Cnb[6,i] * (1 - tau_in[[i]])
    }
    stack_SR <- list()
    for(i in 1:6){
      stack_SR[[i]] <- (image.TOAr[[i]] - path_refl[[i]]) / (tau_in[[i]] * tau_out[[i]])
    }
    image_SR <- do.call(c, stack_SR)
  }
  image_SR <- saveLoadClean(imagestack = image_SR, 
                            stack.names = c("B", "G", "R", "NIR", "SWIR1", "SWIR2"), 
                            file = "image_SR", 
                            overwrite=TRUE)
  return(image_SR)
}  


#' Mask clouds using a quality band
#' @description
#' Masks clouds, cloud shadows and fill pixels of an image using its quality
#' band. Use it after loading the image (e.g. with loadImageSR()) and before
#' METRIC.EB(), so that cloudy pixels are not used as anchor pixels. All layers
#' of the image are masked, including the surface temperature (Band_10).
#' @param path     folder where the quality band is searched when cfmask is
#' missing: a Collection 2 \code{*_QA_PIXEL.TIF} file or, for older products,
#' a \code{*_cfmask.tif} or \code{*_BQA.TIF} file
#' @param image    SpatRaster to mask, e.g. from loadImageSR() or loadImage()
#' @param cfmask   SpatRaster with the quality band. If missing, it is read
#' from path
#' @param keep     For older quality bands only: values to keep (all other
#' pixels are masked). Use 0 for cfmask files, 2720 for Landsat 8 Collection 1
#' BQA files and 672 for Landsat 7 Collection 1 BQA files. If missing, a
#' Collection 2 QA_PIXEL band is decoded by bits: pixels flagged as fill,
#' dilated cloud, cirrus, cloud or cloud shadow (bits 0 to 4) are masked.
#' @param buffer   distance (in map units, usually meters) around clouds and
#' cloud shadows that is also masked. Use 0 for no buffer
#' @return The image (SpatRaster) with the masked pixels set to NA.
#' @author Guillermo Federico Olmedo
#' @family remote sensing support functions
#' @export
cfmask <-  function(path = getwd(), image, cfmask, keep, buffer = 60){
  if(missing(cfmask)){
    file <- list.files(path = path, pattern = "_QA_PIXEL\\.(TIF|tif)$", full.names = TRUE)
    if(length(file) == 0){
      file <- list.files(path = path, pattern = "_cfmask\\.tif$|_(bqa|BQA)\\.(TIF|tif)$",
                         full.names = TRUE)
    }
    if(length(file) != 1){
      stop("Expected one quality band (*_QA_PIXEL.TIF, *_cfmask.tif or ",
           "*_BQA.TIF) in ", path, ", found ", length(file))
    }
    cfmask <- rast(file)
  }
  cfmask <- crop(cfmask, image)
  if(missing(keep)){
    if(!any(grepl("QA_PIXEL", c(names(cfmask), sources(cfmask)), ignore.case = TRUE))){
      keep <- 0  # pre-collection cfmask files
    }
  }
  if(missing(keep)){
    # Collection 2 QA_PIXEL: bit 0 fill; bits 1 to 4 dilated cloud, cirrus,
    # cloud and cloud shadow
    fill <- (cfmask %% 2) == 1
    cloudy <- ((cfmask %/% 2) %% 16) != 0
  } else {
    fill <- is.na(cfmask)
    cloudy <- !(cfmask %in% keep)
  }
  if(buffer > 0){
    cloudy <- terra::buffer(ifel(cloudy, 1, NA), width = buffer)
  }
  clear <- ifel(fill | cloudy, NA, 1)
  if(!isTRUE(global(clear, "max", na.rm = TRUE)[1, 1] == 1)){
    warning("cfmask() masked every pixel: check keep (for a Collection 2 ",
            "QA_PIXEL band leave keep missing; its layer or file name must ",
            "contain QA_PIXEL)")
  }
  image <- mask(image, clear)
  return(image)
}



# cfmask <-  function(path = getwd(), image, cfmask,
#                     buffer = 60){
#   if(missing(cfmask)){
#     file <- list.files(pattern="_sr_cloud.tif$")
#     cfmask <- rast(file)
#   }
#   cfmask <- crop(cfmask, image)
#   values(cfmask)[values(cfmask) <= 16] <- NA
#   values(cfmask)[!is.na(values(cfmask))] <- 1
#   image <- mask(image, cfmask)
#   return(image)
# }