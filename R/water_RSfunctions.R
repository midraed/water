#' Load Landsat data from folder
#' @description
#' This function loads Landsat bands from a specific folder. 
#' @param path  folder where band files are stored
#' @param sat   "L7" for Landsat 7, "L8" for Landsat 8 or "auto" to guess from filenames
#' @param aoi   area of interest to crop images, if waterOptions("autoAoi") == 
#' TRUE will look for any object called aoi on .GlobalEnv
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
loadImage <-  function(path = getwd(), sat="auto", aoi){
  if(sat=="auto"){sat = getSat(path)} #DRY!
  if(sat=="L8"){bands <- c(2:7, 10, 11)}
  if(sat=="L7"){bands <- c(1:5,7, 6)}
  stack1 <- list()
  
  ## Check for more than 1 image on the same folder
  image_list <- list.files(path=path, pattern = paste0("^L[EC]\\d+\\w+\\d+_(B|band)",
                                         bands[1] ,".(TIF|tif)$"))
  if(length(image_list) > 1) {
    image_pattern <- substr(image_list[[1]], 0, nchar(image_list[[1]])-5)
    warning(paste("More than 1 image present on path. Using", 
                  substr(image_pattern, 0, nchar(image_pattern)-2)))
  } else {
    image_pattern <- substr(image_list[[1]], 0, nchar(image_list[[1]])-5)
  }
  
  for(i in 1:length(bands)){
    stack1[i] <- raster(list.files(path=path, 
                                   pattern = paste0(image_pattern, bands[i], "(_VCID_1)?",
                                                    ".(TIF|tif)$"), full.names = T))
  }
  bandnames <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2", "Thermal1")
  if(sat=="L8"){bandnames <- c(bandnames, "Thermal2")}
  raw.image <- do.call(stack, stack1)
  raw.image <- aoiCrop(raw.image, aoi)                               
  raw.image <- saveLoadClean(imagestack = raw.image, 
                             stack.names = bandnames, 
                             file = "imageDN", 
                             overwrite=TRUE)
  return(raw.image) 
}  

#' Calculates Top of atmosphere reflectance
#' @description
#' This function calculates the TOA (Top Of Atmosphere) reflectance considering only the image metadata.
#' @param image.DN      raw image in digital numbers
#' @param sat           "L7" for Landsat 7, "L8" for Landsat 8 or "auto" to guess from filenames 
#' @param ESPA          Logical. If TRUE will look for espa.usgs.gov related products on working folder
#' @param aoi           area of interest to crop images, if waterOptions("autoAoi") == TRUE will look for any object called aoi on .GlobalEnv
#' @param incidence.rel solar incidence angle, considering the relief
#' @param MTL           Landsat Metadata File
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#'
#' LPSO. (2004). Landsat 7 science data users handbook, Landsat Project Science Office, NASA Goddard Space Flight Center, Greenbelt, Md., (http://landsathandbook.gsfc.nasa.gov/) (Feb. 5, 2007) \cr
#' @export
calcTOAr <- function(image.DN, sat="auto", 
                     ESPA=FALSE, aoi, incidence.rel, MTL){
  path = getwd()
  if(sat=="auto"){sat = getSat(path)}
  if(sat=="L8"){bands <- 2:7}
  if(sat=="L7"){bands <- c(1:5,7)}
  if(ESPA==TRUE & sat=="L8"){
    files <- list.files(path = path, pattern = "_toa_band+[2-7].tif$", full.names = T)
    stack1 <- list()
    for(i in 1:6){
      stack1[i] <- raster(files[i])
    }
    image_TOA <- do.call(stack, stack1)
    image_TOA <- aoiCrop(image_TOA, aoi)
    image_TOA <- image_TOA / 10000
  }
  ### Ro TOA L7
  if(sat=="L7"){
    if(missing(MTL)){MTL <- list.files(path = path, pattern = "MTL.txt", full.names = T)}
    MTL <- readLines(MTL, warn=FALSE)
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
    Ro.TOAr <- list()
    for(i in 1:6){
      Ro.TOAr[i] <- (pi * (Gain[i] * image.DN[[i]] + Bias[i])) / (ESUN[i] * cos(incidence.rel) * dr)
    }
    image_TOA <- do.call(stack, Ro.TOAr)
  }
  #### 
  image_TOA <- saveLoadClean(imagestack = image_TOA, 
                             stack.names = c("B", "G", "R", "NIR", "SWIR1", "SWIR2"), 
                             file = "image_TOAr", 
                             overwrite=TRUE)
  return(image_TOA)
}  

#' Calculates surface reflectance for L7
#' @description
#' Calculates surface reflectance from top of atmosphere radiance using the model developed by Tasumi et al. (2008) and Allen et al. (2007), which considers a band-by-band basis.
#' @param image.TOAr      raster stack. top of atmosphere reflectance image
#' @param sat             "L7" for Landsat 7, "L8" for Landsat 8 or "auto" to guess from filenames 
#' @param ESPA            Logical. If TRUE will look for espa.usgs.gov related products on working folder
#' @param aoi             area of interest to crop images, if waterOptions("autoAoi") == TRUE will look for any object called aoi on .GlobalEnv
#' @param incidence.hor   solar incidence angle, considering plain surface
#' @param WeatherStation  Weather Station data
#' @param surface.model   rasterStack with DEM, Slope and Aspect. See surface.model()
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David 
#' @references 
#' Tasumi M.; Allen R.G. and Trezza, R. At-surface albedo from Landsat and MODIS satellites for use in energy balance studies of evapotranspiration Journal of Hydrolog. Eng., 2008, 13, (51-63) \cr
#'
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
# incidence hor from TML?? 
calcSR <- function(image.TOAr, sat="auto", ESPA=FALSE, aoi, incidence.hor, 
                   WeatherStation, surface.model){
  if(class(WeatherStation)== "waterWeatherStation"){
    WeatherStation <- getDataWS(WeatherStation)
  }
  path <- getwd()
  if(sat=="auto"){sat = getSat(path)}
  if(sat=="L8"){bands <- 2:7}
  if(sat=="L7"){bands <- c(1:5,7)}
  if(ESPA==TRUE & sat=="L8"){
    files <- list.files(path = path, pattern = "_sr_band+[2-7].tif$", full.names = T)
    stack1 <- list()
    for(i in 1:6){
      stack1[i] <- raster(files[i])
    }
    image_SR <- do.call(stack, stack1)
    image_SR <- aoiCrop(image_SR, aoi) 
    image_SR <- image_SR / 10000
  }
  if(sat=="L7"){
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
      tau_in[i] <- Cnb[1,i] * exp((Cnb[2,i]*P/(Kt*cos(incidence.hor)))-
                                    ((Cnb[3,i]*W+Cnb[4,i])/cos(incidence.hor)))+Cnb[5,i]
    }
    eta = 0 # eta it's the satellite nadir angle
    for(i in 1:6){
      tau_out[i] <- Cnb[1,i] * exp((Cnb[2,i]*P/(Kt*cos(eta)))-
                                     ((Cnb[3,i]*W+Cnb[4,i])/cos(eta)))+Cnb[5,i]
    }
    path_refl <- list()
    for(i in 1:6){
      path_refl[i] <- Cnb[6,i] * (1 - tau_in[[i]])
    }
    stack_SR <- list()
    for(i in 1:6){
      stack_SR[i] <- (image.TOAr[[i]] - path_refl[[i]]) / (tau_in[[i]] * tau_out[[i]])
    }
    image_SR <- do.call(stack, stack_SR)
  }
  image_SR <- saveLoadClean(imagestack = image_SR, 
                            stack.names = c("B", "G", "R", "NIR", "SWIR1", "SWIR2"), 
                            file = "image_SR", 
                            overwrite=TRUE)
  return(image_SR)
}  
