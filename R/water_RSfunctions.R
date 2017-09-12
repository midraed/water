#' Load Landsat data from folder
#' @description
#' This function loads Landsat bands from a specific folder. 
#' @param path  folder where band files are stored
#' @param sat   "L7" for Landsat 7, "L8" for Landsat 8, "MODIS" for MODIS or 
#' "auto" to guess from filenames
#' @param aoi   area of interest to crop images, if waterOptions("autoAoi") == 
#' TRUE will look for any object called aoi on .GlobalEnv
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @family remote sensing support functions
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
loadImage <-  function(path = getwd(), sat="auto", aoi){
  ## TODO: For L8 i should load the sirfrefl directly, like in MODIS!
  if(sat=="auto"){sat = getSat(path)} #DRY!
  if(sat=="L8"){bands <- c(2:7, 10, 11)}
  if(sat=="L7"){bands <- c(1:5,7, 6)}
  if(sat=="MODIS"){bands <- c(1:7)}
  ## Check for more than 1 image on the same folder
  if(sat=="L8" | sat=="L7"){
    image_list <- list.files(path=path, pattern = paste0("^L[EC]\\d+\\w+\\d+_(b|B|band)",
                                                         bands[1] ,".(TIF|tif)$"))
    if(length(image_list) > 1) {  ## Check if there are more images present on folder
      image_pattern <- substr(image_list[[1]], 0, nchar(image_list[[1]])-5)
      warning(paste("More than 1 image present on path. Using", 
                    substr(image_pattern, 0, nchar(image_pattern)-2)))
    } else {
      image_pattern <- substr(image_list[[1]], 0, nchar(image_list[[1]])-5)
    }
    bandnames <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2", "Thermal1")
    if(sat=="L8"){bandnames <- c(bandnames, "Thermal2")}
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
  for(i in 1:length(bands)){
    stack1[i] <- raster(list.files(path=path, 
                                   pattern = paste0(image_pattern, bands[i], "(_1)?", "(_VCID_1)?",
                                                    ".(TIF|tif)$"), full.names = T))
  }
  if(sat == "MODIS"){
    thermal <- list.files(path=path, pattern = paste0(".LST_Day_1km",
                                                      ".(TIF|tif)$"), full.names = T)[1]
    stack1[8] <- raster(thermal)
    time <- list.files(path=path, pattern = paste0(".Day_view_time",
                                                      ".(TIF|tif)$"), full.names = T)[1]
    stack1[9] <- raster(time)
  }
  raw.image <- do.call(stack, stack1)
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



#' Load Landsat 8 surface reflectance data from folder
#' @description
#' This function loads Landsat bands from a specific folder. 
#' @param path  folder where band files are stored
#' @param aoi   area of interest to crop images, if waterOptions("autoAoi") == 
#' TRUE will look for any object called aoi on .GlobalEnv
#' @author Guillermo Federico Olmedo
#' @family remote sensing support functions
#' @export
loadImageSR <-  function(path = getwd(),  aoi){
  files <- list.files(path = path, pattern = "_sr_band+[2-7].tif$", full.names = T)
  stack1 <- list()
  for(i in 1:6){
    stack1[i] <- raster(files[i])}
  image_SR <- do.call(stack, stack1)
  image_SR <- aoiCrop(image_SR, aoi) 
  image_SR <- image_SR / 10000
  bandnames <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2")
  image_SR <- saveLoadClean(imagestack = image_SR, 
                             stack.names = bandnames, 
                             file = "image_SR", 
                             overwrite=TRUE)
  return(image_SR)}  


#' Calculates radiance
#' @description
#' This function calculates radiance
#' @param image.DN      raw image in digital numbers
#' @param sat           "L7" for Landsat 7, "L8" for Landsat 8 or "auto" to guess from filenames 
#' @param MTL           Landsat Metadata File
#' @author Guillermo Federico Olmedo
#' @author María Victoria Munafó
#' @family remote sensing support functions
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' LPSO. (2004). Landsat 7 science data users handbook, Landsat Project Science Office, NASA Goddard Space Flight Center, Greenbelt, Md., (http://landsathandbook.gsfc.nasa.gov/) (Feb. 5, 2007) \cr
#' @export
calcRadiance <- function(image.DN, sat = "auto", MTL){
  path <- getwd()
  if(sat=="auto"){sat = getSat(path)} #DRY!
  if(sat=="L8"){bands <- c(2:7, 10, 11)}
  if(sat=="L7"){bands <- c(1:5,7, 6)}
  if(missing(MTL)){MTL <- list.files(path = getwd(), pattern = "MTL.txt", full.names = T)}
 
  MTL <- readLines(MTL, warn=FALSE)
  ADD <- vector()
  MULT <- vector()
 
   for( i in 1:length(bands)){
     
     ADDstring <- paste0("RADIANCE_ADD_BAND_", bands[i])
     ADDstring <- grep(ADDstring,MTL,value=TRUE)
     ADD[i] <- as.numeric(regmatches(ADDstring, 
                      regexec(text=ADDstring ,
                         pattern="([-]*)([0-9]{1,5})([.]+)([0-9]+)"))[[1]][1])
    
     MULTstring <- paste0("RADIANCE_MULT_BAND_",bands[i])
     MULTstring <- grep(MULTstring,MTL,value=TRUE)
     MULT[i] <- as.numeric(regmatches(MULTstring, 
                    regexec(text=MULTstring ,
                    pattern="([0-9]{1,5})([.]+)([0-9]+)(E-)([0-9]+)"))[[1]][1])
  
   }
  image <- image.DN * MULT + ADD
  bandnames <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2", "Thermal1")
  if(sat=="L8"){bandnames <- c(bandnames, "Thermal2")}
  image <- saveLoadClean(imagestack = image, 
                            stack.names = bandnames, 
                            file = "image_Rad", 
                            overwrite=TRUE)
  return(image)
}

#' Calculates Top of atmosphere reflectance
#' @description
#' This function calculates the TOA (Top Of Atmosphere) reflectance considering only the image metadata.
#' @param image.DN      raw image in digital numbers
#' @param sat           "L7" for Landsat 7, "L8" for Landsat 8 or "auto" to guess from filenames 
#' @param aoi           area of interest to crop images, if waterOptions("autoAoi") == TRUE will look for any object called aoi on .GlobalEnv
#' @param incidence.rel solar incidence angle, considering the relief
#' @param MTL           Landsat Metadata File
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David
#' @family remote sensing support functions
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#'
#' LPSO. (2004). Landsat 7 science data users handbook, Landsat Project Science Office, NASA Goddard Space Flight Center, Greenbelt, Md., (http://landsathandbook.gsfc.nasa.gov/) (Feb. 5, 2007) \cr
#' @export
calcTOAr <- function(image.DN, sat="auto", 
                     aoi, incidence.rel, MTL){
  path = getwd()
  if(sat=="auto"){sat = getSat(path)}
  if(sat=="L8"){bands <- 2:7}
  if(sat=="L7"){bands <- c(1:5,7)}
  if(sat=="L8"){
    image_TOA <- ((2.0000E-05 * image.DN[[1:6]]) + -0.100000) / incidence.rel  ### There is a small difference with ESPA TOA of less than 0.02
    names(image_TOA) <- names(image.DN[[1:6]])
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
#' @param aoi             area of interest to crop images, if waterOptions("autoAoi") == TRUE will look for any object called aoi on .GlobalEnv
#' @param incidence.hor   solar incidence angle, considering plain surface
#' @param WeatherStation  Weather Station data
#' @param surface.model   rasterStack with DEM, Slope and Aspect. See surface.model()
#' @author Guillermo Federico Olmedo
#' @author Fonseca-Luengo, David 
#' @family remote sensing support functions
#' @references 
#' Tasumi M.; Allen R.G. and Trezza, R. At-surface albedo from Landsat and MODIS satellites for use in energy balance studies of evapotranspiration Journal of Hydrolog. Eng., 2008, 13, (51-63) \cr
#'
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' @export
# incidence hor from TML?? 
calcSR <- function(image.TOAr, sat="auto", aoi, incidence.hor, 
                   WeatherStation, surface.model){
  if(class(WeatherStation)== "waterWeatherStation"){
    WeatherStation <- getDataWS(WeatherStation)
  }
  path <- getwd()
  if(sat=="auto"){sat = getSat(path)}
  if(sat=="L8"){stop("water package does not include a model to calculate surface reflectance 
  for Landsat 8 images. Landsat 8 users should download precalculated surface reflectances from 
  espa website (espa.cr.usgs.gov). ")}
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
