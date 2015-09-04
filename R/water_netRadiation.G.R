
#' Create aoi polygon from topleft and bottomright coordinates
#' @param topleft a vector with topleft x,y coordinates 
#' @param bottomright a vector with bottomright x,y coordinates
#' @param EPSG Coordinate reference system EPSG code
#' @return object of class SpatialPolygons
#' @author Guillermo F Olmedo
#' @examples 
#' tl <- c(493300, -3592700)
#' br <- c(557200, -3700000) 
#' aoi <- createAoi(topleft = tl, bottomright=br, EPSG=32619)
#' plot(aoi)
#' @import raster sp proj4 
#' @export
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
# Maybe i can provide some points and use CHull like in QGIS-Geostat
createAoi <- function(topleft = c(x, y), bottomright= c(x, y), EPSG){
  aoi <- SpatialPolygons(
    list(Polygons(list(Polygon(coords = matrix(
      c(topleft[1],bottomright[1], bottomright[1],topleft[1],topleft[1],
        topleft[2], topleft[2], bottomright[2], 
        bottomright[2],topleft[2]), ncol=2, nrow= 5))), ID=1)))
  if(!missing(EPSG)){aoi@proj4string <- CRS(paste0("+init=epsg:", EPSG))}
  return(aoi)
}

#' Load Landsat data from folder
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
loadImage <-  function(path=getwd(), sat="auto", aoi){
  if(sat=="auto"){sat = getSat(path)} #DRY!
  if(sat=="L8"){bands <- 2:7}
  if(sat=="L7"){bands <- c(1:5,7)}
  files <- list.files(path = path, pattern = "^L[EC]\\d+\\w+\\d+_B\\d{1}.TIF$", 
                   full.names = T)
  files <- substr(files,1,nchar(files)-5)
  stack1 <- list()
  for(i in 1:6){
    stack1[i] <- raster(paste0(files[i], bands[i], ".TIF"))
  }
  raw.image <- do.call(stack, stack1)
  raw.image <- aoiCrop(raw.image, aoi)                               
  raw.image <- saveLoadClean(imagestack = raw.image, 
                               stack.names = c("B", "G", "R", "NIR", "SWIR1", "SWIR2"), 
                               file = paste0(result.folder, "image_DN.tif"), 
                               overwrite=TRUE)
  return(raw.image) 
}  

#' Calculates Top of atmosphere reflectance
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' LPSO. (2004). Landsat 7 science data users handbook, Landsat Project Science Office, NASA Goddard Space Flight Center, Greenbelt, Md., (http://ltpwww.gsfc.nasa.gov/IAS/handbook/handbook_toc.html) (Feb. 5, 2007)
#' @export
calcTOAr <- function(path=getwd(), image.DN, sat="auto", 
                      ESPA=FALSE, aoi, incidence.rel){
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
  }
  ### Ro TOA L7
  if(sat=="L7"){
    L5_ESUN <- c(1957, 1826, 1554, 1036, 215.0, 80.67) #Chandler & Markham 2003
    L7_ESUN <- c(1970, 1842, 1547, 1044, 225.7, 82.06)
    ESUN <- L7_ESUN 
    #L5_Gain <- c(0.762824, 1.442510, 1.039880, 0.872588, 0.119882, 0.055158, 0.065294) # Chandler 2003
    Gain <- c(1.181, 1.210, 0.943, 0.969, 0.191, 0.066)
    #L5_Bias <- c(-1.52, -2.84, -1.17, -1.51, -0.37, 1.2387, -0.15)
    Bias <- c(-7.38071, -7.60984, -5.94252, -6.06929, -1.19122, -0.41650)
    if(missing(image.DN)){image.DN <- loadImage(path = path)}
    DOY <- as.integer(substr(list.files(path = path, 
                                        pattern = "^L[EC]\\d+\\w+\\d+_B\\d{1}.TIF$")[1],14,16))
    d <- sqrt(1/(1+0.033*cos(DOY * 2 * pi/365)))
    Ro.TOAr <- list()
    for(i in 1:6){
      Ro.TOAr[i] <- pi * Gain[i] * image.DN[[i]] + Bias[i] * d^2 / ESUN[i] * cos(incidence.rel)
    }
    image_TOA <- do.call(stack, Ro.TOAr)
  }
  #### 
  image_TOA <- saveLoadClean(imagestack = image_TOA, 
                               stack.names = c("B", "G", "R", "NIR", "SWIR1", "SWIR2"), 
                               file = paste0(result.folder, "image_TOAr.tif"), 
                               overwrite=TRUE)
  return(image_TOA)
}  

#' Calculates surface reflectance for L7
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
# incidence hor from TML?? 
calcSR <- function(path=getwd(), image.TOAr, sat="auto", ESPA=FALSE, format="tif", 
                     aoi, incidence.hor, 
                    WeatherStation, surface.model){
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
    }
  if(sat=="L7"){
    if(missing(image.TOAr)){image.TOAr <- calcTOAr(path = path)}
    P <- 101.3*((293-0.0065 * surface.model$DEM)/293)^5.26
    ea <- (WeatherStation$RH/100)*0.6108*exp((17.27*WeatherStation$Ta)/(WeatherStation$Ta+237.3))
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
      tau_in[i] <- Cnb[1,i] * exp((Cnb[2,i]*P/Kt*cos(incidence.hor))-
                                    ((Cnb[3,i]*W+Cnb[4,i])/cos(incidence.hor)))+Cnb[5,i]
    }
    eta = 0 # eta it's the satellite nadir angle
    for(i in 1:6){
      tau_out[i] <- Cnb[1,i] * exp((Cnb[2,i]*P/Kt*cos(eta))-
                                    ((Cnb[3,i]*W+Cnb[4,i])/cos(eta)))+Cnb[5,i]
    }
    path_refl <- list()
    for(i in 1:6){
      path_refl[i] <- Cnb[6,1] * (1 - tau_in[[i]])
    }
    stack_SR <- list()
    for(i in 1:6){
      stack_SR[i] <- (image.TOAr[[i]] - path_refl[[i]]) / (tau_in[[i]] * tau_out[[i]])
    }
    image_SR <- do.call(stack, stack_SR)
  }
  image_SR <- saveLoadClean(imagestack = image_SR, 
                              stack.names = c("B", "G", "R", "NIR", "SWIR1", "SWIR2"), 
                              file = paste0(result.folder, "image_SR.tif"), 
                              overwrite=TRUE)
  return(image_SR)
}  


#' Check needed SRTM grids from image extent
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @export
# Get links or optionally open web pages... 
# Check if the files are present on path o in a specific SRTM local repo
checkSRTMgrids <-function(raw.image, path = getwd(), format="tif"){
  polyaoi <- SpatialPolygons(
    list(Polygons(list(Polygon(coords = matrix(
      c(xmin(raw.image), xmax(raw.image), xmax(raw.image),
        xmin(raw.image),xmin(raw.image), ymax(raw.image),
        ymax(raw.image), ymin(raw.image), ymin(raw.image),
        ymax(raw.image)), ncol=2, nrow= 5))), ID=1)))
  polyaoi@proj4string <- raw.image@crs
  limits <- proj4::project(xy = matrix(polyaoi@bbox, ncol=2, nrow=2), proj = polyaoi@proj4string, 
                    inverse = TRUE)
  # I have to improve this. It should work ONLY for west and south coordinates.. maybe
  lat_needed <- seq(trunc(limits[3])-1, trunc(limits[4])-1, by=1)
  long_needed <- seq(trunc(limits[1])-1, trunc(limits[2])-1, by = 1)
  grids <- expand.grid(lat_needed, long_needed)
  result <- list()
  link <- "http://earthexplorer.usgs.gov/download/options/8360/SRTM1"
  for(i in 1:nrow(grids)){
    result[[i]] <- paste(link, ifelse(grids[i,1]>0,"N", "S"), abs(grids[i,1]),
                         ifelse(grids[i,2]>0,"E", "W"), "0", abs(grids[i,2]),"V3/", sep="")
  }
  print(paste("You need", nrow(grids), "1deg x 1deg SRTM grids"))
  print("You can get them here:")
  return(unlist(result))
}

#' Create a mosaic with SRTM grid from image extent
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @export
# Should use checkSRTMgrids to get the files list and not use all from the folder...!
# Also look for files on path and local repo
prepareSRTMdata <- function(path=getwd(), format="tif", extent=image.DN){
  files <- list.files(path= path,  pattern=paste("^[sn]\\d{2}_[we]\\d{3}_1arc_v3.", 
                                                 format, "$", sep=""), full.names = T) 
  stack1 <- list()
  for(i in 1:length(files)){
    stack1[[i]] <- raster(files[i])}
  stack1$fun <- mean
  if(length(files)>1){SRTMmosaic <- do.call(mosaic, stack1)}
  if(length(files)==1){SRTMmosaic <- stack1[[1]]}
  destino  <-  projectExtent(extent, extent@crs)
  mosaicp <- projectRaster(SRTMmosaic, destino)
  mosaicp <- saveLoadClean(imagestack = mosaicp, stack.names = "DEM", 
                             file = paste0(result.folder, "DEM.tif"), overwrite=TRUE)
  return(mosaicp)
}

#' Calculates surface model used in METRIC
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
METRICtopo <- function(DEM){
  aspect <- terrain(DEM, opt="aspect") 
  slope <- terrain(DEM, opt="slope") 
  aspect_metric <- aspect-pi  #METRIC expects aspect - 1 pi
  surface.model <- stack(DEM, slope, aspect_metric)
  surface.model <- saveLoadClean(imagestack = surface.model, 
                                   stack.names = c("DEM", "Slope", "Aspect"), 
                                   file = paste0(result.folder, "surface.model.tif"), 
                                   overwrite=TRUE)
  return(surface.model)
}

#' Calculates solar angles
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
### Change to look in metadata for keyword instead of using line #
  solarAngles <- function(path=getwd(), surface.model, MTL){
  if(missing(MTL)){Landsat.MTL <- list.files(path = path, pattern = "MTL.txt", full.names = T)}
  MTL <- readLines(Landsat.MTL, warn=FALSE)
  Elev.line <- grep("SUN_ELEVATION",MTL,value=TRUE)
  sun.elevation <- as.numeric(regmatches(Elev.line, 
                                         regexec(text=Elev.line ,
                                                 pattern="([0-9]{1,5})([.]+)([0-9]+)"))[[1]][1])
  Azim.line <- grep("SUN_AZIMUTH",MTL,value=TRUE)
  sun.azimuth <- as.numeric(regmatches(Azim.line, 
                                       regexec(text=Azim.line ,
                                               pattern="([0-9]{1,5})([.]+)([0-9]+)"))[[1]][1])
  # latitude
  latitude <- surface.model[[1]]
  xy <- SpatialPoints(xyFromCell(latitude, cellFromRowCol(latitude, 1:nrow(latitude), 1)))
  xy@proj4string <- latitude@crs
  lat <- coordinates( spTransform(xy, CRS("+proj=longlat +datum=WGS84")))[,2] 
  values(latitude) <- rep(lat*pi/180,each=ncol(latitude))
  # declination
  DOY <- as.integer(substr(list.files(path = path, 
                                      pattern = "^L[EC]\\d+\\w+\\d+_B\\d{1}.TIF$")[1],
                           14,16))
  declination <- surface.model[[1]]
  values(declination) <- 23.45*pi/180*sin(2*pi*((284+DOY)/36.25))
  # hour angle
  hour.angle <- asin(-1*(cos(sun.elevation)*sin(sun.azimuth)/cos(declination)))
  ## solar incidence angle, for horizontal surface
  incidence.hor <- acos(sin(declination) * sin(latitude) + cos(declination)
                        *cos(latitude)*cos(hour.angle))
  slope <- surface.model$Slope
  aspect <- surface.model$Aspect
  ##solar incidence angle, for sloping surface
  incidence.rel <- acos(sin(declination)*sin(latitude)*cos(slope) 
                        - sin(declination)*cos(latitude)*sin(slope)*cos(aspect)
                        + cos(declination)*cos(latitude)*cos(slope)*cos(hour.angle)
                        + cos(declination)*sin(latitude)*sin(slope)*cos(aspect)*cos(hour.angle)
                        + cos(declination)*sin(aspect)*sin(slope)*sin(hour.angle))
  ## End
  solarAngles <- stack(latitude, declination, hour.angle, incidence.hor, incidence.rel)
  solarAngles <- saveLoadClean(imagestack = solarAngles, 
                                  stack.names = c("latitude", "declination", 
                                                  "hour.angle", "incidence.hor", "incidence.rel"), 
                                  file = paste0(result.folder, "solarAngles.tif"), overwrite=TRUE)
  return(solarAngles)
}

#' Calculates Incoming Solar Radiation
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
incSWradiation <- function(surface.model, solar.angles, WeatherStation){
  ea <- (WeatherStation$RH/100)*0.6108*exp((17.27*WeatherStation$Ta)/(WeatherStation$Ta+237.3))
  tau.sw <- SWtrasmisivity(Kt = 1, ea, surface.model$DEM, solar.angles$incidence.hor)
  DOY <- WeatherStation$DOY
  d <- sqrt(1/(1+0.033*cos(DOY * 2 * pi/365)))
  Rs.inc <- 1367 * cos(solar.angles$incidence.rel) * tau.sw / d^2
  Rs.inc <- saveLoadClean(imagestack = Rs.inc, 
                            file = paste0(result.folder, "Rs.inc.tif"), overwrite=TRUE)
  return(Rs.inc)
}

#' Calculates Broadband Albedo from Landsat data
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
albedo <- function(path=getwd(), image.SR, aoi, coeff="Tasumi", sat="auto", ESPA=FALSE){
  if(sat=="auto"){sat = getSat(path)}
  if(sat=="L8"){bands <- 2:7}
  if(sat=="L7"){bands <- c(1:5,7)}
  if(coeff=="Tasumi"){wb <- c(0.254, 0.149, 0.147, 0.311, 0.103, 0.036) * 10000} 
  # Tasumi 2008
  if(coeff=="Olmedo") {wb <- c(0.246, 0.146, 0.191, 0.304, 0.105, 0.008) * 10000 }
  # Calculated using SMARTS for Kimberly2-noc13 and Direct Normal Irradiance
  if(coeff=="Liang") {wb <- c(0.356, 0, 0.130, 0.373, 0.085, 0.072) * 10000} 
  # Liang 2001
  if(ESPA==TRUE & sat=="L8"){
  files <- list.files(path = path, pattern = "_sr_band+[2-7].tif$", full.names = T)
  albedo <- calc(raster(files[1]), fun=function(x){x *wb[1]})
  for(i in 2:6){
    albedo <- albedo + calc(raster(files[i]), fun=function(x){x *wb[i]})
    removeTmpFiles(h=0.0008) # delete last one... maybe 3 seconds
  }
  albedo <-  albedo/1e+08}
  if(sat=="L7"){
    albedo <- calc(image.SR[[1]], fun=function(x){x *wb[1]/10})
    for(i in 2:6){
      albedo <- albedo + calc(image.SR[[i]], fun=function(x){x *wb[i]/10})
      removeTmpFiles(h=0.0008) # delete last one... maybe 3 seconds
    }
  albedo <-  albedo/1e+06}
  albedo <- aoiCrop(albedo, aoi) 
  if(coeff=="Liang"){
    albedo <- albedo - 0.0018
  }
  albedo <- saveLoadClean(imagestack = albedo, 
                            file = paste0(result.folder, "albedo.tif"), overwrite=TRUE)
  return(albedo)
}

#' Estimate LAI from Landsat Data
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
## Cite Pocas work for LAI from METRIC2010
LAI <- function(method="metric2010", path=getwd(), aoi, L=0.1, ESPA=F, image, sat="auto"){
  if(sat=="auto"){sat = getSat(path)}
  if(sat=="L8" & ESPA==T){
    if(method=="metric" | method=="metric2010" | method=="vineyard" | method=="MCB"){
      removeTmpFiles(h=0)
      toa.red <- raster(list.files(path = path, pattern = "_toa_band4.tif", full.names = T))
      toa.nir <- raster(list.files(path = path, pattern = "_toa_band5.tif", full.names = T))
      toa.4.5 <- stack(toa.red, toa.nir)
      toa.4.5 <- aoiCrop(toa.4.5, aoi)
      toa.4.5 <- toa.4.5 * 0.0001
      }
    if(method=="turner"){
      removeTmpFiles(h=0)
      sr.red <- raster(list.files(path = path, pattern = "_sr_band4.tif", full.names = T))
      sr.nir <- raster(list.files(path = path, pattern = "_sr_band5.tif", full.names = T))
      sr.4.5 <- stack(sr.red, sr.nir)
      sr.4.5 <- aoiCrop(sr.4.5, aoi)}
  }
  if(sat=="L7"){
    if(method=="metric" | method=="metric2010" | method=="vineyard" | method=="MCB"){
      toa.4.5 <- stack(image[[3]], image[[4]])}
    if(method=="turner"){
      sr.4.5 <- stack(image[[3]], image[[4]])
    }
  }  
  if(method=="turner"){
    sr.4.5 <- sr.4.5 * 0.0001
  }
  if(method=="metric2010"){
    SAVI_ID <- (1 + L)*(toa.4.5[[2]] - toa.4.5[[1]])/(L + toa.4.5[[1]] + toa.4.5[[2]])
    LAI <- 11 * SAVI_ID ^3 # for SAVI <= 0.817
    LAI[SAVI_ID > 0.817] <- 6
  }
  if(method=="Bastiaanssen"){
    SAVI_ID <- (1 + L)*(toa.4.5[[2]] - toa.4.5[[1]])/(L + toa.4.5[[1]] + toa.4.5[[2]])
    LAI <- log((0.69-SAVI_ID)/0.59)/0.91 *-1
    LAI[SAVI_ID > 0.817] <- 6
  }
  if(method=="vineyard"){
    NDVI <- (toa.4.5[[2]] - toa.4.5[[1]])/(toa.4.5[[1]] + toa.4.5[[2]])
    LAI <- 4.9 * NDVI -0.46 # Johnson 2003
  }
  ## method carrasco
  if(method=="MCB"){
    NDVI <- (toa.4.5[[2]] - toa.4.5[[1]])/(toa.4.5[[1]] + toa.4.5[[2]])
    LAI <- 1.2 - 3.08*exp(-2013.35*NDVI^6.41) 
  }
  if(method=="turner"){
    NDVI <- (toa.4.5[[2]] - toa.4.5[[1]])/(toa.4.5[[1]] + toa.4.5[[2]])
    LAI <- 0.5724+0.0989*NDVI-0.0114*NDVI^2+0.0004*NDVI^3 # Turner 1999
  }
  LAI[LAI<0]  <- 0
  LAI <- saveLoadClean(imagestack = LAI, stack.names = "LAI", 
                         file = paste0(result.folder, "LAI.tif"), overwrite=TRUE)
  return(LAI)
}

#' Estimates Land Surface Temperature from Landsat Data
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
## Add Sobrino and Qin improvements to LST in ETM+
## Add Rsky estimation from WeatherStation
surfaceTemperature <- function(path=getwd(), sat="auto", LAI, aoi){
  if(sat=="auto"){sat = getSat(path)}
  if(sat=="L8"){
    bright.temp.b10 <- raster(list.files(path = path, 
                                                     pattern = "_toa_band10.tif"))
    bright.temp.b10 <- aoiCrop(bright.temp.b10, aoi) 
    Ts <- bright.temp.b10*0.1
  }
  if(sat=="L7"){
    epsilon_NB <- raster(LAI)
    epsilon_NB <- 0.97 + 0.0033 * LAI  
    epsilon_NB[LAI > 3] <- 0.98
    L_t_6 <-  0.067 * raster(list.files(path = path, pattern = "^L[EC]\\d+\\w+\\d+_B6_VCID_1.TIF$", full.names = T)) + 3.1628
    L7_K1 <- 666.09 
    L7_K2 <- 1282.71 
    Rp <- 0           #Allen estimo en Idaho que el valor medio era 0.91
    tau_NB <- 1       #Allen estimo en Idaho que el valor medio era 0.866
    R_sky <- 1        #Allen estimo en Idaho que el valor medio era 1.32
    Rc <- ((L_t_6 - Rp) / tau_NB) - (1-epsilon_NB)/R_sky
    Ts <- L7_K2 / log((epsilon_NB*L7_K1/Rc)+1)}
  Ts <- saveLoadClean(imagestack = Ts, 
                        file = paste0(result.folder, "Ts.tif"), overwrite=TRUE)
  return(Ts)
}

#' Calculates Long wave outgoing radiation
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
outLWradiation <- function(LAI, Ts){
  surf.emissivity <- 0.95 + 0.01 * LAI
  surf.emissivity[LAI>3] <- 0.98
  Rl.out <- surf.emissivity * 5.67e-8 * Ts^4
  Rl.out <- saveLoadClean(imagestack = Rl.out, 
                            file = paste0(result.folder, "Rs.out.tif"), overwrite=TRUE)
  return(Rl.out)
}

#' Calculates long wave incoming radiation
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
# Add a function to get "cold" pixel temperature so in can be used in the next function
incLWradiation <- function(WeatherStation, DEM, solarAngles){
  ea <- (WeatherStation$RH/100)*0.6108*exp((17.27*WeatherStation$Ta)/
                                             (WeatherStation$Ta+237.3))
  tau.sw <- SWtrasmisivity(Kt = 1, ea, DEM, solar.angles$incidence.hor)
  Ta <-  WeatherStation$Ta+273.15 - (DEM - WeatherStation$elev) * 6.49 / 1000 
  ## Temperature in Kelvin
  #Mountain lapse effects from International Civil Aviation Organization
  ef.atm.emissivity  <- 0.85*(-1*log(tau.sw))^0.09
  Rl.in <- ef.atm.emissivity * 5.67e-8 * Ta^4
  Rl.in <- saveLoadClean(imagestack = Rl.in, 
                           file = paste0(result.folder, "Rl.in.tif"), overwrite=TRUE)
  return(Rl.in)
}

#' Estimates Soil Heat Flux
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
soilHeatFlux <- function(path=getwd(), image, Ts, albedo, Rn, aoi, sat="auto", ESPA=F){
  if(sat=="auto"){sat = getSat(getwd())}
  if(sat=="L8" & ESPA==T){
      removeTmpFiles(h=0)
      sr.red <- raster(list.files(path = path, pattern = "_sr_band4.tif", full.names = T))
      sr.nir <- raster(list.files(path = path, pattern = "_sr_band5.tif", full.names = T))
      sr.4.5 <- stack(sr.red, sr.nir)
      sr.4.5 <- aoiCrop(sr.4.5, aoi)
  }
  if(sat=="L7"){
      sr.4.5 <- stack(image[[3]], image[[4]])
  }
  NDVI <- (sr.4.5[[2]] - sr.4.5[[1]])/(sr.4.5[[1]] + sr.4.5[[2]])
  G <- ((Ts - 273.15)*(0.0038+0.0074*albedo)*(1-0.98*NDVI^4))*Rn
  G <- saveLoadClean(imagestack = G, file = paste0(result.folder, "G.tif"), overwrite=TRUE)
  return(G)
}

#' Calculates Momentum Roughness Length
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
## Create a function to estimate a and b coefficients or the function between Z.om and NDVI
## using some points and tabulated z.om for their covers.
## Perrier by Santos 2012 and Pocas 2014.
momentumRoughnessLength <- function(method="short.crops", LAI, NDVI, 
                                      albedo, a, b, fLAI.Perrier, h.Perrier, 
                                      mountainous=FALSE, surface.model){
  if(method=="short.crops"){
    Z.om <- (0.018*LAI)
  }
  if(method=="custom"){
    Z.om <- exp((a*NDVI/albedo)+b)
  }
  if(method=="Perrier"){
    if(fLAI.Perrier >=0.5){ a <- (2*(1-fLAI.Perrier))^-1 }
    if(fLAI.Perrier <0.5){ a <- 2*fLAI.Perrier }
    Z.om <- ((1-exp(-a*LAI/2))*exp(-a*LAI/2))^h.Perrier
  }
  if(mountainous==TRUE){
    Z.om <- Z.om * (1 + (180/pi*surface.model$Slope - 5)/20)
  }
  Z.om <- saveLoadClean(imagestack = Z.om, 
                          file = paste0(result.folder, "Z.om.tif"), overwrite=TRUE)
  return(Z.om)
}


