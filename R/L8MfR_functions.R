### Add a function to estimate all parameters for some points (and only for those points!)
### Add some meta-functions to insolate processes
### Add METRIC function
### Maybe we need a new class to store all data? 
### Maybe add knitr and generate pdf report on METRIC function
### Export anchor as kml or view in Google
### select anchor points with buffer for not too close pixels...
### Select anchor for multiple criteria with table from Marcos
### Check Rsky on METRIC 2010
### Add three source temperature model..!
### A function to get a template: file.copy(system.file('test.R','MyPackage'), '.')

#################################################################################3
# Common fields for docs
# @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
# @references 
# R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
# Description
# @export
# @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
####################

#' Create aoi polygon from topleft and bottomright coordinates
#' @param topleft a vector with topleft x,y coordinates 
#' @param bottomright a vector with bottomright x,y coordinates
#' @return object of class SpatialPolygons
#' @author Guillermo F Olmedo
#' @examples 
#' tl <- c(493300, -3592700)
#' br <- c(557200, -3700000) 
#' aoi <- create.aoi(topleft = tl, bottomright=br)
#' plot(aoi)
#' @export
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @seealso \code{\link{brocolors}}
#' @keywords hplot
#' ...
#' @importFrom grDevices rgb2hsv
#' @importFrom graphics par plot rect text
# Maybe i can provide some points and use CHull like in QGIS-Geostat
create.aoi <- function(topleft = c(x, y), bottomright= c(x, y)){
  aoi <- SpatialPolygons(
    list(Polygons(list(Polygon(coords = matrix(
      c(topleft[1],bottomright[1], bottomright[1],topleft[1],topleft[1],
        topleft[2], topleft[2], bottomright[2], 
        bottomright[2],topleft[2]), ncol=2, nrow= 5))), ID=1)))
  return(aoi)
}

#' Load Landsat 8 data from a folder 
#' @export
load.image.DN <-  function(path=getwd(), sat="auto", result.folder=NULL){
  if(sat=="auto"){sat = get.sat(path)} #DRY!
  if(sat=="L8"){bands <- 2:7}
  if(sat=="L7"){bands <- c(1:5,7)}
  files <- substr(list.files(path = path, pattern = "^L[EC]\\d+\\w+\\d+_B\\d{1}.TIF$"),0,23)
  stack1 <- list()
  for(i in 1:6){
    stack1[i] <- raster(paste0(path, files[i], bands[i], ".TIF"))
  }
  raw.image <- do.call(stack, stack1)
  if(!is.null(aoi)){
    raw.image <- crop(raw.image,aoi) # Without aoi this should fail on most computers.
  }                                
  raw.image <- save.load.clean(imagestack = raw.image, 
                               stack.names = c("B", "G", "R", "NIR", "SWIR1", "SWIR2"), 
                               file = paste0(result.folder, "image_DN.tif"), 
                               overwrite=TRUE)
  return(raw.image) 
}  


# References L7: LPSO. (2004). Landsat 7 science data users handbook, Landsat Project Science Office, NASA Goddard Space Flight Center, Greenbelt, Md., (http://ltpwww.gsfc.nasa.gov/IAS/handbook/handbook_toc.html) (Feb. 5, 2007)
calc.TOAr <- function(path=getwd(), image.DN, sat="auto", ESPA=FALSE, aoi, result.folder=NULL, incidence.rel){
  if(sat=="auto"){sat = get.sat(path)}
  if(sat=="L8"){bands <- 2:7}
  if(sat=="L7"){bands <- c(1:5,7)}
  if(ESPA==TRUE & sat=="L8"){
    files <- list.files(path = path, pattern = "_toa_band+[2-7].tif$")
    stack1 <- list()
    for(i in 1:6){
      stack1[i] <- raster(paste0(path, files[i]))
    }
    image_TOA <- do.call(stack, stack1)
    if(!is.null(aoi)){
      image_TOA <- crop(image_TOA,aoi) 
    }}
  ### Ro TOA L7
  if(sat=="L7"){
    L5_ESUN <- c(1957, 1826, 1554, 1036, 215.0, 80.67) #Chandler & Markham 2003
    L7_ESUN <- c(1970, 1842, 1547, 1044, 225.7, 82.06)
    ESUN <- L7_ESUN 
    #L5_Gain <- c(0.762824, 1.442510, 1.039880, 0.872588, 0.119882, 0.055158, 0.065294) # Chandler 2003
    Gain <- c(1.181, 1.210, 0.943, 0.969, 0.191, 0.066)
    #L5_Bias <- c(-1.52, -2.84, -1.17, -1.51, -0.37, 1.2387, -0.15)
    Bias <- c(-7.38071, -7.60984, -5.94252, -6.06929, -1.19122, -0.41650)
    if(missing(image.DN)){image.DN <- load.image.DN(path = path)}
    DOY <- as.integer(substr(list.files(path = path, pattern = "^L[EC]\\d+\\w+\\d+_B\\d{1}.TIF$")[1],14,16))
    d <- sqrt(1/(1+0.033*cos(DOY * 2 * pi/365)))
    Ro.TOAr <- list()
    for(i in 1:6){
      Ro.TOAr[i] <- pi * Gain[i] * image.DN[[i]] + Bias[i] * d^2 / ESUN[i] * cos(incidence.rel)
    }
    image_TOA <- do.call(stack, Ro.TOAr)
  }
  #### 
  image_TOA <- save.load.clean(imagestack = image_TOA, 
                               stack.names = c("B", "G", "R", "NIR", "SWIR1", "SWIR2"), 
                               file = paste0(result.folder, "image_TOAr.tif"), 
                               overwrite=TRUE)
  return(image_TOA)
}  

# incidence hor from TML?? 
calc.SR <- function(path=getwd(), image.TOAr, sat="auto", ESPA=FALSE, format="tif", 
                     aoi, result.folder=NULL, incidence.hor, WeatherStation, surface.model){
  if(sat=="auto"){sat = get.sat(path)}
  if(sat=="L8"){bands <- 2:7}
  if(sat=="L7"){bands <- c(1:5,7)}
  if(ESPA==TRUE & sat=="L8"){
    files <- list.files(path = path, pattern = "_sr_band+[2-7].tif$")
    stack1 <- list()
    for(i in 1:6){
      stack1[i] <- raster(paste0(path, files[i]))
    }
    image_SR <- do.call(stack, stack1)
    if(!is.null(aoi)){
      image_SR <- crop(image_SR,aoi) 
    }}
  if(sat=="L7"){
    if(missing(image.TOAr)){image.TOAr <- calc.TOAr(path = path)}
    P <- 101.3*((293-0.0065 * surface.model$DEM)/293)^5.26
    ea.sat <- 0.6108*exp((17.27*WeatherStation$Ta)/(WeatherStation$Ta+237.3))
    ea <- (WeatherStation$RH/100)*ea.sat
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
  image_SR <- save.load.clean(imagestack = image_SR, 
                              stack.names = c("B", "G", "R", "NIR", "SWIR1", "SWIR2"), 
                              file = paste0(result.folder, "image_SR.tif"), 
                              overwrite=TRUE)
  return(image_SR)
}  



# Get links or optionally open web pages... 
checkSRTMgrids <-function(raw.image, path = getwd(), format="tif"){
  polyaoi <- SpatialPolygons(
    list(Polygons(list(Polygon(coords = matrix(
      c(xmin(raw.image), xmax(raw.image), xmax(raw.image),
        xmin(raw.image),xmin(raw.image), ymax(raw.image),
        ymax(raw.image), ymin(raw.image), ymin(raw.image),
        ymax(raw.image)), ncol=2, nrow= 5))), ID=1)))
  polyaoi@proj4string <- raw.image@crs
  limits <- project(xy = matrix(polyaoi@bbox, ncol=2, nrow=2), proj = polyaoi@proj4string, 
                    inverse = TRUE)
  # I have to improve this. It should work ONLY for west and south coordinates.. maybe
  lat_needed <- seq(floor(limits[3])+1, floor(limits[4])+1, by=1)
  long_needed <- seq(floor(limits[1]), floor(limits[2])+1, by = 1)
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

# Should use checkSRTMgrids to get the files list and not use all from the folder...!
prepareSRTMdata <- function(path=getwd(), format="tif", extent=raw.image, result.folder=NULL){
  files <- list.files(path= path,  pattern=paste("^[sn]\\d{2}_[we]\\d{3}_1arc_v3.", 
                                                 format, "$", sep="")) 
  stack1 <- list()
  for(i in 1:length(files)){
    stack1[[i]] <- raster(paste(path, files[i], sep=""))}
  stack1$fun <- mean
  SRTMmosaic <- do.call(mosaic, stack1)
  destino  <-  projectExtent(raw.image, raw.image@crs)
  mosaicp <- projectRaster(SRTMmosaic, destino)
  mosaicp <- save.load.clean(imagestack = mosaicp, stack.names = "DEM", 
                             file = paste0(result.folder, "DEM.tif"), overwrite=TRUE)
  return(mosaicp)
}

METRIC.topo <- function(DEM, result.folder=NULL){
  DEM <- raster(DEM)
  aspect <- terrain(DEM, opt="aspect") 
  slope <- terrain(DEM, opt="slope") 
  aspect_metric <- aspect-pi  #METRIC expects aspect - 1 pi
  surface.model <- stack(DEM, slope, aspect_metric)
  surface.model <- save.load.clean(imagestack = surface.model, 
                                   stack.names = c("DEM", "Slope", "Aspect"), 
                                   file = paste0(result.folder, "surface.model.tif"), 
                                   overwrite=TRUE)
  return(surface.model)
}


### Change to look in metadata for keyword instead of using line #
solar.angles <- function(L8MTL, raw.image, slope, aspect, result.folder=NULL){
  test <- scan(L8MTL, character(0), sep = "\n") 
  ## Here uses stringr,.. only 1 line.. doesn't justify adding a dependency!
  sun.azimuth <- as.numeric(str_extract(test[68], 
                                        pattern = "([0-9]{1,5})([.]+)([0-9]+)"))*pi/180
  sun.elevation <- as.numeric(str_extract(test[69], 
                                          pattern = "([0-9]{1,5})([.]+)([0-9]+)"))*pi/180
  # latitude
  latitude <- raw.image[[1]]
  xy <- SpatialPoints(xyFromCell(latitude, cellFromRowCol(latitude, 1:nrow(latitude), 1)))
  xy@proj4string <- latitude@crs
  lat <- coordinates( spTransform(xy, CRS("+proj=longlat +datum=WGS84")))[,2] 
  values(latitude) <- rep(lat*pi/180,each=ncol(latitude))
  # declination
  DOY <- strptime(str_extract(test[21], 
                              pattern = "([0-9]{4})([-]+)([0-9]{2})([-]+)([0-9]{2})"), 
                  "%Y-%m-%d")$yday+1
  declination <- raw.image[[1]]
  values(declination) <- 23.45*pi/180*sin(2*pi*((284+DOY)/36.25))
  # hour angle
  hour.angle <- asin(-1*(cos(sun.elevation)*sin(sun.azimuth)/cos(declination)))
  ## solar incidence angle, for horizontal surface
  incidence.hor <- acos(sin(declination) * sin(latitude) + cos(declination)
                        *cos(latitude)*cos(hour.angle))
  ##solar incidence angle, for sloping surface
  incidence.rel <- acos(sin(declination)*sin(latitude)*cos(slope) 
                        - sin(declination)*cos(latitude)*sin(slope)*cos(aspect)
                        + cos(declination)*cos(latitude)*cos(slope)*cos(hour.angle)
                        + cos(declination)*sin(latitude)*sin(slope)*cos(aspect)*cos(hour.angle)
                        + cos(declination)*sin(aspect)*sin(slope)*sin(hour.angle))
  ## End
  solar.angles <- stack(latitude, declination, hour.angle, incidence.hor, incidence.rel)
  solar.angles <- save.load.clean(imagestack = solar.angles, 
                                  stack.names = c("latitude", "declination", 
                                                  "hour.angle", "incidence.hor", "incidence.rel"), 
                                  file = paste0(result.folder, "solar.angles.tif"), overwrite=TRUE)
  return(solar.angles)
}

sw.trasmisivity <- function(Kt = 1, ea, dem, incidence.hor, result.folder=NULL){
  P <- 101.3*((293-0.0065 * dem)/293)^5.26
  W <- 0.14 * ea * P + 2.1
  sw.t <- 0.35 + 0.627 * exp((-0.00149 * P / Kt * cos(incidence.hor))-0.075*(W / cos(incidence.hor))^0.4)
  sw.t <- save.load.clean(imagestack = sw.t, 
                          file = paste0(result.folder, "sw.t.tif"), overwrite=TRUE)
  return(sw.t)
}

incoming.solar.radiation <- function(incidence.rel, tau.sw, DOY, result.folder=NULL){
  d <- sqrt(1/(1+0.033*cos(DOY * 2 * pi/365)))
  Rs.inc <- 1367 * cos(incidence.rel) * tau.sw / d^2
  Rs.inc <- save.load.clean(imagestack = Rs.inc, 
                            file = paste0(result.folder, "Rs.inc.tif"), overwrite=TRUE)
  return(Rs.inc)
}

albedo <- function(path=getwd(), aoi, coeff="Tasumi", result.folder=NULL){
  if(coeff=="Tasumi"){wb <- c(0.254, 0.149, 0.147, 0.311, 0.103, 0.036) * 10000} 
  # Tasumi 2008
  if(coeff=="Olmedo") {wb <- c(0.246, 0.146, 0.191, 0.304, 0.105, 0.008) * 10000 }
  # Calculated using SMARTS for Kimberly2-noc13 and Direct Normal Irradiance
  if(coeff=="Liang") {wb <- c(0.356, 0, 0.130, 0.373, 0.085, 0.072) * 10000} 
  # Liang 2001
  files <- list.files(path = path, pattern = "_sr_band+[2-7].tif$")
  albedo <- calc(raster(paste(path, files[1], sep="")), fun=function(x){x *wb[1]})
  for(i in 2:6){
    albedo <- albedo + calc(raster(paste(path, files[i], sep="")), 
                            fun=function(x){x *wb[i]})
    removeTmpFiles(h=0.0008) # delete last one... maybe 3 seconds
  }
  albedo <-  albedo/1e+08
  if(!missing(aoi)){
    albedo <- crop(albedo,aoi) # Without aoi this should fail on most computers.
  }                                
  if(coeff=="Liang"){
    albedo <- albedo - 0.0018
  }
  albedo <- save.load.clean(imagestack = albedo, 
                            file = paste0(result.folder, "albedo.tif"), overwrite=TRUE)
  return(albedo)
}

## Cite Pocas work for LAI from METRIC2010
LAI.from.L8 <- function(method="metric2010", path=getwd(), aoi, L=0.1, result.folder=NULL){
  if(method=="metric" | method=="metric2010" | method=="vineyard" | method=="carrasco-benavides"){
    toa.red <- raster(paste(path, list.files(path = path, pattern = "_toa_band4.tif"), sep=""))
    toa.nir <- raster(paste(path, list.files(path = path, pattern = "_toa_band5.tif"), sep=""))
    toa.4.5 <- stack(toa.red, toa.nir)
    if(!missing(aoi)){
      toa.4.5 <- crop(toa.4.5,aoi) # Without aoi this should fail on most computers.
    }    
    toa.4.5 <- toa.4.5 * 0.0001
  }
  if(method=="turner"){
    toa.red <- raster(paste(path, list.files(path = path, pattern = "_sr_band4.tif"), sep=""))
    toa.nir <- raster(paste(path, list.files(path = path, pattern = "_sr_band5.tif"), sep=""))
    toa.4.5 <- stack(toa.red, toa.nir) # It says toa, but they are the sr images
    if(!missing(aoi)){
      toa.4.5 <- crop(toa.4.5,aoi) # Without aoi this should fail on most computers.
    }    
    toa.4.5 <- toa.4.5 * 0.0001
  }
  if(method=="metric2010"){
    SAVI_ID <- (1 + L)*(toa.4.5[[2]] - toa.4.5[[1]])/(L + toa.4.5[[1]] + toa.4.5[[2]])
    LAI <- 11 * SAVI_ID ^3 # for SAVI <= 0.817
    LAI[SAVI_ID > 0.817] <- 6
  }
  if(method=="metric"){
    SAVI_ID <- (1 + L)*(toa.4.5[[2]] - toa.4.5[[1]])/(L + toa.4.5[[1]] + toa.4.5[[2]])
    LAI <- 11 * SAVI_ID ^3 # for SAVI <= 0.817
    LAI[SAVI_ID > 0.817] <- 6
  }
  if(method=="vineyard"){
    NDVI <- (toa.4.5[[2]] - toa.4.5[[1]])/(toa.4.5[[1]] + toa.4.5[[2]])
    LAI <- 4.9 * NDVI -0.46 # Johnson 2003
  }
  ## method carrasco
  if(method=="carrasco-benavides"){
    NDVI <- (toa.4.5[[2]] - toa.4.5[[1]])/(toa.4.5[[1]] + toa.4.5[[2]])
    LAI <- 1.2 - 3.08*exp(-2013.35*NDVI^6.41) 
  }
  if(method=="turner"){
    NDVI <- (toa.4.5[[2]] - toa.4.5[[1]])/(toa.4.5[[1]] + toa.4.5[[2]])
    LAI <- 0.5724+0.0989*NDVI-0.0114*NDVI^2+0.0004*NDVI^3 # Turner 1999
  }
  removeTmpFiles(h=0)
  LAI <- save.load.clean(imagestack = LAI, stack.names = "LAI", 
                         file = paste0(result.folder, "LAI.tif"), overwrite=TRUE)
  return(LAI)
}

surface.temperature <- function(path=getwd(), aoi, result.folder=NULL){
  bright.temp.b10 <- raster(paste(path, list.files(path = path, 
                                                   pattern = "_toa_band10.tif"), sep=""))
  if(!missing(aoi)){
    bright.temp.b10 <- crop(bright.temp.b10,aoi) 
  }
  Ts <- bright.temp.b10*0.1
  Ts <- save.load.clean(imagestack = Ts, 
                        file = paste0(result.folder, "Ts.tif"), overwrite=TRUE)
  return(Ts)
}

outgoing.lw.radiation <- function(path=getwd(), LAI, aoi, result.folder=NULL){
  Ts <- surface.temperature(path=path, aoi=aoi)
  surf.emissivity <- 0.95 + 0.01 * LAI ## And when LAI 3 or more = 0.98
  Rl.out <- surf.emissivity * 5.67e-8 * Ts^4
  Rl.out <- save.load.clean(imagestack = Rl.out, 
                            file = paste0(result.folder, "Rs.out.tif"), overwrite=TRUE)
  return(Rl.out)
}

# Add a function to get "cold" pixel temperature so in can be used in the next function
incoming.lw.radiation <- function(air.temperature, DEM, sw.trasmisivity, aoi, result.folder=NULL){
  Ta <-  air.temperature - (DEM - 702) * 6.49 / 1000 ## Temperature in Kelvin
  if(!missing(aoi)){
    Ta <- crop(Ta,aoi) 
  }
  ef.atm.emissivity  <- 0.85*(-1*log(sw.trasmisivity))^0.09
  Rl.in <- ef.atm.emissivity * 5.67e-8 * Ta^4
  Rl.in <- save.load.clean(imagestack = Rl.in, 
                           file = paste0(result.folder, "Rl.in.tif"), overwrite=TRUE)
  return(Rl.in)
}

soil.heat.flux1 <- function(path=getwd(), albedo, Rn, aoi, result.folder=NULL){
  toa.red <- raster(paste(path, list.files(path = path, 
                                           pattern = "_sr_band4.tif"), sep="")[1])
  toa.nir <- raster(paste(path, list.files(path = path, 
                                           pattern = "_sr_band5.tif"), sep="")[1])
  toa.4.5 <- stack(toa.red*0.0001, toa.nir*0.0001) ## chuncks
  if(!missing(aoi)){
    toa.4.5 <- crop(toa.4.5,aoi) # Without aoi this should fail on most computers.
  }
  NDVI <- (toa.4.5[[2]] - toa.4.5[[1]])/(toa.4.5[[1]] + toa.4.5[[2]])
  Ts <- surface.temperature(path=path, aoi=aoi)
  G <- ((Ts - 273.15)*(0.0038+0.0074*albedo)*(1-0.98*NDVI^4))*Rn
  G <- save.load.clean(imagestack = G, file = paste0(result.folder, "G.tif"), overwrite=TRUE)
  return(G)
}

## Create a function to estimate a and b coefficients or the function between Z.om and NDVI
## using some points and tabulated z.om for their covers.
## Perrier by Santos 2012 and Pocas 2014.
momentum.roughness.length <- function(method="short.crops", path=getwd(), LAI, NDVI, 
                                      albedo, a, b, fLAI.Perrier, h.Perrier, 
                                      mountainous=FALSE, surf.model, result.folder=NULL){
  if(method=="short.crops"){
    Z.om <- (0.018*LAI)
  }
  if(method=="custom"){
    Z.om <- exp((a*NDVI/albedo)+b)
  }
  if(method=="Perrier"){
    if(fLAI.Perrier >=0.5){ a <- (2*(1-fLAI.Perrier))^-1 }
    if(fLAI.Perrier <0.5){ a <- 2*fLAI.Perrier }
    Z.om <- ((1-exp(-a*LAI/2))*exp(-a*LAI/2))^h
  }
  if(mountainous==TRUE){
    Z.om <- Z.om * (1 + (180/pi*surf.model$Slope - 5)/20)
  }
  Z.om <- save.load.clean(imagestack = Z.om, 
                          file = paste0(result.folder, "Z.om.tif"), overwrite=TRUE)
  return(Z.om)
}

#' Calculates ET using Penman Monteith hourly formula
#' @param WeatherStation a data frame with all the needed fields (see example)
#' @param hours time of the day in hours in 24hs format
#' @param DOY day of year
#' @param long.z longitude for local time
#' @return ET hourly in mm.h-1
#' @author Guillermo F Olmedo
#' @examples 
#' WeatherStation  <- data.frame(wind=4.72,
#'                               RH=59, 
#'                               Ta=24.3,
#'                               Gr.Rad=675, 
#'                               height=2.2, 
#'                               lat=-35.37, 
#'                               long=71.5946, 
#'                               elev=124)
#'   ETo.PM.hourly(WeatherStation, hours=10.5, DOY=363, long.z=71.635)
#' @export
#' @references 
#' Allen 2005 ASCE
ETo.PM.hourly <- function(WeatherStation, hours, DOY, long.z=WeatherStation$long){
  TaK <- WeatherStation$Ta + 273.16
  Rs <- WeatherStation$Gr.Rad * 3600 / 1e6
  P <- 101.3*((293-0.0065*WeatherStation$elev)/293)^5.26
  psi <- 0.000665*P
  Delta <- 2503 * exp((17.27*WeatherStation$Ta)/
                        (WeatherStation$Ta+237.3))/((WeatherStation$Ta+237.3)^2)
  ea.sat <- 0.6108*exp((17.27*WeatherStation$Ta)/(WeatherStation$Ta+237.3))
  ea <- (WeatherStation$RH/100)*ea.sat
  DPV <- ea.sat - ea
  dr <- 1 + 0.033*(cos(2*pi*DOY/365))
  delta <- 0.409*sin((2*pi*DOY/365)-1.39)
  phi <- WeatherStation$lat*(pi/180)
  b <- 2*pi*(DOY-81)/364
  Sc <- 0.1645*sin(2*b)-0.1255*cos(b)-0.025*sin(b)
  hour.angle <- (pi/12)*((hours+0.06667*(WeatherStation$long*pi/180-long.z*pi/180)+Sc)-12)
  w1 <- hour.angle-((pi)/24)
  w2 <- hour.angle+((pi)/24)
  hour.angle.s <- acos(-tan(phi)*tan(delta))
  w1c <- ifelse(w1< -hour.angle.s, -hour.angle.s, 
                ifelse(w1>hour.angle.s, hour.angle.s, ifelse(w1>w2, w2, w1)))
  w2c <- ifelse(w2< -hour.angle.s, -hour.angle.s, 
                ifelse(w2>hour.angle.s, hour.angle.s, w2))
  Beta <- asin((sin(phi)*sin(delta)+cos(phi)*cos(delta)*cos(hour.angle)))
  Ra <- ifelse(Beta <= 0, 1e-45, ((12/pi)*4.92*dr)*
                 (((w2c-w1c)*sin(phi)*sin(delta))+(cos(phi)*cos(delta)*(sin(w2)-sin(w1)))))
  Rso <- (0.75+2e-5*WeatherStation$elev)*Ra
  Rs.Rso <- ifelse(Rs/Rso<=0.3, 0, ifelse(Rs/Rso>=1, 1, Rs/Rso))
  fcd <- ifelse(1.35*Rs.Rso-0.35<=0.05, 0.05, 
                ifelse(1.35*Rs.Rso-0.35<1, 1.35*Rs.Rso-0.35,1))
  Rn.a <- ((1-0.23)*Rs) - (2.042e-10*fcd*(0.34-0.14*(ea^0.5))*TaK^4)
  G.day <- Rn.a * 0.1
  wind.2 <- WeatherStation$wind *(4.87/(log(67.8*WeatherStation$height-5.42)))
  ETo.hourly <- ((0.408*Delta*(Rn.a-G.day))+(psi*(37/TaK)*wind.2*(DPV)))/
    (Delta+(psi*(1+(0.24*wind.2))))
  return(ETo.hourly)
}


#########################################################################

