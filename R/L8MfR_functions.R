### ADD writeRaster to all the functions
### Add a function to estimate all parameters for some points (and only for those points!)
### Add some meta-functions to insolate processes
### Add METRIC function
### Maybe we need a new class to store all data? 
### Maybe a class for weather stations... check previous work on the field
### Maybe add knitr and generate pdf report on METRIC function
### Add Perrier equation for zom
### Maybe change from this.names to thisNames
### Export anchor as kml or view in Google
### select anchor points with buffer for not too close pixels...
### Select anchor for multiple criteria with table from Marcos
### Check Rsky on METRIC 2010
### Add dependency with evapotranspiration and used to measure ETp
### Add three source temperature model..!


#################################################################################3
# Common fields for docs
# @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
# @references 
# R. G. Allen, M. Tasumi, and R. Trezza, “Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model,” Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007


####################

save.load.clean <- function(imagestack, stack.names=NULL, file,  ...){
  if(missing(file)){file <- paste(deparse(substitute(raw.image)),".tif", sep="")}
  writeRaster(imagestack, filename = file, ...)
  stack <- stack(file)
  names(stack) <- stack.names
  removeTmpFiles(h=0)
  return(stack)
}


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
#' @param landsat.band one of the landsat 8 bands to import
#' @param aoi area of interest to crop the raster after loading
#' @return object of class rasterStack
#' @exportClass object of class RasterStack
#' @export
# Better get the images from a folder...! (like albedo)
# Check when creates temp files... on raster() or in stack()
load_L8data <-  function(landsat.band, aoi){
  if(grepl("[LC[:digit:]]+LGN00_B", landsat.band, perl=TRUE) == TRUE){
  B2 <- raster(sub("B[[:digit:]]", "B2", landsat.band))
  B3 <- raster(sub("B[[:digit:]]", "B3", landsat.band))
  B4 <- raster(sub("B[[:digit:]]", "B4", landsat.band))
  B5 <- raster(sub("B[[:digit:]]", "B5", landsat.band))
  B6 <- raster(sub("B[[:digit:]]", "B6", landsat.band))
  B7 <- raster(sub("B[[:digit:]]", "B7", landsat.band))
  raw.image <- stack(B2,B3,B4,B5,B6,B7)
  if(missing(aoi)){
    print("Bands 2,3,4,5,6,7 full scene loaded successfully")
    raw.image <- save.load.clean(imagestack = raw.image, file = "raw.image.tif", overwrite=TRUE)
    return(raw.image)}
  else
    raw.image <- crop(raw.image, aoi)
    print("Bands 2,3,4,5,6,7 loaded successfully and cropped with aoi")
    #plotRGB(raw.image, 3,2,1, stretch="lin",  main="RGB 4,3,2")
    raw.image <- save.load.clean(imagestack = raw.image, stack.names = c("B", "G", "R", "NIR", "SWIR1", "SWIR2"), file = "raw.image.tif", overwrite=TRUE)
    return(raw.image)
   }
  else
  print("ERROR: I expected something like landsat.band = LC82320832013319LGN00_BX.TIF")
  return(NULL)
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
prepareSRTMdata <- function(path=getwd(), format="tif", extent=raw.image){
  files <- list.files(path= path,  pattern=paste("^[sn]\\d{2}_[we]\\d{3}_1arc_v3.", format, "$", sep="")) 
  stack1 <- list()
  for(i in 1:length(files)){
    stack1[[i]] <- raster(paste(path, files[i], sep=""))}
  stack1$fun <- mean
  SRTMmosaic <- do.call(mosaic, stack1)
  destino  <-  projectExtent(raw.image, raw.image@crs)
  mosaicp <- projectRaster(SRTMmosaic, destino)
  mosaicp <- save.load.clean(imagestack = mosaicp, stack.names = "DEM", file = "DEM.tif", overwrite=TRUE)
  return(mosaicp)
}

METRIC.topo <- function(DEM){
  DEM <- raster(DEM)
  aspect <- terrain(DEM, opt="aspect") 
  slope <- terrain(DEM, opt="slope") 
  aspect_metric <- aspect-pi  #METRIC expects aspect - 1 pi
  surface.model <- stack(DEM, slope, aspect_metric)
  surface.model <- save.load.clean(imagestack = surface.model, stack.names = c("DEM", "Slope", "Aspect"), file = "surface.model.tif", overwrite=TRUE)
  return(surface.model)
}

solar.angles <- function(L8MTL, raw.image, slope, aspect){
   test <- scan(L8MTL, character(0), sep = "\n") ## Here uses stringr,.. only 1 line.. doesn't justify adding a dependency!
   sun.azimuth <- as.numeric(str_extract(test[68], pattern = "([0-9]{1,5})([.]+)([0-9]+)"))*pi/180
   sun.elevation <- as.numeric(str_extract(test[69], pattern = "([0-9]{1,5})([.]+)([0-9]+)"))*pi/180
   # latitude
   latitude <- raw.image[[1]]
   xy <- SpatialPoints(xyFromCell(latitude, cellFromRowCol(latitude, 1:nrow(latitude), 1)))
   xy@proj4string <- latitude@crs
   lat <- coordinates( spTransform(xy, CRS("+proj=longlat +datum=WGS84")))[,2] 
   values(latitude) <- rep(lat*pi/180,each=ncol(latitude))
   # declination
   DOY <- strptime(str_extract(test[21], pattern = "([0-9]{4})([-]+)([0-9]{2})([-]+)([0-9]{2})"), "%Y-%m-%d")$yday+1
   declination <- raw.image[[1]]
   values(declination) <- 23.45*pi/180*sin(2*pi*((284+DOY)/36.25))
   # hour angle
   hour.angle <- asin(-1*(cos(sun.elevation)*sin(sun.azimuth)/cos(declination)))
   ## solar incidence angle, for horizontal surface
   incidence.hor <- acos(sin(declination) * sin(latitude) + cos(declination)*cos(latitude)*cos(hour.angle))
   ##solar incidence angle, for sloping surface
   incidence.rel <- acos(sin(declination)*sin(latitude)*cos(slope) 
                     - sin(declination)*cos(latitude)*sin(slope)*cos(aspect)
                     + cos(declination)*cos(latitude)*cos(slope)*cos(hour.angle)
                     + cos(declination)*sin(latitude)*sin(slope)*cos(aspect)*cos(hour.angle)
                     + cos(declination)*sin(aspect)*sin(slope)*sin(hour.angle))
   ## End
   solar.angles <- stack(latitude, declination, hour.angle, incidence.hor, incidence.rel)
   solar.angles <- save.load.clean(imagestack = solar.angles, stack.names = c("latitude", "declination", "hour.angle", "incidence.hor", "incidence.rel"), file = "solar.angles.tif", overwrite=TRUE)
   return(solar.angles)
}

sw.trasmisivity <- function(Kt = 1, ea, dem, incidence.hor){
  P <- 101.3*((293-0.0065 * dem)/293)^5.26
  W <- 0.14 * ea * P + 2.1
  sw.t <- 0.35 + 0.627 * exp((-0.00149 * P / Kt * cos(incidence.hor))-0.075*(W / cos(incidence.hor))^0.4)
  sw.t <- save.load.clean(imagestack = sw.t, file = "sw.t.tif", overwrite=TRUE)
  return(sw.t)
}

incoming.solar.radiation <- function(incidence.rel, tau.sw, DOY){
  d <- sqrt(1/(1+0.033*cos(DOY * 2 * pi/365)))
  Rs.inc <- 1367 * cos(incidence.rel) * tau.sw / d^2
  Rs.inc <- save.load.clean(imagestack = Rs.inc, file = "Rs.inc.tif", overwrite=TRUE)
  return(Rs.inc)
}

albedo <- function(path=getwd(), aoi, coeff="Tasumi"){
    if(coeff=="Tasumi"){wb <- c(0.254, 0.149, 0.147, 0.311, 0.103, 0.036) * 10000} # Tasumi 2008
    if(coeff=="Olmedo") {wb <- c(0.246, 0.146, 0.191, 0.304, 0.105, 0.008) * 10000 }# Calculated using SMARTS for Kimberly2-noc13 and Direct Normal Irradiance
    if(coeff=="Liang") {wb <- c(0.356, 0, 0.130, 0.373, 0.085, 0.072) * 10000} # Liang 2001
    srb2 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band2.tif"), sep="")[1]), fun=function(x){x *wb[1]})
    srb3 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band3.tif"), sep="")[1]), fun=function(x){x *wb[2]})
    srb4 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band4.tif"), sep="")[1]), fun=function(x){x *wb[3]})
    srb5 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band5.tif"), sep="")[1]), fun=function(x){x *wb[4]})
    srb6 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band6.tif"), sep="")[1]), fun=function(x){x *wb[5]})
    srb7 <- calc(raster(paste(path, list.files(path = path, pattern = "_sr_band7.tif"), sep="")[1]), fun=function(x){x *wb[6]})
    l8.albedo <-  stack(srb2, srb3, srb4, srb5, srb6, srb7)/1e+08
    if(!missing(aoi)){
      l8.albedo <- crop(l8.albedo,aoi) # Without aoi this should fail on most computers.
      }                                
    if(coeff=="Liang"){
      l8.albedo <- l8.albedo - 0.0018
    }
    l8.albedo <- stackApply(l8.albedo, indices = c(1,1,1,1,1,1), fun=sum)
    l8.albedo <- save.load.clean(imagestack = l8.albedo, file = "l8.albedo.tif", overwrite=TRUE)
    return(l8.albedo)
}

LAI.from.L8 <- function(method="metric", path=getwd(), aoi, L=0.1){
  if(method=="metric" | method=="vineyard" | method=="carrasco-benavides"){
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
  LAI <- save.load.clean(imagestack = LAI, stack.names = "LAI", file = "LAI.tif", overwrite=TRUE)
  return(LAI)
}

### Correcto to only B10
surface.temperature <- function(path=getwd(), aoi){
  bright.temp.b10 <- raster(paste(path, list.files(path = path, pattern = "_toa_band10.tif"), sep=""))
  if(!missing(aoi)){
    bright.temp.b10 <- crop(bright.temp.b10,aoi) 
  }
  Ts <- bright.temp.b10*0.1
  Ts <- save.load.clean(imagestack = Ts, file = "Ts.tif", overwrite=TRUE)
  return(Ts)
}

outgoing.lw.radiation <- function(path=getwd(), LAI, aoi){
  Ts <- surface.temperature(path=path, aoi=aoi)
  surf.emissivity <- 0.95 + 0.01 * LAI ## And when LAI 3 or more = 0.98
  Rl.out <- surf.emissivity * 5.67e-8 * Ts^4
  Rl.out <- save.load.clean(imagestack = Rs.out, file = "Rs.out.tif", overwrite=TRUE)
  return(Rl.out)
}

# Add a function to get "cold" pixel temperature so in can be used in the next function
incoming.lw.radiation <- function(air.temperature, DEM, sw.trasmisivity, aoi){
  Ta <-  air.temperature - (DEM - 702) * 6.49 / 1000 ## Temperature in Kelvin
    if(!missing(aoi)){
    Ta <- crop(Ta,aoi) 
  }
  ef.atm.emissivity  <- 0.85*(-1*log(sw.trasmisivity))^0.09
  Rl.in <- ef.atm.emissivity * 5.67e-8 * Ta^4
  Rl.in <- save.load.clean(imagestack = Rl.in, file = "Rl.in.tif", overwrite=TRUE)
  return(Rl.in)
}

soil.heat.flux1 <- function(path=getwd(), albedo, Rn, aoi){
  toa.red <- raster(paste(path, list.files(path = path, pattern = "_sr_band4.tif"), sep="")[1])
  toa.nir <- raster(paste(path, list.files(path = path, pattern = "_sr_band5.tif"), sep="")[1])
  toa.4.5 <- stack(toa.red*0.0001, toa.nir*0.0001) ## chuncks
  if(!missing(aoi)){
    toa.4.5 <- crop(toa.4.5,aoi) # Without aoi this should fail on most computers.
  }
  NDVI <- (toa.4.5[[2]] - toa.4.5[[1]])/(toa.4.5[[1]] + toa.4.5[[2]])
  Ts <- surface.temperature(path=path, aoi=aoi)
  G <- ((Ts - 273.15)*(0.0038+0.0074*albedo)*(1-0.98*NDVI^4))*Rn
  G <- save.load.clean(imagestack = G, file = "G.tif", overwrite=TRUE)
  return(G)
}

## Create a function to estimate a and b coefficients or the function between Z.om and NDVI
## using some points and tabulated z.om for their covers.
## Also I should add the function for olives sparse trees used by Santos 2012
## Or Pocas 2014.
momentum.roughness.length <- function(method, path=getwd(), LAI, NDVI, albedo, a, b,mountainous=FALSE, surf.model){
  if(method=="short.crops"){
    print="using method for short agricultural crops (Tasumi 2003)"
    Z.om <- (0.018*LAI)
  }
  if(method=="custom"){
    Z.om <- exp((a*NDVI/albedo)+b)
  }
  if(mountainous==TRUE){
    Z.om <- Z.om * (1 + (180/pi*surf.model$Slope - 5)/20)
  }
  Z.om <- save.load.clean(imagestack = Z.om, file = "Z.om.tif", overwrite=TRUE)
  return(Z.om)
}

air.density <- function(DEM, Ts){
  P <- 101.3*((293-0.0065 * dem)/293)^5.26
  air.density <- 1000 * P / (1.01*(Ts)*287)
}


#' Calculates sensible heat flux for METRIC Model
#' @param rho.air - air density
#' @param dT - dT from anchor points
#' @param r.ah - rugosity
#' @return Returns sensible heat flux
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
sensible.heat.flux <- function(rho.air, dT, r.ah){
  H <- rho.air*1007*dT/r.ah
  H <- save.load.clean(imagestack = H, file = "H.tif", overwrite=TRUE)
  return(H)
}






#########################################################################

