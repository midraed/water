# Maybe i can provide some points and use CHull like in QGIS-Geostat
create.aoi <- function(topleft = c(x, y), bottomright= c(x, y)){
  asd <- SpatialPolygons(
    list(Polygons(list(Polygon(coords = matrix(
      c(topleft[1],bottomright[1], bottomright[1],topleft[1],topleft[1],
        topleft[2], topleft[2], bottomright[2], 
        bottomright[2],topleft[2]), ncol=2, nrow= 5))), ID=1)))
  return(asd)
}

# Better get the images from a folder...! (like albedo)
load_L8data <-  function(landsat.band, aoi=NULL){
  require(raster)
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
    return(raw.image)}
  else
    raw.image <- crop(raw.image, aoi)
    print("Bands 2,3,4,5,6,7 loaded successfully and cropped with aoi")
    return(raw.image)
   }
  else
  print("ERROR: I expected something like LC82320832013319LGN00_BX.TIF")
  return(NULL)
}

# Get links or optionally open web pages... 
checkSRTMgrids <-function(raw.image, path = getwd(), format="tif"){
  require(raster)
  require(proj4)
  asd <- SpatialPolygons(
    list(Polygons(list(Polygon(coords = matrix(
      c(xmin(raw.image), xmax(raw.image), xmax(raw.image),
        xmin(raw.image),xmin(raw.image), ymax(raw.image),
        ymax(raw.image), ymin(raw.image), ymin(raw.image),
        ymax(raw.image)), ncol=2, nrow= 5))), ID=1)))
  asd@proj4string <- raw.image@crs
  limits <- project(xy = matrix(asd@bbox, ncol=2, nrow=2), proj = asd@proj4string, 
          inverse = TRUE)
  # I have to improve this. It should ONLY work for west and south coordinates.. maybe
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
prepareSRTMdata <- function(path=getwd(), format="tif", slope=TRUE, aspect=TRUE, extent=raw.image){
  files <- list.files(path= path,  pattern=paste("^[sn]\\d{2}_[we]\\d{3}_1arc_v3.", format, "$", sep="")) 
  stack1 <- list()
  for(i in 1:length(files)){
    stack1[[i]] <- raster(paste(path, files[i], sep=""))}
  stack1$fun <- mean
  SRTMmosaic <- do.call(mosaic, stack1)
  destino  <-  projectExtent(raw.image, raw.image@crs)
  mosaicp <- projectRaster(SRTMmosaic, destino)
  aspect <- terrain(mosaicp, opt="aspect") 
  slope <- terrain(mosaicp, opt="slope") 
  aspect_metric <- aspect-pi  #METRIC expects aspect - 1 pi
  surface.model <- stack(mosaicp, slope, aspect_metric)
  names(surface.model) <- c("DEM", "Slope", "Aspect")
  removeTmpFiles(h=0)
  return(surface.model)
}

solar.angles <- function(L8MTL, raw.image, slope, aspect){
   require(stringr)
   test <- scan(L8MTL, character(0), sep = "\n")
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
#  hour.angle2 <- asin((sin(sun.elevation)-(sin(declination)*sin(latitude)))
#         /(cos(declination)*cos(latitude)))
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
   names(solar.angles) <- c("latitude", "declination", "hour.angle", "incidence.hor", "incidence.rel") 
   return(solar.angles)
}

sw.trasmisivity <- function(Kt = 1, ea, dem, incidence.hor){
  P <- 101.3*((293-0.0065 * dem)/293)^5.26
  W <- 0.14 * ea * P + 2.1
  0.35 + 0.627 * exp((-0.00149 * P / Kt * cos(incidence.hor))-0.075*(W / cos(incidence.hor))^0.4)
}

incoming.solar.radiation <- function(incidence.rel, tau.sw, DOY){
  d <- sqrt(1/(1+0.033*cos(DOY * 2 * pi/365)))
  1367 * cos(incidence.rel) * tau.sw / d^2
}

albedo <- function(path=getwd(), aoi= NULL){
    srb2 <- raster(paste(path, list.files(path = path, pattern = "_sr_band2.tif"), sep=""))
    srb3 <- raster(paste(path, list.files(path = path, pattern = "_sr_band3.tif"), sep=""))
    srb4 <- raster(paste(path, list.files(path = path, pattern = "_sr_band4.tif"), sep=""))
    srb5 <- raster(paste(path, list.files(path = path, pattern = "_sr_band5.tif"), sep=""))
    srb6 <- raster(paste(path, list.files(path = path, pattern = "_sr_band6.tif"), sep=""))
    srb7 <- raster(paste(path, list.files(path = path, pattern = "_sr_band7.tif"), sep=""))
    l8.albedo <-  stack(srb2, srb3, srb4, srb5, srb6, srb7)
    if(!missing(aoi)){
      l8.albedo <- crop(l8.albedo,aoi) # Without aoi this should fail on most computers.
      }                                # I have to read every row or chunk and then do the calc
    values(l8.albedo) <- values(l8.albedo[[1]])*0.0001*0.254+values(l8.albedo[[2]])*0.0001*0.149+values(l8.albedo[[3]])*0.0001*0.147+
                         values(l8.albedo[[4]])*0.0001*0.311+values(l8.albedo[[5]])*0.0001*0.103+values(l8.albedo[[6]])*0.0001*0.036
    return(l8.albedo[[1]])
}

LAI.metric <- function(path=getwd(), L=0.1, aoi=NULL){
  toa.red <- raster(paste(path, list.files(path = path, pattern = "_sr_band4.tif"), sep=""))
  toa.nir <- raster(paste(path, list.files(path = path, pattern = "_sr_band5.tif"), sep=""))
  toa.4.5 <- stack(toa.red, toa.nir)
  if(!missing(aoi)){
    toa.4.5 <- crop(toa.4.5,aoi) # Without aoi this should fail on most computers.
  }    
  SAVI_ID <- (1 + L)*(toa.4.5[[2]]*0.0001 - toa.4.5[[1]]*0.0001)/(L + toa.4.5[[1]]*0.0001 + toa.4.5[[2]]*0.0001)
  LAI <- log((0.69-SAVI_ID)/0.59)/0.91 *-1
  removeTmpFiles(h=0)
  return(LAI)
}

outgoing.lw.radiation <- function(path=getwd(), LAI, aoi=NULL){
  bright.temp.b10 <- raster(paste(path, list.files(path = path, pattern = "_toa_band10.tif"), sep=""))
  bright.temp.b11 <- raster(paste(path, list.files(path = path, pattern = "_toa_band11.tif"), sep=""))
  if(!missing(aoi)){
    bright.temp.b10 <- crop(bright.temp.b10,aoi) 
    bright.temp.b11 <- crop(bright.temp.b11,aoi) 
  }
  Ts <- mean(bright.temp.b10, bright.temp.b11)*0.1 ## this isn't split-window corrections
  surf.emissivity <- 0.95 + 0.01 * LAI ## And when LAI 3 or more = 0.98
  return(surf.emissivity * 5.67e-8 * Ts^4)
}

# Add a function to get "cold" pixel temperature so in can be used in the next function
incoming.lw.radiation <- function(air.temperature, DEM, sw.trasmisivity, aoi=NULL){
  Ta <-  air.temperature - (DEM - 702) * 6.49 / 1000 ## Temperature in Kelvin
    if(!missing(aoi)){
    Ta <- crop(Ta,aoi) 
  }
  plot(Ta, main="near surface air temperature in K")
  ef.atm.emissivity <- epsilon_a <- 0.85*(-1*log(sw.trasmisivity))^0.09
  return(ef.atm.emissivity * 5.67e-8 * Ta^4)
}

soil.heat.flux1 <- function(path=getwd(), albedo, Rn, aoi=NULL){
  bright.temp.b10 <- raster(paste(path, list.files(path = path, pattern = "_toa_band10.tif"), sep=""))
  bright.temp.b11 <- raster(paste(path, list.files(path = path, pattern = "_toa_band11.tif"), sep=""))
  toa.red <- raster(paste(path, list.files(path = path, pattern = "_sr_band4.tif"), sep=""))
  toa.nir <- raster(paste(path, list.files(path = path, pattern = "_sr_band5.tif"), sep=""))
  toa.4.5 <- stack(toa.red*0.0001, toa.nir*0.0001) ## chuncks
  if(!missing(aoi)){
    bright.temp.b10 <- crop(bright.temp.b10,aoi) 
    bright.temp.b11 <- crop(bright.temp.b11,aoi)
    toa.4.5 <- crop(toa.4.5,aoi) # Without aoi this should fail on most computers.
  }
  NDVI <- (toa.4.5[[2]] - toa.4.5[[1]])/(toa.4.5[[1]] + toa.4.5[[2]])
  Ts <- mean(bright.temp.b10, bright.temp.b11)*0.1 ## this isn't split-window corrections
  G <- ((Ts - 273.15)*(0.0038+0.0074*albedo)*(1-0.98*NDVI^4))*Rn
  removeTmpFiles(h=0)
  return(G)
}

## Create a function to estimate a and b coefficients or the function between Z.om and NDVI
## Also I should add the function for olives sparse trees used by Santos 2012
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
  return(Z.om)
}

aerodynamic.transport <- function(Z.om, wind, height.ws=2, Z.omw = 0.0018, z1=0.1, z2=2, mountainous=FALSE, surf.model){
    u200 <- wind * log(200/Z.omw)/log(height.ws/Z.omw)
    if(mountainous==TRUE){
      u200 <- u200 * (1+0.1*(surf.model$DEM-height.ws))
    }
    friction.velocity <- 0.41 * u200 / log(200/Z.om)
    r.ah <- log(z2/z1)/friction.velocity*0.41
}

hot.and.cold <- function(method="random", n=1, path=getwd(), ETr, Rn, G, r.ah, DEM, LAI, aoi){
  bright.temp.b10 <- raster(paste(path, list.files(path = path, pattern = "_toa_band10.tif"), sep=""))
  bright.temp.b11 <- raster(paste(path, list.files(path = path, pattern = "_toa_band11.tif"), sep=""))
  if(!missing(aoi)){
    bright.temp.b10 <- crop(bright.temp.b10,aoi) 
    bright.temp.b11 <- crop(bright.temp.b11,aoi) 
  }
  Ts <- mean(bright.temp.b10, bright.temp.b11)*0.1
  Ts_datum <- Ts - (DEM - 702) * 6.49 / 1000
  P <- 101.3*((293-0.0065 * DEM)/293)^5.26
  air.density <- 1000 * P / (1.01*(Ts)*287)
  if(method=="random"){
    if(!max(values(Ts))>=310){
      warning(paste("Ts max value is", round(max(values(Ts)),2), "and I expected Ts>=310"))
      hot <- sample(which(values(Ts)==max(values(Ts))),1)  
    } else hot <- sample(which(values(Ts>310)),n) 
    if(!max(values(LAI))>=4){
      warning(paste("LAI max value is", round(max(values(LAI)),2), "and I expected LAI>=4"))
      cold <- sample(which(values(LAI)==max(values(LAI))),1)  
    } else cold <- sample(which(values(LAI>4)),n)  
  }
  dT_hot <- (Rn[hot] - G[hot])*r.ah[hot]/(air.density[hot]*1007)
  lambda <- (2.501-0.00236*(Ts-273.15))  # En el paper dice por 1e6
  LE_cold <- 1.05 * ETr * lambda[cold]
  H_cold <- Rn[cold]-G[cold]-LE_cold
  dT_cold <- H_cold*r.ah[cold]/(air.density[cold]*1007)
  a <- (dT_hot-dT_cold)/(Ts_datum[hot]-Ts_datum[cold])
  b <- (dT_hot-a)/Ts_datum[hot]
  dT <- a+b*Ts_datum
  # Prepare and populate result list
  result <- list()
  result$dT <- dT
  result$hot.and.cold <- data.frame(pixel=integer(), Ts=double(), 
                                    LAI=double(), factor(levels = c("hot", "cold")))
  for(i in 1:n){result$hot.and.cold[i, ] <- c(hot, Ts[hot], LAI[hot], "hot")}
  for(i in 1:n){result$hot.and.cold[i+n, ] <- c(cold, Ts[cold], LAI[cold], "cold")}
  #removeTmpFiles(h=0)
  return(result)
}






### Add a function to estimate all parameters for some points

#########################################################################
save(create.aoi, load_L8data, checkSRTMgrids, prepareSRTMdata, solar.angles,
     sw.trasmisivity, incoming.solar.radiation, albedo,
     LAI.metric, outgoing.lw.radiation, incoming.lw.radiation,
     soil.heat.flux1, momentum.roughness.length, aerodynamic.transport,
     hot.and.cold,
     file="L8METRICforR/MfR_functions.RData")
