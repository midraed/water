get.sat <- function(path){
  band <- substr(list.files(path=path, pattern=paste0("^L[EC]\\d+\\w+\\d+_B2.TIF")), 0,3)
  if(length(band)==0){
    print(paste("ERROR: I expected something like landsat.band = LC82320832013319LGN00_BX.TIF in ", path))
    return()
  }
  if(band =="LC8"){return("L8")}
  if(band =="LE7"){return("L7")}
  if(band != "LE7" & band != "LC8"){
    print("Can't establish sat from band names")
    return()}
}


save.load.clean <- function(imagestack, stack.names=NULL, file, ...){
  if(missing(file)){file <- paste(deparse(substitute(raw.image)),".tif", sep="")}
  writeRaster(imagestack, filename = file, ...)
  stack <- stack(file)
  names(stack) <- stack.names
  removeTmpFiles(h=0)
  return(stack)
}


#' Calculates short wave transmisivity
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
sw.trasmisivity <- function(Kt = 1, ea, dem, incidence.hor, result.folder=NULL){
  P <- 101.3*((293-0.0065 * dem)/293)^5.26
  W <- 0.14 * ea * P + 2.1
  sw.t <- 0.35 + 0.627 * exp((-0.00149 * P / Kt * 
                                cos(incidence.hor))-0.075*(W / cos(incidence.hor))^0.4)
  return(sw.t)
}

aoi.crop <- function(raster, aoi){
  if(!missing(aoi)){
    raster <- crop(raster,aoi)
    return(raster)
  }
  if(missing(aoi) & exists(x = "aoi", envir=.GlobalEnv)){
    aoi <- get(x = "aoi", envir=.GlobalEnv)
    raster <- crop(raster,aoi)
    return(raster)
  }
  return(raster)
}
