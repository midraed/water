#' Plot method for waterLSEB S3 class
#' @param x       waterLSEB object.
#' @param ...     further arguments passed to or from other methods. 
#' @author Guillermo Federico Olmedo
#' @export
#' @family LSEB objects related functions
#' @method plot waterLSEB
plot.waterLSEB <- function(x, ...){
  plot(x$EB[[1:4]])
}


#' Print method for waterLSEB S3 class
#' @param x        waterLSEB object.
#' @param ...      further arguments passed to or from other methods.
#' @author Guillermo Federico Olmedo
#' @export
#' @family LSEB objects related functions
#' @method print waterLSEB
print.waterLSEB <- function(x, ...){
  print(x$EB, ...)
  print(x$WeatherStation, ...)
  print(x$anchors, ...)
  print(x$anchors@data, ...)
  print(x$methods, ...)
}


#' writeRaster method for waterLSEB S3 class
#' @param x        waterLSEB object.
#' @param ...      additional parameters to pass to writeRaster()
#' @author  María Victoria Munafó
#' @export
#' @family LSEB objects related functions
#' @method writeRaster waterLSEB

writeRaster.waterLSEB <- function(x, ...){
  raster::writeRaster(x$EB, ...)
}


