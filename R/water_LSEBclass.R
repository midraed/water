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
  print(x$methods, ...)
  invisible(x)
}


#' writeRaster method for waterLSEB objects
#' @description Writes the energy balance layers (\code{x$EB}) of a
#' METRIC.EB() result to a file.
#' @param x        waterLSEB object, from METRIC.EB().
#' @param filename output file name, e.g. "EB.tif"
#' @param ...      additional arguments passed to terra::writeRaster(), e.g.
#' overwrite
#' @author  María Victoria Munafó
#' @family LSEB objects related functions
#' @name writeRaster.waterLSEB
#' @aliases writeRaster,waterLSEB,character-method
#' @exportMethod writeRaster
#' @importFrom methods setOldClass setMethod
NULL

setOldClass("waterLSEB")

setMethod("writeRaster", signature(x = "waterLSEB", filename = "character"),
          function(x, filename, ...){
            invisible(terra::writeRaster(x$EB, filename, ...))
          })


