#' Plot method for waterLSEB S3 class
#' @param x       waterLSEB object. 
#' @author Guillermo Federico Olmedo
#' @export
#' @family LSEB objects related functions
#' @method plot waterLSEB
plot.waterLSEB <- function(x){
  plot(x$EB[[1:4]])
}


#' Print method for waterLSEB S3 class
#' @param x        waterLSEB object.
#' @author Guillermo Federico Olmedo
#' @export
#' @family LSEB objects related functions
#' @method print waterLSEB
print.waterLSEB <- function(x){
  print(x$EB)
  print(x$WeatherStation)
  print(x$anchors)
  print(x$anchors@data)
  print(x$methods)
}


