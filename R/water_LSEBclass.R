#' Plot method for waterLSEB S3 class
#' @param x       waterLSEB object. 
#' @author Guillermo Federico Olmedo
#' @export
#' @family LSEB objects related functions
#' @method plot waterLSEB
plot.waterLESB <- function(x){
  plot(x$EB[[1:4]])
}


#' Print method for waterLSEB S3 class
#' @param x        waterLSEB object.
#' @author Guillermo Federico Olmedo
#' @export
#' @family LSEB objects related functions
#' @method print waterLSEB
print.waterWeatherStation <- function(x){
  print(x)
}


