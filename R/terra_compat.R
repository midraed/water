# Compatibilidad mínima para migrar de raster a terra

raster <- function(x, ...) {
  if (missing(x)) return(terra::rast())
  if (inherits(x, "SpatRaster")) return(x)
  terra::rast(x, ...)
}

stack <- function(...) {
  args <- list(...)
  if (length(args) == 1 && is.list(args[[1]])) {
    args <- args[[1]]
  }
  rasters <- lapply(args, terra::rast)
  do.call(terra::c, rasters)
}

brick <- stack

extent <- function(x) {
  terra::ext(x)
}

projectExtent <- function(x, crs) {
  terra::rast(ext = terra::ext(x), crs = crs)
}

projectRaster <- function(x, y) {
  terra::project(x, y)
}

calc <- function(x, fun, ...) {
  terra::app(x, fun = fun, ...)
}

cellStats <- function(x, stat, ...) {
  as.numeric(terra::global(x, fun = stat, ...)[1, 1])
}

nlayers <- function(x) {
  terra::nlyr(x)
}

removeTmpFiles <- function(h = 0) {
  terra::tmpFiles(old = TRUE, remove = TRUE)
  invisible(TRUE)
}

rasterOptions <- function(...) {
  invisible(list(...))
}

raster_which_min <- function(x) {
  as.integer(terra::which.minmax(x, "min"))
}

raster_which_max <- function(x) {
  as.integer(terra::which.minmax(x, "max"))
}
