#' Landsat 7 scene from east of Talca, Chile
#'
#' A Landsat 7 ETM+ subset (Level-1, digital numbers) stored as a terra
#' \code{PackedSpatRaster}. Unpack it with \code{rast(L7_Talca)}. Layers are, in
#' this order, B, G, R, NIR, SWIR1, Thermal1 (band 6, low gain) and SWIR2
#' (calcTOAr() and calcRadiance() select the bands by name).
#'
#' Metadata for this image is also provided as system file. You can
#' load it with:  system.file("extdata", "L7.MTL.txt", package="water")
#'
#' Data available from the U.S. Geological Survey.
#'
#' @format A \code{PackedSpatRaster}; \code{rast(L7_Talca)} is a SpatRaster with
#' 417 rows, 508 columns and 7 layers, 30 m resolution, WGS 84 / UTM zone 19S.
#' @source \url{https://www.usgs.gov/}
#' @examples
#' image.DN <- rast(L7_Talca)
#' names(image.DN)
"L7_Talca"


#' SRTM DEM from east of Talca, Chile
#'
#' A SRTM digital elevation model (m) on the same grid as \code{L7_Talca},
#' stored as a terra \code{PackedSpatRaster}. Unpack it with
#' \code{rast(DEM_Talca)}.
#'
#' Original data is float. Example data is integer.
#'
#' Data available from the U.S. Geological Survey.
#'
#' @format A \code{PackedSpatRaster}; \code{rast(DEM_Talca)} is a SpatRaster with
#' 417 rows, 508 columns and one layer (DEM, 131 to 643 m), WGS 84 / UTM zone 19S.
#' @source \url{https://www.usgs.gov/}
#' @examples
#' DEM <- rast(DEM_Talca)
#' plot(DEM)
"DEM_Talca"