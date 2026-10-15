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


#' Example files of the water package
#'
#' Files in \code{system.file("extdata", package = "water")}, used in the
#' examples and vignettes:
#' \describe{
#'   \item{apples.csv}{Weather station in an apple orchard near Talca, Chile
#'   (latitude -35.42222, longitude -71.38639, 201 m, sensors at 2.2 m), with
#'   records every 15 minutes on 15 February 2013, in local time
#'   (America/Santiago). Columns: Date (dd/mm/yyyy), Time, Rad (solar
#'   radiation, W/m2), wind_speed (km/h), wind_dir (degrees), RH (\%), temp
#'   (Celsius) and pp (rain, mm). Read it with \code{cf = c(1, 0.2777778, 1)}
#'   to convert the wind speed to m/s.}
#'   \item{INTA.csv}{Hourly weather station records near Mendoza, Argentina
#'   (latitude -33.00513, longitude -68.86469, 927 m, sensors at 2 m) on 9
#'   February 2016, in local time (America/Argentina/Mendoza). Columns:
#'   datetime (yyyy/mm/dd HH:MM), temp (Celsius), RH (\%), pp (rain, mm),
#'   radiation (W/m2) and wind (m/s).}
#'   \item{L7.MTL.txt}{Metadata of the Landsat 7 scene of \code{L7_Talca}
#'   (15 February 2013).}
#'   \item{LC82320832016040LGN00_band2.tif to _band11.tif, LC82320832016040LGN00_MTL.txt}{Landsat 8
#'   Level-1 digital numbers (bands 2 to 7, 10 and 11, with pre-Collection
#'   file names) of path 232, row 83, 9 February 2016, near Mendoza, and its
#'   metadata. Read them with \code{loadImage(sat = "L8")}.}
#'   \item{L8_C2}{Folder with the Collection 2 Level-2 product of the same scene,
#'   on the same grid: surface reflectance (SR_B2 to SR_B7), surface
#'   temperature (ST_B10), quality band (QA_PIXEL) and metadata (MTL). Read it
#'   with \code{loadImageSR()} and mask clouds with \code{cfmask()}.}
#' }
#' @name water-extdata
#' @aliases extdata apples.csv INTA.csv
#' @keywords datasets
#' @examples
#' list.files(system.file("extdata", package = "water"))
NULL