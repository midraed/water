.onLoad <- function(libname, pkgname) {
  op <- options()
  op.water <- list(
    waterOverwrite = TRUE,
    waterWriteResults = FALSE,
    waterOutputFolder = ".",
    waterSRTMrepo = NULL
  )
  toset <- !(names(op.water) %in% names(op))
  if(any(toset)) options(op.water[toset])
  invisible()
}



#' Global options for water package
#' @description
#' Sets or shows the global options of the water package. Inspired by
#' rasterOptions() of the raster package (Robert J. Hijmans).
#' @param overwrite     Logical. If TRUE and writeResults is TRUE, the final
#' outputs overwrite files of earlier sessions; a file written earlier in the
#' same session is never overwritten (an object may still read it): a new name
#' is used, with a warning (so, in one session, EB.tif keeps the first result). If FALSE, a time stamp is appended to the file name
#' (name_YYYYmmddHHMMSS.tif).
#' @param writeResults  Logical. If TRUE (default FALSE) results are written to
#' disk as GeoTIFF files: the final outputs (EB from METRIC.EB(), ET24 from
#' ET24h(), DEM from getDEM() and prepareSRTMdata()) to outputFolder, and the
#' intermediate results to terra's temporary files, with unique names. This is
#' slower, but avoids out-of-memory problems with large images. Values are
#' stored as 32-bit floats. The temporary files are deleted when the R session
#' ends (they can be large for full scenes).
#' @param outputFolder  Name of a folder to save files, relative to the working
#' directory.
#' @param SRTMrepo      A folder where SRTM grids are stored. Currently not
#' used by prepareSRTMdata(); use its path argument
#' @param autoAoi       Deprecated and ignored, with a warning. An object
#' called aoi in the workspace is no longer used: pass aoi to loadImage(),
#' loadImageSR(), albedo() or the METRIC functions. See \code{\link{createAoi}}
#' @param default       Logical. If TRUE will revert all options to default
#' values
#' @return
#' list of the current options (invisibly). If no arguments are provided the options are printed.
#' @family support functions
#' @references
#' Robert J. Hijmans (2015). raster: Geographic Data Analysis and Modeling. R
#' package version 2.4-18. https://CRAN.R-project.org/package=raster
#' @author Guillermo Federico Olmedo
#' @export
#' @examples
#' # show the options
#' waterOptions()
waterOptions <- function (overwrite, writeResults, outputFolder,
                          SRTMrepo, autoAoi, default = FALSE) 
{
  setOverwrite <- function(overwrite) {
    if (is.logical(overwrite)) {
      options(waterOverwrite = overwrite)
    }
    else {
      warning(paste("Could not set overwrite. It must be a logical value"))
    }
  }
  setWriteResults <- function(write) {
    if (is.logical(write)) {
      options(waterWriteResults = write)
    }
    else {
      warning(paste("Could not set writeResults. It must be a logical value"))
    }
  }
  setOutputFolder <- function(tmpdir) {
    if (!missing(tmpdir)) {
      tmpdir <- trimws(tmpdir)
      if (tmpdir != "") {
        lastchar = substr(tmpdir, nchar(tmpdir), nchar(tmpdir))
        if (lastchar != "/" & lastchar != "\\") {
          tmpdir <- paste0(tmpdir, "/")
        }
        options(waterOutputFolder = tmpdir)
      }
    }
  }
  setSRTMrepo <- function(srtmDir) {
    if (!missing(srtmDir)) {
      srtmDir <- trimws(srtmDir)
      if (srtmDir != "") {
        lastchar = substr(srtmDir, nchar(srtmDir), nchar(srtmDir))
        if (lastchar != "/" & lastchar != "\\") {
          srtmDir <- paste0(srtmDir, "/")
        }
        options(waterSRTMrepo = srtmDir)
      }
    }
  }
  cnt <- 0
  if (default) {
    cnt <- 1
    options(waterOverwrite = TRUE)
    options(waterWriteResults = FALSE)
    options(waterOutputFolder = ".")
    options(waterSRTMrepo = NULL)
  }
  if (!missing(overwrite)) {
    setOverwrite(overwrite)
    cnt <- cnt + 1
  } else {overwrite <- getOption("waterOverwrite")}
  if (!missing(writeResults)) {
    setWriteResults(writeResults)
    cnt <- cnt + 1
  } else {writeResults <- getOption("waterWriteResults")}
  if (!missing(outputFolder)) {
    setOutputFolder(outputFolder)
    cnt <- cnt + 1
  } else {outputFolder <- getOption("waterOutputFolder")}
  if (!missing(SRTMrepo)) {
    setSRTMrepo(SRTMrepo)
    cnt <- cnt + 1
  } else {SRTMrepo <- getOption("waterSRTMrepo")}
  if (!missing(autoAoi)) {
    warning("waterOptions(autoAoi) is no longer used: an object called aoi in ",
            "the workspace is never used. Pass aoi to loadImage(), ",
            "loadImageSR(), albedo() or the METRIC functions")
    cnt <- cnt + 1
  }
  lst <- list(overwrite = overwrite, writeResults = writeResults,
              outputFolder = outputFolder, SRTMrepo= SRTMrepo)
  if (cnt == 0) {
    cat("overwrite     :", lst$overwrite, "\n")
    cat("writeResults  :", lst$writeResults, "\n")
    cat("outputFolder  :", lst$outputFolder, "\n")
    cat("SRTMrepo      :", lst$SRTMrepo, "\n")
  }
  invisible(lst)
}


