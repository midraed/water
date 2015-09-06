.onLoad <- function(libname, pkgname) {
  # Wickham it's the man!
  op <- options()
  op.water <- list(
    waterOverwrite = TRUE,
    waterWriteResults = TRUE,
    waterDestFolder = ".",
    waterSRTMrepo = NULL,
    waterAutoAoi = TRUE
  )
  toset <- !(names(op.water) %in% names(op))
  if(any(toset)) options(op.water[toset])
  packageStartupMessage("This package writes function's results to working directory. You can change")
  packageStartupMessage("output folder, or completely disable this feature using waterOptions()")

  invisible()
}

#' Global options for water package
#' @description 
#' This function is based on raster::rasterOptions by Robert Hijmans. 
#' @return 
#' list of the current options (invisibly). If no arguments are provided the options are printed.
#' @references 
#' Robert J. Hijmans (2015). raster: Geographic Data Analysis and Modeling. R
#' package version 2.4-18. http://CRAN.R-project.org/package=raster
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @export
waterOptions <- function (overwrite, writeResults, destinationFolder,
                          SRTMrepo, autoAoi, default = FALSE) 
{
#   setFiletype <- function(format) {
#     if (.isSupportedFormat(format)) {
#       options(rasterFiletype = format)
#     }
#     else {
#       warning(paste("Cannot set filetype to unknown or unsupported file format:", 
#                     format, ". See writeFormats()"))
#     }
#   }
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
  setDestFolder <- function(tmpdir) {
    if (!missing(tmpdir)) {
      tmpdir <- trim(tmpdir)
      if (tmpdir != "") {
        lastchar = substr(tmpdir, nchar(tmpdir), nchar(tmpdir))
        if (lastchar != "/" & lastchar != "\\") {
          tmpdir <- paste0(tmpdir, "/")
        }
        options(waterDestFolder = tmpdir)
      }
    }
  }
  setSRTMrepo <- function(srtmDir) {
    if (!missing(srtmDir)) {
      srtmDir <- trim(srtmDir)
      if (srtmDir != "") {
        lastchar = substr(srtmDir, nchar(srtmDir), nchar(srtmDir))
        if (lastchar != "/" & lastchar != "\\") {
          srtmDir <- paste0(srtmDir, "/")
        }
        options(waterSRTMrepo = srtmDir)
      }
    }
  }
  setAutoAoi <- function(autoAoi) {
    if (is.logical(autoAoi)) {
      options(waterAutoAoi = autoAoi)
    }
    else {
      warning(paste("Could not set autoAoi. It must be a logical value"))
    }
  }
#   setToDisk <- function(todisk) {
#     if (is.logical(todisk)) {
#       options(rasterToDisk = todisk)
#     }
#     else {
#       warning(paste("todisk argument must be a logical value"))
#     }
#   }
#   depracatedWarnings <- function(x) {
#     if (is.logical(x)) {
#       if (is.na(x)) {
#         x <- TRUE
#       }
#       options(rasterDepracatedWarnings = x)
#     }
#   }
  cnt <- 0
  if (default) {
    cnt <- 1
    options(waterOverwrite = TRUE)
    options(waterWriteResults = TRUE)
    options(waterDestFolder = ".")
    options(waterSRTMrepo = NULL)
    options(waterAutoAoi = TRUE)
    v <- utils::packageDescription("water")[["Version"]]
  }
  if (!missing(overwrite)) {
    setOverwrite(overwrite)
    cnt <- cnt + 1
  }
  if (!missing(writeResults)) {
    setWriteResults(writeResults)
    cnt <- cnt + 1
  }
  if (!missing(destinationFolder)) {
    setDestFolder(destinationFolder)
    cnt <- cnt + 1
  }
  if (!missing(SRTMrepo)) {
    setSRTMrepo(SRTMrepo)
    cnt <- cnt + 1
  }
  if (!missing(autoAoi)) {
    setAutoAoi(autoAoi)
    cnt <- cnt + 1
  }
  
  ## Continue here:
  
  lst <- list(format = .filetype(), overwrite = .overwrite(), 
              datatype = .datatype(), tmpdir = tmpDir(create = FALSE), 
              tmptime = .tmptime(), progress = .progress(), timer = .timer(), 
              chunksize = .chunksize(), maxmemory = .maxmemory(), 
              todisk = .toDisk(), setfileext = .setfileext(), tolerance = .tolerance(), 
              standardnames = .standardnames(), depwarning = .depracatedwarnings(), 
              addheader = .addHeader())
  save <- FALSE
  if (save) {
    v <- utils::packageDescription("water")[["Version"]]
    fn <- paste(options("startup.working.directory"), "/rasterOptions_", 
                v, sep = "")
    oplst <- NULL
    oplst <- c(oplst, paste("rasterFiletype='", lst$format, 
                            "'", sep = ""))
    oplst <- c(oplst, paste("rasterOverwrite=", lst$overwrite, 
                            sep = ""))
    oplst <- c(oplst, paste("rasterDatatype='", lst$datatype, 
                            "'", sep = ""))
    oplst <- c(oplst, paste("rasterTmpDir='", lst$tmpdir, 
                            "'", sep = ""))
    oplst <- c(oplst, paste("rasterTmpTime='", lst$tmptime, 
                            "'", sep = ""))
    r <- try(write(unlist(oplst), fn), silent = TRUE)
    cnt <- 1
  }
  #overwrite, writeResults, destinationFolder, SRTMrepo, autoAoi
  if (cnt == 0) {
    cat("overwrite     :", lst$overwrite, "\n")
    cat("writeResults  :", lst$writeResults, "\n")
    cat("destFolder    :", lst$destinationFolder, "\n")
    cat("SRTMrepo      :", lst$SRTMrepo, "\n")
    cat("autoAoi       :", lst$autoAoi, "\n")
  }
  invisible(lst)
}
