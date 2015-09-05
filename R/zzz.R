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
waterOptions <- function (format, overwrite, writeResults, destinationFolder,
                          SRTMrepo, autoAoi, default = FALSE) 
{
  setFiletype <- function(format) {
    if (.isSupportedFormat(format)) {
      options(rasterFiletype = format)
    }
    else {
      warning(paste("Cannot set filetype to unknown or unsupported file format:", 
                    format, ". See writeFormats()"))
    }
  }
  setOverwrite <- function(overwrite) {
    if (is.logical(overwrite)) {
      options(rasterOverwrite = overwrite)
    }
    else {
      warning(paste("Could not set overwrite. It must be a logical value"))
    }
  }
  setDestFolder <- function(tmpdir) {
    if (!missing(tmpdir)) {
      tmpdir <- trim(tmpdir)
      if (tmpdir != "") {
        lastchar = substr(tmpdir, .nchar(tmpdir), .nchar(tmpdir))
        if (lastchar != "/" & lastchar != "\\") {
          tmpdir <- paste(tmpdir, "/", sep = "")
        }
        options(rasterTmpDir = tmpdir)
      }
    }
  }
  setToDisk <- function(todisk) {
    if (is.logical(todisk)) {
      options(rasterToDisk = todisk)
    }
    else {
      warning(paste("todisk argument must be a logical value"))
    }
  }
  depracatedWarnings <- function(x) {
    if (is.logical(x)) {
      if (is.na(x)) {
        x <- TRUE
      }
      options(rasterDepracatedWarnings = x)
    }
  }
  cnt <- 0
  if (default) {
    cnt <- 1
    options(rasterFiletype = "raster")
    options(rasterOverwrite = FALSE)
    options(rasterDatatype = "FLT8S")
    options(rasterProgress = "none")
    options(rasterTimer = FALSE)
    options(rasterTmpDir = tmpDir(create = FALSE))
    options(rasterTmpTime = 24 * 7)
    options(rasterToDisk = FALSE)
    options(rasterSetFileExt = TRUE)
    options(rasterChunkSize = 1e+06)
    options(rasterMaxMemory = 1e+07)
    options(rasterTolerance = 0.1)
    options(rasterStandardNames = TRUE)
    options(rasterDepracatedWarnings = TRUE)
    options(rasterAddHeader = "")
    v <- utils::packageDescription("raster")[["Version"]]
  }
  if (!missing(format)) {
    setFiletype(format)
    cnt <- cnt + 1
  }
  if (!missing(overwrite)) {
    setOverwrite(overwrite)
    cnt <- cnt + 1
  }
  if (!missing(datatype)) {
    setDataType(datatype)
    cnt <- cnt + 1
  }
  if (!missing(progress)) {
    setProgress(progress)
    cnt <- cnt + 1
  }
  if (!missing(timer)) {
    setTimer(timer)
    cnt <- cnt + 1
  }
  if (!missing(tmpdir)) {
    setTmpdir(tmpdir)
    cnt <- cnt + 1
  }
  if (!missing(tmptime)) {
    setTmpTime(tmptime)
    cnt <- cnt + 1
  }
  if (!missing(todisk)) {
    setToDisk(todisk)
    cnt <- cnt + 1
  }
  if (!missing(setfileext)) {
    setFileExt(setfileext)
    cnt <- cnt + 1
  }
  if (!missing(maxmemory)) {
    setMaxMemorySize(maxmemory)
    cnt <- cnt + 1
  }
  if (!missing(chunksize)) {
    setChunksize(chunksize)
    cnt <- cnt + 1
  }
  if (!missing(tolerance)) {
    setTolerance(tolerance)
    cnt <- cnt + 1
  }
  if (!missing(standardnames)) {
    setStandardNames(standardnames)
    cnt <- cnt + 1
  }
  if (!missing(depracatedwarnings)) {
    depracatedWarnings(depracatedwarnings)
    cnt <- cnt + 1
  }
  if (!missing(addheader)) {
    addHeader(addheader)
    cnt <- cnt + 1
  }
  lst <- list(format = .filetype(), overwrite = .overwrite(), 
              datatype = .datatype(), tmpdir = tmpDir(create = FALSE), 
              tmptime = .tmptime(), progress = .progress(), timer = .timer(), 
              chunksize = .chunksize(), maxmemory = .maxmemory(), 
              todisk = .toDisk(), setfileext = .setfileext(), tolerance = .tolerance(), 
              standardnames = .standardnames(), depwarning = .depracatedwarnings(), 
              addheader = .addHeader())
  save <- FALSE
  if (save) {
    v <- utils::packageDescription("raster")[["Version"]]
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
    oplst <- c(oplst, paste("rasterProgress='", lst$progress, 
                            "'", sep = ""))
    oplst <- c(oplst, paste("rasterTimer=", lst$timer, sep = ""))
    oplst <- c(oplst, paste("rasterChunkSize=", lst$chunksize, 
                            sep = ""))
    oplst <- c(oplst, paste("rasterMaxMemory=", lst$maxmemory, 
                            sep = ""))
    oplst <- c(oplst, paste("rasterSetFileExt=", lst$setfileext, 
                            sep = ""))
    oplst <- c(oplst, paste("rasterTolerance=", lst$tolerance, 
                            sep = ""))
    oplst <- c(oplst, paste("rasterStandardNames=", lst$standardnames, 
                            sep = ""))
    oplst <- c(oplst, paste("rasterDepracatedWarnings=", 
                            lst$depwarning, sep = ""))
    oplst <- c(oplst, paste("rasterAddHeader=", lst$addheader, 
                            sep = ""))
    r <- try(write(unlist(oplst), fn), silent = TRUE)
    cnt <- 1
  }
  if (cnt == 0) {
    cat("format        :", lst$format, "\n")
    cat("datatype      :", lst$datatype, "\n")
    cat("overwrite     :", lst$overwrite, "\n")
    cat("progress      :", lst$progress, "\n")
    cat("timer         :", lst$timer, "\n")
    cat("chunksize     :", lst$chunksize, "\n")
    cat("maxmemory     :", lst$maxmemory, "\n")
    cat("tmpdir        :", lst$tmpdir, "\n")
    cat("tmptime       :", lst$tmptime, "\n")
    cat("setfileext    :", lst$setfileext, "\n")
    cat("tolerance     :", lst$tolerance, "\n")
    cat("standardnames :", lst$standardnames, "\n")
    cat("warn depracat.:", lst$depwarning, "\n")
    if (lst$addheader == "") {
      cat("header        : none\n")
    }
    else {
      cat("header        :", lst$addheader, "\n")
    }
    if (lst$todisk) {
      cat("todisk        : TRUE\n")
    }
  }
  invisible(lst)
}
