getSat <- function(path = getwd(), MTL = NULL){
  # 1. SPACECRAFT_ID from the MTL file, same in every Landsat product format
  if(is.null(MTL)){MTL <- list.files(path = path, pattern = "MTL.txt$", full.names = T)}
  MTL <- MTL[file.exists(MTL)]
  if(length(MTL)){
    num <- unique(vapply(MTL, function(f){
      id <- grep("SPACECRAFT_ID", readLines(f, warn = FALSE), value = TRUE)
      sub("\\D*(\\d+)\\D*", "\\1", id[1])
    }, character(1), USE.NAMES = FALSE))
  } else {
    # 2. Landsat file names: LC8..., LE7... (pre-collection) or LC08_, LE07_
    files <- list.files(path = path, pattern = "^L[COTE]0?[789].*\\.(TIF|tif)$")
    num <- unique(sub("^L[COTE]0?([789]).*", "\\1", files))
  }
  if(length(num) > 1){
    stop("Files from more than one Landsat satellite in ", path, ". Set sat.")
  }
  if(length(num) == 1 && !is.na(num) && num %in% c("7", "8", "9")){
    return(paste0("L", num))
  }
  stop("Can't establish satellite from MTL file or file names in ", path,
       ". Set sat = \"L7\", \"L8\", \"L9\" or \"MODIS\".")
}

# Numeric values from a Landsat MTL file, e.g. mtlValues(MTL, "K1_CONSTANT_BAND_10")
mtlValues <- function(MTL, keys){
  lines <- readLines(MTL[1], warn = FALSE)
  values <- vapply(keys, function(key){
    line <- grep(paste0("^\\s*", key, " = "), lines, value = TRUE)[1]
    as.numeric(gsub("\"", "", sub(".*= *", "", line)))
  }, numeric(1))
  if(anyNA(values)){
    stop("Not found in MTL file: ", paste(keys[is.na(values)], collapse = ", "))
  }
  values
}


saveLoadClean <- function(imagestack, stack.names=NULL, file, ...){
  if(!getOption("waterWriteResults")){
    names(imagestack) <- stack.names
    return(imagestack)
  }
  tstamp <- NULL
  if(!getOption("waterOverwrite")){
    tstamp <- paste0("_", strftime(Sys.time(), format = "%Y%m%d%H%M%S"))}
  tmpdir <- getOption("waterOutputFolder")
  lastchar = substr(tmpdir, nchar(tmpdir), nchar(tmpdir))
  if (lastchar != "/" & lastchar != "\\") {
    tmpdir <- paste0(tmpdir, "/")
  }
  if(missing(file)){file <- paste0(deparse(substitute(imagestack)),".tif")}
  resultfile <- paste0(tmpdir, file, tstamp, ".tif")
  writeRaster(imagestack, filename = resultfile, ...)
  message(paste("Result saved as", paste0(tmpdir, file, tstamp, ".tif"), "    OK"))
  stack <- rast(resultfile)
  names(stack) <- stack.names
  terra::tmpFiles(old = TRUE, remove = TRUE)
  return(stack)
}


aoiCrop <- function(raster, aoi){
  if(!missing(aoi)){
    raster <- crop(raster,aoi)
    return(raster)
  }
  if(missing(aoi) & exists(x = "aoi", envir=.GlobalEnv)){
    aoi <- get(x = "aoi", envir=.GlobalEnv)
    raster <- crop(raster,aoi)
    return(raster)
  }
  return(raster)
}


