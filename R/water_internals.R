getSat <- function(path = getwd(), MTL = NULL){
  # 1. SPACECRAFT_ID from the MTL file, same in every Landsat product format
  if(is.null(MTL)){
    MTL <- list.files(path = path, pattern = "MTL\\.txt$", full.names = T)
  } else {
    MTL <- findMTL(MTL)
  }
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

# Bands by name (B, G, R, NIR, SWIR1, SWIR2, Thermal1, ...) when the image has
# them, so the layer order does not matter; otherwise the first layers in order
landsatBands <- function(image, bands){
  if(all(bands %in% names(image))) image[[bands]] else image[[seq_along(bands)]]
}

# One element of x at random; sample(x, 1) samples from 1:x when length(x) == 1
sampleOne <- function(x){
  x[sample.int(length(x), 1)]
}

# Path of the Landsat metadata (MTL) file: MTL itself, checked, or the only
# *MTL.txt file in path when MTL is NULL
findMTL <- function(MTL = NULL, path = getwd(), required = TRUE){
  if(is.null(MTL)){
    MTL <- list.files(path = path, pattern = "MTL\\.txt$", full.names = TRUE)
    if(length(MTL) > 1){
      stop("More than one MTL file in ", path, " (",
           paste(basename(MTL), collapse = ", "), "). Pass the one to use with MTL")
    }
    if(length(MTL) == 0){
      if(required){stop("No Landsat MTL file found in ", path, ". Pass its path with MTL")}
      return(NULL)
    }
  }
  if(!is.character(MTL) || length(MTL) != 1){
    stop("MTL must be the path to the Landsat *_MTL.txt file (a character ",
         "string); pass the file name, not its contents")
  }
  if(!file.exists(MTL)){stop("MTL file not found: ", MTL)}
  MTL
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
  terra::tmpFiles(current = FALSE, old = TRUE, remove = TRUE)
  return(stack)
}


aoiCrop <- function(raster, aoi){
  if(!missing(aoi)){
    raster <- crop(raster,aoi)
    return(raster)
  }
  if(missing(aoi) && isTRUE(getOption("waterAutoAoi")) &&
     exists(x = "aoi", envir = .GlobalEnv, inherits = FALSE)){
    aoi <- get(x = "aoi", envir=.GlobalEnv)
    raster <- crop(raster,aoi)
    return(raster)
  }
  return(raster)
}

# Hint for failed anchor selections
anchorHint <- function(albedo){
  alb <- stats::median(values(albedo), na.rm = TRUE)
  paste0("Median albedo is ", round(alb, 2),
         if(isTRUE(alb > 0.5)) " (too high: check the surface reflectance scaling, e.g. use loadImageSR())",
         ". Mask clouds, use an area with both well-watered full-cover vegetation and ",
         "dry bare soil, or try anchors.method = \"flexible\", a smaller n, a larger ",
         "WSbuffer, or provide anchors.")
}
