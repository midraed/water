# Satellite for sat = "auto": from the MTL file or, for the loaders, from the
# MTL file or the Landsat file names in path. The working directory is never
# searched
getSat <- function(path = NULL, MTL = NULL){
  if(is.null(MTL) && is.null(path)){
    stop("sat = \"auto\" needs the MTL file: pass MTL, or set sat (\"L7\", ",
         "\"L8\", \"L9\")", call. = FALSE)
  }
  # 1. SPACECRAFT_ID from the MTL file, same in every Landsat product format
  if(is.null(MTL)){
    MTL <- list.files(path = path, pattern = "MTL\\.txt$", full.names = T)
  } else {
    MTL <- findMTL(MTL)
  }
  where <- if(is.null(path)) MTL else path
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
    stop("Files from more than one Landsat satellite in ", where, ". Set sat.")
  }
  if(length(num) == 1 && !is.na(num) && num %in% c("7", "8", "9")){
    return(paste0("L", num))
  }
  stop("Can't establish satellite from MTL file or file names in ", where,
       ". Set sat = \"L7\", \"L8\", \"L9\" or \"MODIS\".")
}

# Bands by name (B, G, R, NIR, SWIR1, SWIR2, Thermal1, ...) when the image has
# them, so the layer order does not matter; otherwise the first layers in order
landsatBands <- function(image, bands){
  if(all(bands %in% names(image))) image[[bands]] else image[[seq_along(bands)]]
}

# x must be exactly one of choices: no partial matching (match.arg() takes
# G.method = "T" as "Tasumi"), and a hint for typos and case changes
checkChoice <- function(x, choices, name){
  ok <- is.character(x) && length(x) == 1 && !is.na(x)
  if(ok && x %in% choices){return(x)}
  hint <- NULL
  if(ok){
    d <- utils::adist(tolower(x), tolower(choices))[1, ]
    if(min(d) <= 2){hint <- choices[which.min(d)]}
  }
  stop("Unknown ", name, " = ",
       if(ok) paste0("\"", x, "\"") else paste(deparse(x), collapse = " "),
       ". Use one of: ", paste0("\"", choices, "\"", collapse = ", "),
       if(!is.null(hint)) paste0(" (did you mean \"", hint, "\"?)"), call. = FALSE)
}

# Choices of the method arguments
waterChoices <- list(
  sat = c("auto", "L7", "L8", "L9", "MODIS"),
  alb.coeff = c("Tasumi", "Olmedo", "Liang"),
  LAI.method = c("metric2010", "metric", "MCB", "vineyard", "turner"),
  LST.method = c("SC", "SW"),
  G.method = c("Tasumi", "Bastiaanssen"),
  Zom.method = c("short.crops", "custom", "Perrier"),
  anchors.method = c("flexible", "best", "random")
)

# anchors.method, with the old CITRA names mapped (deprecated)
checkAnchorsMethod <- function(anchors.method){
  old <- c("CITRA-MCB" = "best", "CITRA-MCBbc" = "best", "CITRA-MCBr" = "random")
  if(is.character(anchors.method) && length(anchors.method) == 1 &&
     anchors.method %in% names(old)){
    warning("anchor method names have changed. Old names (CITRA-MCBx) are ",
            "deprecated. Options now include 'best', 'random' and 'flexible'")
    anchors.method <- old[[anchors.method]]
  }
  checkChoice(anchors.method, waterChoices$anchors.method, "anchors.method")
}

# Method arguments of METRIC.EB(), METRIC.Rn() and METRIC.G(), checked before
# any raster work
checkMetricArgs <- function(sat, alb.coeff, LAI.method, LST.method,
                            G.method = "Tasumi", Zom.method = "short.crops"){
  checkChoice(sat, waterChoices$sat, "sat")
  checkChoice(alb.coeff, waterChoices$alb.coeff, "alb.coeff")
  checkChoice(LAI.method, waterChoices$LAI.method, "LAI.method")
  checkChoice(LST.method, waterChoices$LST.method, "LST.method")
  checkChoice(G.method, waterChoices$G.method, "G.method")
  checkChoice(Zom.method, waterChoices$Zom.method, "Zom.method")
  invisible(TRUE)
}

# The split window needs the two thermal bands of Landsat 8 and 9
checkLSTmethod <- function(LST.method, sat, used, name = "LST.method"){
  if(used && LST.method == "SW" && sat == "L7"){
    stop(name, " = \"SW\" (split window) needs two thermal bands: only for ",
         "Landsat 8 and 9", call. = FALSE)
  }
  invisible(TRUE)
}

# extraParameters of METRIC.EB(): a named numeric vector with known names, and
# the ones needed by Zom.method
checkExtraParameters <- function(extraParameters, Zom.method){
  allowed <- c("a", "b", "fLAI", "h", "deltaTemp", "minDist", "WSbuffer")
  if(length(extraParameters) > 0){
    if(!is.numeric(extraParameters) || is.null(names(extraParameters)) ||
       any(names(extraParameters) == "")){
      stop("extraParameters must be a named numeric vector, e.g. ",
           "c(fLAI = 0.5, h = 1.8)", call. = FALSE)
    }
    unknown <- setdiff(names(extraParameters), allowed)
    if(length(unknown)){
      stop("Unknown extraParameters: ", paste(unknown, collapse = ", "),
           ". Allowed: ", paste(allowed, collapse = ", "), call. = FALSE)
    }
  }
  needed <- switch(Zom.method, custom = c("a", "b"), Perrier = c("fLAI", "h"), NULL)
  if(length(needed) && !all(needed %in% names(extraParameters))){
    stop("Zom.method = \"", Zom.method, "\" needs extraParameters ",
         paste(needed, collapse = " and "), call. = FALSE)
  }
  invisible(TRUE)
}

# One element of x at random; sample(x, 1) samples from 1:x when length(x) == 1
sampleOne <- function(x){
  x[sample.int(length(x), 1)]
}

# Path of the Landsat metadata (MTL) file: MTL itself, checked, or the only
# *MTL.txt file in path when MTL is NULL and path is given. The working
# directory is never searched
findMTL <- function(MTL = NULL, path = NULL, required = TRUE){
  if(!is.null(MTL) && length(MTL) == 0){
    if(required){stop("MTL is empty: no Landsat MTL file was given")}
    return(NULL)
  }
  if(is.null(MTL) && is.null(path)){
    if(required){
      stop("MTL is required: pass the path of the Landsat *_MTL.txt file ",
           "(MTL = ...)", call. = FALSE)
    }
    return(NULL)
  }
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

# Numeric values from a Landsat MTL file, e.g. mtlValues(MTL, "K1_CONSTANT_BAND_10").
# Level-2 MTL files repeat some keys in several groups (e.g. REFLECTANCE_MULT_BAND_4
# for the SR scaling and for TOA): such keys must be read from one group
mtlValues <- function(MTL, keys, group = NULL){
  lines <- readLines(MTL[1], warn = FALSE)
  if(!is.null(group)){
    start <- grep(paste0("^\\s*GROUP = ", group, "\\s*$"), lines)
    end <- grep(paste0("^\\s*END_GROUP = ", group, "\\s*$"), lines)
    if(length(start) != 1 || length(end) != 1){
      stop("Group ", group, " not found in MTL file")
    }
    lines <- lines[start:end]
  }
  values <- vapply(keys, function(key){
    line <- grep(paste0("^\\s*", key, " = "), lines, value = TRUE)
    if(length(line) > 1){
      stop(key, " appears more than once in the MTL file; give its group")
    }
    as.numeric(gsub("\"", "", sub(".*= *", "", line[1])))
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
  terra::tmpFiles(current = FALSE, orphan = TRUE, remove = TRUE)
  return(stack)
}


# Crop to aoi only when it is passed (an object called aoi in the workspace is
# never used)
aoiCrop <- function(raster, aoi){
  if(missing(aoi)){return(raster)}
  cropToAoi(raster, aoi)
}

# Crop x to aoi, reprojecting a polygon aoi to the CRS of x when both have one
cropToAoi <- function(x, aoi){
  if(inherits(aoi, c("sf", "sfc"))){aoi <- terra::vect(aoi)}
  if(inherits(aoi, "SpatVector") && terra::crs(aoi) != "" && terra::crs(x) != "" &&
     !terra::same.crs(aoi, x)){
    aoi <- terra::project(aoi, x)
  }
  crop(x, aoi)
}

# TRUE if R knows the time zone (unknown ones only give a warning and use UTC)
knownTimezone <- function(tz){
  tryCatch({as.POSIXlt(Sys.time(), tz = tz); TRUE}, warning = function(w) FALSE)
}

isURL <- function(x){grepl("^(https?|ftp)://", x)}

# Longitude (decimal degrees, negative west) of the centre of the time zone of
# the times x: 15 degrees per hour of UTC offset, including daylight saving time
tzLongitude <- function(x){
  15 * as.POSIXlt(x)$gmtoff / 3600
}

# 1 x 1 degree DEM tiles (latitude and longitude of the south-west corner)
# covering the extent of an image
demTiles <- function(image, margin = 0){
  area <- ext(terra::project(terra::as.polygons(ext(image), crs = terra::crs(image)),
                             "EPSG:4326"))
  expand.grid(lat = floor(area[3] - margin):floor(area[4] + margin),
              lon = floor(area[1] - margin):floor(area[2] + margin))
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
