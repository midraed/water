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
    MTL <- list.files(path = path, pattern = "MTL\\.txt$", full.names = TRUE)
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

# Bands by name (B, G, R, NIR, SWIR1, SWIR2, Thermal1, ...), so the layer
# order does not matter. Layers without any standard name (e.g. a stack made
# by hand) are taken by position pos, with a warning; a stack with standard
# names but without one of the bands is an error
landsatBands <- function(image, bands, pos = seq_along(bands)){
  if(all(bands %in% names(image))){return(image[[bands]])}
  standard <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2", "SWIR3", "Thermal1", "Thermal2")
  if(!any(standard %in% names(image))){
    warning("image layers are not named ", paste(bands, collapse = ", "),
            ": using layers ", paste(pos, collapse = ", "), " by position (name ",
            "them, e.g. names(x) <- c(\"B\", \"G\", \"R\", \"NIR\", \"SWIR1\", ",
            "\"SWIR2\"))", call. = FALSE)
    return(image[[pos]])
  }
  stop("image has no layer named ",
       paste0("\"", setdiff(bands, names(image)), "\"", collapse = ", "),
       " (layers: ", paste(names(image), collapse = ", "), ")", call. = FALSE)
}

# Cells of the hot and cold anchors of calcH(), checked: columns and labels
# ("Hot" and "Cold" are accepted), inside the image, no NA inputs, and hot
# anchors at least min.dT warmer than the cold ones (for method = "lm", in
# every pair); the area must contain both wet and dry extremes
anchorCells <- function(anchors, Ts, Ts.datum, Rn, G, Z.om, DEM, method = "mean",
                        min.dT = 5){
  if(inherits(anchors, "sf")){
    if(!"type" %in% names(anchors)){
      stop("anchors must have a type column (\"hot\" or \"cold\")", call. = FALSE)
    }
    if(!is.na(sf::st_crs(anchors)) && terra::crs(Ts) != ""){
      anchors <- sf::st_transform(anchors, terra::crs(Ts))
    }
    xy <- sf::st_coordinates(anchors)
  } else {
    if(!is.data.frame(anchors) || !all(c("X", "Y", "type") %in% names(anchors))){
      stop("anchors must have columns X, Y and type (\"hot\" or \"cold\")", call. = FALSE)
    }
    xy <- cbind(anchors$X, anchors$Y)
  }
  type <- tolower(trimws(as.character(anchors$type)))
  wrong <- setdiff(unique(type), c("hot", "cold"))
  if(length(wrong)){
    stop("anchors type must be \"hot\" or \"cold\" (found ",
         paste0("\"", wrong, "\"", collapse = ", "), ")", call. = FALSE)
  }
  if(!all(c("hot", "cold") %in% type)){
    stop("anchors need at least one hot and one cold pixel", call. = FALSE)
  }
  cells <- terra::cellFromXY(Ts, xy)
  if(anyNA(cells)){
    stop("anchor(s) ", paste(which(is.na(cells)), collapse = ", "), " are outside ",
         "the image", call. = FALSE)
  }
  at <- function(r){if(inherits(r, "SpatRaster")) as.numeric(r[cells][[1]]) else rep(r, length(cells))}
  missing.input <- Reduce(`|`, lapply(list(Ts, Rn, G, Z.om, DEM), function(r) is.na(at(r))))
  if(any(missing.input)){
    stop("Ts, Rn, G, Z.om or DEM is NA at anchor(s) ",
         paste(which(missing.input), collapse = ", "), " (masked pixel?): choose ",
         "other anchors", call. = FALSE)
  }
  hot <- cells[type == "hot"]
  cold <- cells[type == "cold"]
  Td <- at(Ts.datum)
  dT <- if(method == "lm"){
    n <- min(length(hot), length(cold))
    Td[type == "hot"][1:n] - Td[type == "cold"][1:n]
  } else mean(Td[type == "hot"]) - mean(Td[type == "cold"])
  which.pair <- function(i){if(method == "lm") paste0(" (pair ", i, ")") else ""}
  i <- which.min(dT)
  if(dT[i] <= 0){
    stop("The hot anchors are not warmer than the cold anchors", which.pair(i),
         ": difference ", round(dT[i], 2), " K. Are \"hot\" and \"cold\" swapped?",
         call. = FALSE)
  }
  if(dT[i] < min.dT){
    stop("The hot and cold anchors differ by only ", round(dT[i], 2), " K",
         which.pair(i), " (minimum ", min.dT, " K): the area must contain both ",
         "well-watered full-cover vegetation and dry bare soil. Use a larger ",
         "area (tens of km) and crop the results, or choose the anchors",
         call. = FALSE)
  }
  list(hot = hot, cold = cold, dT = dT)
}

# Stop when the H iteration of calcH() diverges: a and b of dT = a Ts + b,
# or the friction velocity u at the anchors, not finite (u also not positive)
checkIteration <- function(a = 0, b = 0, u = 1){
  if(!all(is.finite(c(a, b, u))) || any(u <= 0)){
    stop("The H iteration diverged: check the weather data at the overpass ",
         "(temp, RH, wind), the roughness (Z.om) and the anchors", call. = FALSE)
  }
  invisible(TRUE)
}

# Message when the scene has clouds (CLOUD_COVER in the MTL file) and none
# of the images was masked with cfmask(), which tags its output
checkCloudMask <- function(MTL, ...){
  images <- Filter(Negate(is.null), list(...))
  cover <- tryCatch(suppressWarnings(mtlValues(MTL, "CLOUD_COVER")),
                    error = function(e) NA)
  if(!length(images) || !isTRUE(cover > 0)){return(invisible(FALSE))}
  tagged <- vapply(images, function(x){
    m <- terra::metags(x)
    "water_cfmask" %in% (if(is.data.frame(m)) m$name else names(m))
  }, logical(1))
  if(!any(tagged)){
    message("The scene has ", cover, "% clouds (CLOUD_COVER in the MTL file) and ",
            "the image was not masked with cfmask() (or its mark was lost, e.g. ",
            "by selecting layers afterwards): clouds and their shadows can be ",
            "taken as cold anchor pixels")
  }
  invisible(!any(tagged))
}

# Warning when the station recorded more than limit mm of rain in the days
# before the overpass (if its records cover them): the hot pixels may not be
# dry, and calcH() assumes ET = 0 there unless ETrF.hot is set
checkRecentRain <- function(WeatherStation, days = 3, limit = 5){
  if(!inherits(WeatherStation, "waterWeatherStation") || is.null(WeatherStation$at.sat)){
    return(invisible(NA))
  }
  rain <- WeatherStation$alldata$rain
  t <- as.POSIXct(WeatherStation$alldata$datetime)
  overpass <- as.POSIXct(WeatherStation$at.sat$datetime)
  if(all(is.na(rain)) || min(t) > overpass - days * 86400){return(invisible(NA))}
  total <- sum(rain[t > overpass - days * 86400 & t <= overpass], na.rm = TRUE)
  if(total > limit){
    warn("The station recorded ", round(total, 1), " mm of rain in the ", days,
         " days before the image: the hot anchor pixels may not be dry (ET = 0 is ",
         "assumed). Consider calcH.control(ETrF.hot = ) (Allen et al., 2007, Eq. 47)")
  }
  invisible(total)
}

# Albedo coefficients for alb.coeff = "auto": Olmedo for Landsat 8 and 9,
# Tasumi et al. (2008) for Landsat 7 and MODIS
albedoCoeff <- function(coeff, sat){
  if(coeff != "auto"){return(coeff)}
  switch(sat, L7 = "Tasumi", L8 = , L9 = "Olmedo", MODIS = "Tasumi", "Olmedo")
}

# The inputs of the METRIC functions must be on the grid of the first one
# (extent, resolution and CRS)
checkGrids <- function(...){
  r <- Filter(Negate(is.null), list(...))
  if(length(r) < 2){return(invisible(TRUE))}
  for(n in names(r)[-1]){
    if(!terra::compareGeom(r[[1]], r[[n]], stopOnError = FALSE)){
      stop(n, " is not on the grid of ", names(r)[1], " (extent, resolution or ",
           "CRS differ, within terra's tolerance): use e.g. terra::project(", n,
           ", ", names(r)[1], "[[1]])", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# aoi of the functions that do not use it
unusedAoi <- function(missing.aoi){
  if(!missing.aoi){
    warn("aoi is not used by this function: crop the inputs beforehand (e.g. ",
         "loadImageSR(aoi = ))")
  }
}

# DEM tiles from their URLs, cropped to area. A tile that does not exist
# (HTTP 404: the ocean) is skipped; any other error (e.g. no connection) stops
readDEMtiles <- function(url, area){
  parts <- list()
  missing <- integer()
  for(i in seq_along(url)){
    tile <- tryCatch(crop(rast(url[i]), area), error = function(e) e)
    if(inherits(tile, "error")){
      if(grepl("404|does not exist|No such file", conditionMessage(tile))){
        missing <- c(missing, i)
      } else {
        stop("The DEM tile ", url[i], " could not be read (internet connection?): ",
             conditionMessage(tile), call. = FALSE)
      }
    } else {
      parts[[length(parts) + 1]] <- tile
    }
  }
  list(parts = parts, missing = missing)
}

# Surface temperature must be in Kelvin. Below 200 K is Celsius (all cells),
# or unmasked fill or cloud tops (some cells: a warning, unless they are many)
checkKelvin <- function(Ts, name = "Ts"){
  r <- as.numeric(global(Ts[[1]], "range", na.rm = TRUE))
  if(!all(is.finite(r)) || r[1] >= 200){return(invisible(TRUE))}
  if(r[2] < 200){
    stop(name, " must be in Kelvin: the values look like Celsius (maximum ",
         round(r[2], 1), "); add 273.15", call. = FALSE)
  }
  n <- global(Ts[[1]] < 200, "sum", na.rm = TRUE)[1, 1]
  valid <- global(!is.na(Ts[[1]]), "sum")[1, 1]
  text <- paste0(name, " has ", n, " cells (", signif(100 * n / valid, 2),
                 "%) below 200 K (minimum ", round(r[1], 1), " K): unmasked fill ",
                 "or clouds, or wrong scaling")
  if(n > 0.01 * valid){
    stop(text, ". Mask them (cfmask(); loadImage() and loadImageSR() mask the ",
         "fill)", call. = FALSE)
  }
  warning(text, ". Mask clouds with cfmask()", call. = FALSE)
  invisible(FALSE)
}

# Run expr without the "below 200 K" warning of checkKelvin(), already given
# once (METRIC.EB() checks Ts before calcAnchors() and calcH())
muffleKelvin <- function(expr){
  withCallingHandlers(expr, warning = function(w){
    if(grepl("below 200 K", conditionMessage(w))){invokeRestart("muffleWarning")}
  })
}

# x must be exactly one of choices: no partial matching (match.arg() takes
# G.method = "T" as "Tasumi"), and a hint for typos and case changes
checkChoice <- function(x, choices, name){
  ok <- is.character(x) && length(x) == 1 && !is.na(x)
  if(ok && x %in% choices){return(x)}
  hint <- NULL
  if(ok){
    d <- utils::adist(tolower(x), tolower(choices))[1, ]
    # up to 2 edits, fewer for short names ("L8" is not a typo of "L2")
    if(min(d) <= min(2, nchar(x) %/% 3)){hint <- choices[which.min(d)]}
  }
  stop("Unknown ", name, " = ",
       if(ok) paste0("\"", x, "\"") else paste(deparse(x), collapse = " "),
       ". Use one of: ", paste0("\"", choices, "\"", collapse = ", "),
       if(!is.null(hint)) paste0(" (did you mean \"", hint, "\"?)"), call. = FALSE)
}

# Choices of the method arguments
waterChoices <- list(
  sat = c("auto", "L7", "L8", "L9", "MODIS"),
  alb.coeff = c("auto", "Tasumi", "Olmedo", "Liang"),
  LAI.method = c("metric2010", "metric", "MCB", "vineyard", "turner"),
  LST.method = c("auto", "L2", "METRIC", "SW"),
  G.method = c("Tasumi", "Bastiaanssen"),
  Zom.method = c("short.crops", "custom", "Perrier"),
  anchors.method = c("percentile", "flexible", "best", "random")
)

# anchors.method, with the old CITRA names mapped (deprecated)
checkAnchorsMethod <- function(anchors.method){
  old <- c("CITRA-MCB" = "best", "CITRA-MCBbc" = "best", "CITRA-MCBr" = "random")
  if(is.character(anchors.method) && length(anchors.method) == 1 &&
     anchors.method %in% names(old)){
    warn("anchor method names have changed. Old names (CITRA-MCBx) are ",
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
  removedSC(LST.method, "LST.method")
  checkChoice(LST.method, waterChoices$LST.method, "LST.method")
  checkChoice(G.method, waterChoices$G.method, "G.method")
  checkChoice(Zom.method, waterChoices$Zom.method, "Zom.method")
  invisible(TRUE)
}

# The split window needs the two thermal bands of Landsat 8 and 9; the
# Level-2 surface temperature is read only for Landsat 8 and 9
checkLSTmethod <- function(LST.method, sat, used, name = "LST.method"){
  if(used && LST.method == "SW" && sat == "L7"){
    stop(name, " = \"SW\" (split window) needs two thermal bands: only for ",
         "Landsat 8 and 9", call. = FALSE)
  }
  if(used && LST.method == "L2" && sat == "L7"){
    stop(name, " = \"L2\" (Level-2 surface temperature) is only for Landsat 8 ",
         "and 9; use \"METRIC\", or pass surf.temp", call. = FALSE)
  }
  invisible(TRUE)
}

# "SC" (single channel) was removed in version 1.0
removedSC <- function(method, name){
  if(identical(method, "SC")){
    stop(name, " = \"SC\" was removed in version 1.0. Use \"auto\" (the ",
         "Level-2 surface temperature for Landsat 8 and 9, the METRIC equations ",
         "otherwise), \"L2\" or \"METRIC\"", call. = FALSE)
  }
}

# The surface temperature to use: "surf.temp" when it is supplied; "auto" is
# "L2" for Landsat 8 and 9 with Band_10 in image.SR, "METRIC" otherwise
resolveLST <- function(LST.method, sat, image.SR = NULL, has.surf.temp = FALSE){
  if(has.surf.temp){return("surf.temp")}
  hasB10 <- !is.null(image.SR) && "Band_10" %in% names(image.SR)
  if(LST.method == "auto"){
    LST.method <- if(sat %in% c("L8", "L9") && hasB10) "L2" else "METRIC"
  }
  if(LST.method == "L2" && !hasB10){
    stop("LST.method = \"L2\" needs the Level-2 surface temperature: ",
         "image.SR$Band_10 from loadImageSR() (with the *_ST_B10.TIF file), or ",
         "pass surf.temp", call. = FALSE)
  }
  LST.method
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

# Named columns of read.WSdata() as the positional vector c(date, time,
# radiation, wind, RH, temp, rain); datetime sets date and time, rain is
# optional and pp is an old name of rain
columnsByName <- function(columns){
  roles <- c("date", "time", "radiation", "wind", "RH", "temp", "rain")
  # what the vector meant when only its order was used
  old <- unname(c(columns, rep(NA, 7))[1:7])
  columns.in <- columns
  nm <- names(columns)
  nm[nm == "pp"] <- "rain"
  unknown <- nm[is.na(nm) | !nm %in% c(roles, "datetime")]
  if(length(unknown)){
    stop("Unknown names in columns: ", paste0("\"", unknown, "\"", collapse = ", "),
         ". Use date, time (or datetime), radiation, wind, RH, temp and rain",
         call. = FALSE)
  }
  if(anyDuplicated(nm)){
    stop("Names repeated in columns: ", paste(unique(nm[duplicated(nm)]), collapse = ", "),
         call. = FALSE)
  }
  if("datetime" %in% nm){
    if(any(c("date", "time") %in% nm)){
      stop("columns: use datetime, or date and time, not both", call. = FALSE)
    }
    dt <- columns[nm == "datetime"]
    columns <- c(dt, dt, columns[nm != "datetime"])
    nm <- c("date", "time", nm[nm != "datetime"])
  }
  names(columns) <- nm
  needed <- setdiff(roles[1:6], nm)
  if(length(needed)){
    stop("columns has no ", paste0("\"", needed, "\"", collapse = ", "), call. = FALSE)
  }
  positional <- unname(columns[roles])  # NA for a missing rain
  if(!"datetime" %in% names(columns.in) && !identical(as.numeric(old), as.numeric(positional))){
    warn("columns are matched by name (before version 1.0 only their ",
            "order was used): check that the names are right")
  }
  positional
}

# Plausible weather station values (after cf): impossible values stop, unusual
# ones warn with the likely unit problem. hours: time span of the records
checkWSvalues <- function(radiation, wind, RH, temp, rain, hours){
  rng <- function(x){
    if(all(is.na(x))) c(NA, NA) else range(x, na.rm = TRUE)
  }
  span <- function(r, units){paste0("from ", signif(r[1], 4), " to ", signif(r[2], 4), " ", units)}
  t <- rng(temp)
  if(!anyNA(t)){
    if(t[1] > 200){
      stop("Temperature ", span(t, "C"), ": in Kelvin? Use Celsius ",
           "(subtract 273.15)", call. = FALSE)
    }
    if(t[1] < -60 || t[2] > 60){
      stop("Temperature ", span(t, "C"), " is outside -60 to 60 C: in Fahrenheit? ",
           "Use Celsius, (F - 32) / 1.8 (see cf)", call. = FALSE)
    }
    if(t[2] > 45){
      warn("Temperature up to ", signif(t[2], 4), " C: check that it is in Celsius")
    }
  }
  h <- rng(RH)
  if(!anyNA(h)){
    if(h[1] < 0 || h[2] > 110){
      stop("RH ", span(h, "%"), " is outside 0 to 110%", call. = FALSE)
    }
    if(h[2] <= 1){
      stop("RH ", span(h, "%"), ": as a fraction? Use percent (multiply by 100)",
           call. = FALSE)
    }
  }
  w <- rng(wind)
  if(!anyNA(w)){
    if(w[1] < 0 || w[2] > 75){
      stop("Wind speed ", span(w, "m/s"), " is outside 0 to 75 m/s", call. = FALSE)
    }
    if(w[2] > 20){
      warn("Wind speed up to ", signif(w[2], 4), " m/s: in km/h? Use ",
              "cf = c(1, 0.2777778, 1)")
    }
  }
  s <- rng(radiation)
  if(!anyNA(s)){
    if(s[1] < -50 || s[2] > 2000){
      stop("Solar radiation ", span(s, "W/m2"), " is outside -50 to 2000 W/m2: ",
           "check the units (W/m2) and cf", call. = FALSE)
    }
    if(s[2] > 1400){
      warn("Solar radiation up to ", signif(s[2], 4), " W/m2, above the ",
              "usual clear-sky maximum: check the units (W/m2) and the sensor")
    }
    if(s[2] < 50 && isTRUE(hours >= 12)){
      warn("Solar radiation is at most ", signif(s[2], 4), " W/m2 over ",
              round(hours), " hours: in MJ m-2 h-1? Use W/m2 (multiply by ",
              "277.8, e.g. with cf)")
    }
  }
  r <- rng(rain)
  if(!anyNA(r) && r[1] < 0){
    stop("Rain from ", signif(r[1], 4), " mm: negative values", call. = FALSE)
  }
  invisible(TRUE)
}

# Station radiation at the overpass (Rs) against clear sky (Rso), both W/m2.
# METRIC needs a clear image, so Rs/Rso should be near 1 (ASCE-EWRI 2005 uses
# about +-10% on clear days); 0.8 leaves room for haze. Not checked at low sun
checkOverpassRadiation <- function(Rs, Rso, overpass.radiation = "measured"){
  if(!isTRUE(Rso >= 200) || is.na(Rs)){return(invisible(NA))}
  ratio <- Rs / Rso
  if(ratio < 0.8){
    warn("Radiation at the overpass is ", round(Rs), " W/m2, ", round(100 * ratio),
            "% of the clear-sky value (", round(Rso), " W/m2): clouds near the ",
            "station, or a faulty or miscalibrated pyranometer. ",
            if(overpass.radiation == "measured") {
              "Check the sensor, or use overpass.radiation = \"clear-sky\""
            } else "The clear-sky value is used (overpass.radiation = \"clear-sky\")")
  }
  if(ratio > 1.1){
    warn("Radiation at the overpass is ", round(Rs), " W/m2, ", round(100 * ratio),
            "% of the clear-sky value (", round(Rso), " W/m2): check the units ",
            "and the calibration of the pyranometer")
  }
  invisible(ratio)
}

# Kt, the clearness (turbidity) coefficient of the atmosphere: 0 < Kt <= 1
checkKt <- function(Kt){
  if(!is.numeric(Kt) || length(Kt) != 1 || !isTRUE(Kt > 0 && Kt <= 1)){
    stop("Kt must be a number above 0 and up to 1 (1 for clean air, 0.5 for ",
         "extremely turbid air)", call. = FALSE)
  }
  invisible(Kt)
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


# Results of the package functions. With waterOptions(writeResults = TRUE)
# they are written to disk (for images too large for memory): intermediate
# results to terra's temporary files, with unique names, and the final outputs
# (finalOutputs) to waterOptions(outputFolder). A final file written earlier
# in this session is not overwritten, as an object may still read it: a new
# name is used, with a warning
finalOutputs <- c("EB", "ET24", "DEM")
written <- new.env()
saveLoadClean <- function(imagestack, stack.names=NULL, file, ...){
  # a single layer without a name is named as its file (e.g. "Rn", "G", "H")
  if(is.null(stack.names) && !missing(file) && nlyr(imagestack) == 1){
    stack.names <- file
  }
  if(!getOption("waterWriteResults")){
    names(imagestack) <- stack.names
    return(imagestack)
  }
  if(missing(file)){file <- deparse(substitute(imagestack))}
  if(file %in% finalOutputs){
    tmpdir <- getOption("waterOutputFolder")
    lastchar = substr(tmpdir, nchar(tmpdir), nchar(tmpdir))
    if (lastchar != "/" & lastchar != "\\") {
      tmpdir <- paste0(tmpdir, "/")
    }
    tstamp <- NULL
    if(!getOption("waterOverwrite")){
      tstamp <- paste0("_", strftime(Sys.time(), format = "%Y%m%d%H%M%S"))}
    resultfile <- paste0(tmpdir, file, tstamp, ".tif")
    if(file.exists(resultfile) &&
       (!getOption("waterOverwrite") || exists(normalizePath(resultfile), envir = written))){
      new <- tempfile(pattern = paste0(file, tstamp, "_"), tmpdir = sub("/$", "", tmpdir),
                      fileext = ".tif")
      warning(resultfile, " exists (written earlier in this session, or a time ",
              "stamp in the same second): saved as ", new, call. = FALSE)
      resultfile <- new
    }
    writeRaster(imagestack, filename = resultfile, overwrite = TRUE)
    assign(normalizePath(resultfile), TRUE, envir = written)
    message("Result saved as ", resultfile)
  } else {
    resultfile <- tempfile(pattern = paste0(file, "_"), tmpdir = terra::terraOptions(print = FALSE)$tempdir,
                           fileext = ".tif")
    writeRaster(imagestack, filename = resultfile)
  }
  stack <- rast(resultfile)
  names(stack) <- stack.names
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
  # on some systems R accepts unknown names silently (as UTC)
  if(length(OlsonNames()) && !tz %in% OlsonNames()){return(FALSE)}
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

# warning() without the call of the internal function that gives it
warn <- function(...){
  warning(..., call. = FALSE)
}

# Acquisition date and time (UTC, POSIXct) from the lines of an MTL file;
# the fractional seconds of SCENE_CENTER_TIME are dropped
mtlDateTime <- function(MTL){
  date.line <- grep("DATE_ACQUIRED", MTL, value = TRUE)
  time.line <- grep("SCENE_CENTER_TIME", MTL, value = TRUE)
  sat.date <- regmatches(date.line, regexec(text = date.line,
                         pattern = "([0-9]{4})(-)([0-9]{2})(-)([0-9]{2})"))[[1]][1]
  sat.time <- regmatches(time.line, regexec(text = time.line,
                         pattern = "([0-9]{2}:[0-9]{2}:[0-9]{2})"))[[1]][1]
  out <- as.POSIXct(paste(sat.date, sat.time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  if(length(out) != 1 || is.na(out)){
    stop("DATE_ACQUIRED and SCENE_CENTER_TIME not found in the MTL file", call. = FALSE)
  }
  out
}
