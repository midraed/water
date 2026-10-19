#' Plot method for waterLSEB S3 class
#' @param x       waterLSEB object.
#' @param ...     further arguments passed to or from other methods.
#' @return No return value, called for side effects (plots the NetRadiation,
#' SoilHeat, SensibleHeat and LatentHeat layers of x$EB).
#' @author Guillermo Federico Olmedo
#' @export
#' @family LSEB objects related functions
#' @method plot waterLSEB
plot.waterLSEB <- function(x, ...){
  plot(x$EB[[1:4]])
}


#' Print method for waterLSEB S3 class
#' @param x        waterLSEB object.
#' @param ...      further arguments passed to or from other methods.
#' @return x, invisibly (called for its side effect of printing the energy
#' balance layers, the weather station, the anchor pixels and the methods).
#' @author Guillermo Federico Olmedo
#' @export
#' @family LSEB objects related functions
#' @method print waterLSEB
print.waterLSEB <- function(x, ...){
  print(x$EB, ...)
  print(x$WeatherStation, ...)
  print(x$anchors, ...)
  print(x$methods, ...)
  if(!is.null(x$calibration)){
    cat("\nCalibration:\n")
    summary.waterLSEB(x)
  }
  invisible(x)
}

#' Summary of the calibration of a METRIC energy balance
#' @description Prints the quality checks of the calibration of a METRIC.EB()
#' result: the anchor pixels, the iteration of H and the share of pixels with
#' unlikely values of the reference ET fraction (ETrF).
#' @param object   waterLSEB object, from METRIC.EB().
#' @param ...      not used.
#' @details In METRIC, the cold anchors have ETrF = ETp.coef (1.05) and few
#' pixels should be above about 1.1 (Allen et al., 2007): many pixels above
#' 1.1 mean that the cold anchors are too warm. Pixels with H > Rn - G have
#' negative ETrF, and LE is set to 0 there.
#' @return The calibration list (object$calibration), invisibly.
#' @author Guillermo Federico Olmedo
#' @export
#' @family LSEB objects related functions
#' @method summary waterLSEB
summary.waterLSEB <- function(object, ...){
  k <- object$calibration
  pct <- function(x){paste0(round(100 * x, 1), "%")}
  cat("Anchors:", sum(object$anchors$type == "hot"), "hot and",
      sum(object$anchors$type == "cold"), "cold")
  if(!is.null(k$candidates)){
    cat(" (from", k$candidates[["hot"]], "hot and", k$candidates[["cold"]],
        "cold candidates)")
  }
  cat("\nTs of the anchors: hot", round(k$Ts.hot, 2), "K, cold", round(k$Ts.cold, 2),
      "K, difference", round(k$dT.anchors, 2), "K\n")
  cat("dT = ", signif(k$a, 4), " * Ts_datum + ", signif(k$b, 4), "; u200 ",
      round(k$u200, 2), " m/s; ", if(isTRUE(all(k$converged))) "converged" else
        "NOT converged", " after ", paste(k$iterations, collapse = ", "),
      " iterations\n", sep = "")
  cat("ETrF at the overpass: 5%", round(k$ETrF[1], 2), " median", round(k$ETrF[2], 2),
      " 95%", round(k$ETrF[3], 2), " max", round(k$ETrF[4], 2), "\n")
  cat("Pixels with ETrF > 1.1:", pct(k$ETrF.above.1.1), " ETrF < 0:", pct(k$ETrF.below.0),
      " H > Rn - G:", pct(k$H.above.Rn.G), "\n")
  invisible(k)
}

#' Plot the anchor pixels of a METRIC energy balance
#' @description Plots the hot (red) and cold (blue) anchor pixels of a
#' METRIC.EB() result, and the weather station (black), over a false-colour
#' composite of the image (NIR, R, G) or over the surface temperature.
#' @param x         waterLSEB object, from METRIC.EB().
#' @param image.SR  Optional. Reflectance image with layers NIR, R and G, e.g.
#' from loadImageSR(). If missing, the surface temperature is plotted.
#' @param ...       further arguments passed to terra::plotRGB() or plot().
#' @return No return value, called for its side effect (a plot).
#' @author Guillermo Federico Olmedo
#' @export
#' @family LSEB objects related functions
plotAnchors <- function(x, image.SR, ...){
  if(!inherits(x, "waterLSEB")){stop("x must be a METRIC.EB() result", call. = FALSE)}
  if(!missing(image.SR)){
    rgb <- landsatBands(image.SR, c("NIR", "R", "G"))
    terra::plotRGB(rgb, stretch = "lin", ...)
  } else {
    plot(x$EB$surfaceTemperature, main = "Anchor pixels over Ts (K)", ...)
  }
  anchors <- sf::st_transform(x$anchors, terra::crs(x$EB))
  type <- tolower(as.character(anchors$type))
  xy <- sf::st_coordinates(anchors)
  graphics::points(xy[type == "hot", , drop = FALSE], col = "red", pch = 3)
  graphics::points(xy[type == "cold", , drop = FALSE], col = "blue", pch = 4)
  ws <- sf::st_as_sf(x$WeatherStation$location, coords = c("long", "lat"), crs = 4326)
  graphics::points(sf::st_coordinates(sf::st_transform(ws, terra::crs(x$EB))), pch = 13)
  invisible()
}

# Quality checks of the calibration for METRIC.EB(): Ts of the anchors and the
# reference ET fraction at the overpass (before LE is set to 0)
calibrationQC <- function(Rn, G, H, Ts, Ts.anchors, anchors, WeatherStation){
  type <- tolower(trimws(as.character(anchors$type)))
  Ta <- terra::extract(Ts.anchors, terra::vect(sf::st_transform(anchors, terra::crs(Ts))))[, 2]
  lambda <- (2.501 - 0.00236 * (Ts - 273.15)) * 1e6
  ETrF <- 3600 * (Rn - G - H) / lambda / hourlyET(WeatherStation, ET = "ETr")
  share <- function(r){global(r, "mean", na.rm = TRUE)[1, 1]}
  list(Ts.hot = mean(Ta[type == "hot"]), Ts.cold = mean(Ta[type == "cold"]),
       dT.anchors = mean(Ta[type == "hot"]) - mean(Ta[type == "cold"]),
       ETrF = stats::setNames(as.numeric(global(ETrF, stats::quantile,
                                                  probs = c(0.05, 0.5, 0.95, 1), na.rm = TRUE)),
                              c("5%", "50%", "95%", "max")),
       ETrF.above.1.1 = share(ETrF > 1.1), ETrF.below.0 = share(ETrF < 0),
       H.above.Rn.G = share(H > Rn - G))
}

# One message with the calibration; a warning when ETrF > 1.1 on more than
# 10% of the pixels
reportCalibration <- function(k){
  message("Calibration: anchors hot ", round(k$Ts.hot, 1), " K, cold ",
          round(k$Ts.cold, 1), " K (difference ", round(k$dT.anchors, 1),
          " K); ETrF > 1.1 on ", round(100 * k$ETrF.above.1.1, 1),
          "% of the pixels; u200 ", round(k$u200, 1), " m/s; ",
          if(isTRUE(all(k$converged))) "converged" else "not converged",
          ". See summary()")
  if(isTRUE(k$ETrF.above.1.1 > 0.1)){
    warn("ETrF is above 1.1 on ", round(100 * k$ETrF.above.1.1, 1), "% of the ",
         "pixels: the cold anchors may be too warm. Check them with plotAnchors(), ",
         "or use anchors.method = \"percentile\"")
  }
  invisible(k)
}


#' writeRaster method for waterLSEB objects
#' @description Writes the energy balance layers (\code{x$EB}) of a
#' METRIC.EB() result to a file.
#' @param x        waterLSEB object, from METRIC.EB().
#' @param filename output file name, e.g. "EB.tif"
#' @param ...      additional arguments passed to terra::writeRaster(), e.g.
#' overwrite
#' @return The SpatRaster written to the file (x$EB), invisibly, as returned by
#' terra::writeRaster().
#' @author María Victoria Munafó
#' @family LSEB objects related functions
#' @name writeRaster.waterLSEB
#' @aliases writeRaster,waterLSEB,character-method
#' @exportMethod writeRaster
#' @importFrom methods setOldClass setMethod
NULL

setOldClass("waterLSEB")

setMethod("writeRaster", signature(x = "waterLSEB", filename = "character"),
          function(x, filename, ...){
            invisible(terra::writeRaster(x$EB, filename, ...))
          })


