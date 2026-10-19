#' Daily reference ET for each day of a weather station series
#' @description Applies \code{\link{dailyET}} to every day (or to the days
#' requested) of the records of a weather station, to be used as the
#' reference ET of \code{\link{ETperiod}}.
#' @param WeatherStation  \code{waterWeatherStation} object, see
#' \code{\link{read.WSdata}}
#' @param dates     Dates to compute, as Date or "YYYY-MM-DD" text. By default
#' all the days with records in the station data.
#' @param ET        "ETr" (alfalfa reference, default) or "ETo" (grass).
#' @param ...       further arguments of \code{\link{dailyET}}
#' @details Days that do not have all their hourly records give NA (and no
#' warning): \code{\link{ETperiod}} stops on them unless \code{fill.gaps = TRUE}.
#' @return data.frame with the columns date and ET (mm/day)
#' @family ET functions
#' @export
#' @examples
#' csvfile <- system.file("extdata", "apples.csv", package="water")
#' WeatherStation <- read.WSdata(WSdata = csvfile, date.format = "%d/%m/%Y",
#' lat=-35.42222, long= -71.38639, elev=201, height= 2.2, cf=c(1,0.2777778,1),
#' tz = "America/Santiago")
#' dailyETseries(WeatherStation, dates = c("2013-02-14", "2013-02-15"))
dailyETseries <- function(WeatherStation, dates = NULL, ET = "ETr", ...){
  checkChoice(ET, c("ETr", "ETo"), "ET")
  if(!inherits(WeatherStation, "waterWeatherStation")){
    stop("WeatherStation must be a waterWeatherStation object (see read.WSdata())",
         call. = FALSE)
  }
  if(is.null(dates)){
    t <- WeatherStation$hourly$datetime
    # a record ending at 00:00 belongs to the day before (ASCE-EWRI 2005)
    dates <- unique(as.Date(format(t - 1, "%Y-%m-%d")))
  }
  dates <- sort(as.Date(dates))
  quiet <- function(expr) withCallingHandlers(expr,
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage"))
  et <- vapply(as.character(dates), function(d){
    r <- tryCatch(quiet(dailyET(WeatherStation, ET = ET, date = d, ...)),
                  error = function(e) NA_real_)
    as.numeric(r)
  }, numeric(1))
  data.frame(date = dates, ET = unname(et))
}

#' Evapotranspiration of a period from the ETrF of several images
#' @description Seasonal or monthly ET from the reference ET fraction (ETrF) of
#' several satellite dates, as in Allen et al. (2007), Eqs. 58 and 59: ETrF is
#' interpolated for each day of the period, multiplied by the daily reference
#' ET of the day and added up.
#' @param ETrF       SpatRaster with one layer per image date (the ETrF layer
#' of \code{\link{ET24h}} of each image, on the same grid)
#' @param dates      Dates of the images (Date or "YYYY-MM-DD" text), one for
#' each layer of ETrF
#' @param ETr.daily  Daily reference ET (mm/day) of every day of the period: a
#' data.frame with the columns date and ET, as returned by
#' \code{\link{dailyETseries}}
#' @param period     Start and end of the period (two dates, both included).
#' Default: from the first to the last image date. Days before the first and
#' after the last valid image of a pixel use its closest ETrF
#' @param C.rad      Correction for slope and aspect (Eq. 55): a number or a
#' SpatRaster, e.g. the C.rad layer of \code{\link{ET24h}}; 1 on flat terrain.
#' Only one value can be used for the period.
#' @param method     "linear" (default) or "spline" (natural spline, for
#' pixels with 3 or more valid dates) interpolation of ETrF between dates
#' @param fill.gaps  If TRUE, days without ETr.daily (NA or missing) are
#' filled by linear interpolation. If FALSE (default) they are an error.
#' @details The ETrF of a pixel that is masked (NA, e.g. clouds) in one date is
#' interpolated between the closest valid dates, and the number of valid dates
#' is returned in layer n.valid. A pixel without any valid date is NA. ETrF is
#' used as it comes from ET24h(): values below 0 or above 1.1 are not clipped;
#' mask or clip them before if needed.
#' @return SpatRaster with layers ET.period (mm in the period), ETrF.period
#' (Eq. 59, ET.period divided by the reference ET of the period) and n.valid
#' @family ET functions
#' @export
#' @references
#' Allen, R.G., Tasumi, M., Trezza, R., 2007. Satellite-based energy balance for
#' mapping evapotranspiration with internalized calibration (METRIC) - Model.
#' Journal of Irrigation and Drainage Engineering 133, 380-394.
#' @examples
#' f <- terra::rast(nrows = 2, ncols = 2, nlyrs = 2, vals = c(0.2, 0.4, 0.6, NA, 0.8, 0.6, 0.3, 0.5))
#' etr <- data.frame(date = as.Date("2013-02-01") + 0:9, ET = 6)
#' ETperiod(f, dates = c("2013-02-02", "2013-02-08"), ETr.daily = etr)
ETperiod <- function(ETrF, dates, ETr.daily, period = NULL, C.rad = 1,
                     method = c("linear", "spline"), fill.gaps = FALSE){
  method <- match.arg(method)
  if(!inherits(ETrF, "SpatRaster")){stop("ETrF must be a SpatRaster", call. = FALSE)}
  dates <- as.Date(dates)
  if(terra::nlyr(ETrF) != length(dates) || anyNA(dates) || anyDuplicated(dates)){
    stop("dates must be different dates, one for each layer of ETrF", call. = FALSE)
  }
  if(!is.data.frame(ETr.daily) || !all(c("date", "ET") %in% names(ETr.daily))){
    stop("ETr.daily must be a data.frame with the columns date and ET (see dailyETseries())",
         call. = FALSE)
  }
  if(is.null(period)){period <- range(dates)}
  period <- as.Date(period)
  if(length(period) != 2 || anyNA(period) || period[1] > period[2]){
    stop("period must be the start and the end date", call. = FALSE)
  }
  days <- seq(period[1], period[2], by = "day")
  et <- ETr.daily$ET[match(days, as.Date(ETr.daily$date))]
  if(anyNA(et)){
    if(!fill.gaps){
      stop(sum(is.na(et)), " day(s) of the period without ETr.daily (first: ",
           days[which(is.na(et))[1]], "). Provide them or use fill.gaps = TRUE",
           call. = FALSE)
    }
    ok <- !is.na(et)
    if(sum(ok) < 1){stop("ETr.daily has no valid day in the period", call. = FALSE)}
    et <- stats::approx(as.numeric(days[ok]), et[ok], as.numeric(days), rule = 2)$y
  }
  ord <- order(dates)
  ETrF <- ETrF[[ord]]
  dn <- as.numeric(dates[ord])
  dd <- as.numeric(days)
  pixel <- function(x){
    ok <- !is.na(x)
    n <- sum(ok)
    if(n == 0){return(c(NA_real_, NA_real_, 0))}
    f <- if(n == 1) rep(x[ok], length(dd)) else if(method == "spline" && n >= 3){
      stats::spline(dn[ok], x[ok], xout = dd, method = "natural")$y
    } else stats::approx(dn[ok], x[ok], dd, rule = 2)$y
    c(sum(f * et), sum(et), n)
  }
  r <- terra::app(ETrF, pixel)
  out <- c(r[[1]] * C.rad, r[[1]] * C.rad / r[[2]], r[[3]])
  names(out) <- c("ET.period", "ETrF.period", "n.valid")
  out
}
