#' Calculates hourly reference ET with the ASCE Penman-Monteith equation
#' @description
#' Calculates the hourly (or instantaneous) reference evapotranspiration with
#' the ASCE standardized Penman-Monteith equation, for a short (ETo, grass) or a
#' tall (ETr, alfalfa) reference crop.
#' @param WeatherStation     data.frame with the needed fields (see example),
#' or a waterWeatherStation object with the conditions at the satellite
#' overpass; in that case they are used and ET.instantaneous is set to TRUE
#' @param hours              time of the day in decimal hours (24 h clock),
#' e.g. 10.5 for 10:30
#' @param DOY                day of year
#' @param long.z             longitude of the centre of the time zone of the
#' weather station clock, in decimal degrees (negative west), for the solar
#' time correction. By default it is computed from the UTC offset of the
#' station records (15 degrees per hour, e.g. -45 for UTC-3, including daylight
#' saving time): from WeatherStation$long.z (see \code{\link{getDataWS}}) or from the time
#' zone of WeatherStation$datetime. If the time zone is unknown, the station
#' longitude is used (no correction)
#' @param ET.instantaneous   Logical. TRUE to calculate instantaneous ET
#' instead of hourly ET. See Details.
#' @param ET       "ETo" for short crops, similar to clipped, cool-season
#' grass; or "ETr" (default) for tall crops, similar to 0.5 m tall full-cover
#' alfalfa.
#' @param height             weather station sensors height in meters
#' @param lat                latitude of weather station in decimal degrees. 
#' Negative values for south latitude
#' @param long               longitude of weather station in decimal degrees. 
#' Negative values for west longitude
#' @param elev               elevation of weather station in meters
#' @details
#' The only difference on instantaneous ET is how the hour is interpreted. With
#' FALSE, and for example at 11:00, ET is calculated between 10:00 and 11:00;
#' with TRUE, ET is calculated at 11:00.
#'
#' The equations follow ASCE-EWRI (2005): Rs/Rso is limited to 0.3-1 (Eq.
#' 45), the extraterrestrial radiation is integrated between sunrise and
#' sunset (Eqs. 48 and 56-58), the daytime and nighttime coefficients Cd and G
#' are chosen by Rn > 0 (Table 1), and T + 273 is used in the aerodynamic term
#' and T + 273.16 in the net long-wave radiation. For a single period at night
#' (sun angle below 0.3 rad) the cloudiness function cannot be computed and
#' Rs/Rso = 0.3 is used; dailyET() carries the value of the last sunny period
#' over the night instead (Eq. 46).
#' @return Numeric. Hourly (or, with ET.instantaneous = TRUE, instantaneous)
#' reference ET (mm/h).
#' @family ET functions
#' @author Guillermo Federico Olmedo
#' @examples 
#' WeatherStation  <- data.frame(wind=4.72,
#'                               RH=59, 
#'                               temp=24.3,
#'                               radiation=675, 
#'                               height=2.2, 
#'                               lat=-35.37, 
#'                               long=-71.5946,
#'                               elev=124)
#'   hourlyET(WeatherStation, hours=10.5, DOY=363, long.z=-71.635)
#' @export
#' @references
#' ASCE-EWRI (2005). The ASCE Standardized Reference Evapotranspiration
#' Equation. Report of the ASCE-EWRI Task Committee on Standardization of
#' Reference Evapotranspiration.
hourlyET <- function(WeatherStation, hours, DOY, long.z=NULL,
                     ET.instantaneous=FALSE, ET="ETr", height=2, lat, long, elev){
  checkChoice(ET, c("ETr", "ETo"), "ET")
  if(inherits(WeatherStation, "waterWeatherStation")){
    if(!is.null(WeatherStation$at.sat)){
      WeatherStation <- getDataWS(WeatherStation)
      ET.instantaneous = TRUE
      DOY <- WeatherStation$DOY
      hours <-  WeatherStation$hours}
  } else {
    if("datetime" %in% names(WeatherStation)){
      date <- as.POSIXlt(WeatherStation$datetime, format="%Y-%m-%d %H:%M:%S")
      DOY=date$yday+1
      hours=date$hour + date$min/60 + date$sec/3600}
    if(!missing(lat)){WeatherStation$lat <- lat}
    if(!missing(long)){WeatherStation$long <- long}
    if(!missing(elev)){WeatherStation$elev <- elev}
    if(!missing(height)){WeatherStation$height <- height}
    if(length(WeatherStation$DOY) != 0 & length(WeatherStation$hours) != 0 &
       missing(DOY) & missing(hours)){
      DOY <- WeatherStation$DOY
      hours <- WeatherStation$hours
    } 
  }
  if(is.null(long.z)){
    long.z <- WeatherStation$long.z
    if(is.null(long.z) && inherits(WeatherStation$datetime, "POSIXt")){
      long.z <- tzLongitude(WeatherStation$datetime)
    }
    if(is.null(long.z)){long.z <- WeatherStation$long}  # time zone unknown
  }
  if(ET.instantaneous==TRUE){hours <- hours +0.5}
  asceHourly(WeatherStation, hours = hours, DOY = DOY, long.z = long.z, ET = ET)$ET
}

# ASCE-EWRI (2005) standardized reference ET for the hours (periods) ending at
# hours, with the station data in the data.frame WS (temp, RH, wind,
# radiation, lat, long, elev, height). fcd: cloudiness function to use, e.g.
# carried over the night by dailyET() (Eq. 46); if NULL it is computed from
# Rs/Rso of each period. Returns a list with ET (mm/h), fcd and beta (sun
# angle at the middle of the period, radians)
asceHourly <- function(WS, hours, DOY, long.z, ET = "ETr", fcd = NULL){
  temp <- WS$temp
  Rs <- WS$radiation * 3600 / 1e6
  P <- 101.3*((293-0.0065*WS$elev)/293)^5.26                  # Eq. 34
  psi <- 0.000665*P                                           # Eq. 35
  Delta <- 2503 * exp((17.27*temp)/(temp+237.3))/((temp+237.3)^2)  # Eq. 36
  ea.sat <- 0.6108*exp((17.27*temp)/(temp+237.3))             # Eq. 37
  ea <- (WS$RH/100)*ea.sat
  DPV <- ea.sat - ea
  sun <- asceRa(DOY = DOY, hours = hours, lat = WS$lat, long = WS$long, long.z = long.z)
  Rso <- (0.75+2e-5*WS$elev)*sun$Ra                           # Eq. 47
  if(is.null(fcd)){
    # Rs/Rso limited to 0.3-1 (Eq. 45); at night (Rso = 0) the minimum
    ratio <- ifelse(Rso > 0, Rs/Rso, 0.3)
    fcd <- 1.35*pmin(pmax(ratio, 0.3), 1) - 0.35
  }
  Rnl <- 2.042e-10*fcd*(0.34-0.14*sqrt(ea))*(temp + 273.16)^4  # Eq. 44
  Rn <- (1-0.23)*Rs - Rnl                                     # Eqs. 42-43
  # daytime when Rn > 0 (Table 1)
  day <- Rn > 0
  if(ET == "ETr"){
    G <- Rn * ifelse(day, 0.04, 0.2); Cn <- 66; Cd <- ifelse(day, 0.25, 1.7)
  } else {
    G <- Rn * ifelse(day, 0.1, 0.5); Cn <- 37; Cd <- ifelse(day, 0.24, 0.96)
  }
  wind.2 <- WS$wind *(4.87/(log(67.8*WS$height-5.42)))         # Eq. 33
  ET.hourly <- ((0.408*Delta*(Rn-G))+(psi*(Cn/(temp+273))*wind.2*(DPV)))/
    (Delta+(psi*(1+(Cd*wind.2))))                             # Eq. 1
  list(ET = ET.hourly, fcd = fcd, beta = sun$beta)
}

# Extraterrestrial radiation Ra (MJ m-2 h-1, ASCE-EWRI 2005 Eq. 48) for the
# hour ending at hours (decimal hours of the station clock), and the sun
# angle beta (radians, Eq. 62) at its middle; long.z is the longitude of the
# centre of the time zone of the clock
asceRa <- function(DOY, hours, lat, long, long.z){
  dr <- 1 + 0.033*(cos(2*pi*DOY/365))
  delta <- 0.409*sin((2*pi*DOY/365)-1.39)
  phi <- lat*(pi/180)
  b <- 2*pi*(DOY-81)/364
  Sc <- 0.1645*sin(2*b)-0.1255*cos(b)-0.025*sin(b)
  hour.angle <- (pi/12)*((hours-0.5+0.06667*(long-long.z)+Sc)-12)
  w1 <- hour.angle-((pi)/24)
  w2 <- hour.angle+((pi)/24)
  hour.angle.s <- acos(-tan(phi)*tan(delta))
  # integration limits within sunrise and sunset (Eqs. 56-58)
  w1c <- ifelse(w1< -hour.angle.s, -hour.angle.s,
                ifelse(w1>hour.angle.s, hour.angle.s, w1))
  w2c <- ifelse(w2< -hour.angle.s, -hour.angle.s,
                ifelse(w2>hour.angle.s, hour.angle.s, w2))
  w1c <- pmin(w1c, w2c)
  beta <- asin((sin(phi)*sin(delta)+cos(phi)*cos(delta)*cos(hour.angle)))
  Ra <- ((12/pi)*4.92*dr)*(((w2c-w1c)*sin(phi)*sin(delta))+
                             (cos(phi)*cos(delta)*(sin(w2c)-sin(w1c))))
  list(Ra = pmax(Ra, 0), beta = beta)
}

# Clear-sky solar radiation Rso (MJ m-2 h-1, ASCE-EWRI 2005 Eq. 47) for the
# hour ending at hours (decimal hours of the station clock)
clearSkyRso <- function(DOY, hours, lat, long, long.z, elev){
  (0.75+2e-5*elev)*asceRa(DOY = DOY, hours = hours, lat = lat, long = long,
                          long.z = long.z)$Ra
}

#' Calculates daily ET from a surface energy balance and Weather Station
#' @description
#' Calculates the daily (24 h) actual evapotranspiration: the instantaneous ET
#' (from LE = Rn - G - H) divided by the instantaneous alfalfa reference ET
#' (ETr, see hourlyET()) at the satellite overpass gives the reference ET
#' fraction (ETrF), which is multiplied by the daily reference ET.
#' @param Rn             Net radiation (W/m2), see \code{\link{netRadiation}}; or
#' the result of METRIC.EB(), which supplies Rn, G, H, Ts and the weather
#' station. See Value
#' @param G              Soil Heat Flux (W/m2). See \code{\link{soilHeatFlux}}
#' @param H              Sensible Heat Flux (W/m2). See \code{\link{calcH}}
#' @param Ts             Land surface temperature (K). See \code{\link{surfaceTemperature}}
#' @param WeatherStation WeatherStation data at the satellite overpass.
#' Can be a waterWeatherStation object calculated using read.WSdata and MTL file
#' (with read.WSdata(..., overpass.radiation = "clear-sky"), the clear-sky
#' radiation is used, as in calcH()). Optional with a METRIC.EB() result
#' @param ETr.daily      daily reference ET (mm/day) for the day of the image,
#' e.g. from dailyET() (whose default is the alfalfa reference, ETr). If
#' missing, dailyET(WeatherStation) is used (the day of the overpass)
#' @param C.rad          correction term used in sloping terrain to correct for
#' variation in 24 h versus instantaneous energy availability. See Allen (2007)
#' @details
#' ETrF is a crop coefficient on the alfalfa reference basis. To compare it
#' with a crop coefficient on the grass reference basis (Kc, FAO-56), multiply
#' it by the ratio of the daily alfalfa and grass reference ET,
#' dailyET(WeatherStation) / dailyET(WeatherStation, ET = "ETo").
#'
#' LE = Rn - G - H is not limited, so ETinst and ETrF are negative where H >
#' Rn - G (usually hot, dry pixels, or where the calibration does not fit);
#' layer HgtRnG flags them. METRIC.EB() sets LE to 0 there.
#' @return With Rn, G, H and Ts: a SpatRaster with daily (24 h) actual
#' evapotranspiration (mm/day) in a layer named ET24; negative values are set
#' to 0. With a METRIC.EB() result: a SpatRaster with layers ETinst
#' (instantaneous ET at the overpass, mm/h), ETrF (reference ET fraction,
#' unitless), ET24 (mm/day, negative values set to 0) and HgtRnG (1 where H >
#' Rn - G, else 0). ET24 is also plotted.
#' @family ET functions
#' @author Guillermo Federico Olmedo
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @importFrom grDevices colorRampPalette
#' @export
ET24h <- function(Rn, G, H, Ts, WeatherStation, ETr.daily, C.rad=1){
  full <- inherits(Rn, "waterLSEB")  # a METRIC.EB() result
  if(full){
    EB <- Rn
    if(missing(WeatherStation)){WeatherStation <- EB$WeatherStation}
    Rn <- EB$EB$NetRadiation; G <- EB$EB$SoilHeat; H <- EB$EB$SensibleHeat
    Ts <- EB$EB$surfaceTemperature
  }
  if(missing(ETr.daily)){
    if(!inherits(WeatherStation, "waterWeatherStation")){
      stop("ETr.daily is missing: pass it (e.g. from dailyET()) or a ",
           "waterWeatherStation object as WeatherStation", call. = FALSE)
    }
    ETr.daily <- dailyET(WeatherStation)
  }
  if(inherits(WeatherStation, "waterWeatherStation")){
    WeatherStation <- getDataWS(WeatherStation)
  }
  LE = Rn - G - H
  ET.inst <- 3600*LE/((2.501 - 0.00236 * (Ts - 273.15)) * (1e6))
  ETr.hourly <- hourlyET(WeatherStation, WeatherStation$hours, WeatherStation$DOY,
                         ET.instantaneous = TRUE, ET = "ETr")
  ETr.Fr <- ET.inst/ETr.hourly
  ET.24 <- ETr.Fr * ETr.daily * C.rad
  ET.24[ET.24 < 0]  <- 0
  #ET.24[ET.24 > quantile(ET.24, 0.9)] <- quantile(ET.24, 0.9)
  rgb.palette <- grDevices::colorRampPalette(c("red3","snow2","blue"),  space = "rgb")
  plot(ET.24, col = rgb.palette(50), range = c(0, ceiling(ETr.daily*1.5)),
       fill_range = TRUE, xlab = "Easting", ylab = "Northing",
       main = "24-Hour Evapotranspiration (mm/day)", maxcell = ncell(ET.24) * 0.3)
  if(!full){
    return(saveLoadClean(imagestack = ET.24, file = "ET24", overwrite=TRUE))
  }
  result <- c(ET.inst, ETr.Fr, ET.24, H > Rn - G)
  saveLoadClean(imagestack = result, stack.names = c("ETinst", "ETrF", "ET24", "HgtRnG"),
                file = "ET24", overwrite = TRUE)
}


#' Calculates daily reference ET as the sum of hourly Penman-Monteith values
#' @description
#' Sums the hourly reference ET (ASCE standardized Penman-Monteith equation,
#' see \code{\link{hourlyET}}) of the hourly weather station values of one day,
#' in the time zone of the station records: the day of the satellite overpass
#' or another day.
#' @param WeatherStation waterWeatherStation object, see \code{\link{read.WSdata}}.
#' Required (other inputs stop)
#' @param DOY      Not used
#' @param height   weather station sensors height in meters. If missing, taken
#' from WeatherStation
#' @param lat      latitude in decimal degrees of the weather station. If
#' missing, taken from WeatherStation
#' @param long     longitude in decimal degrees of the weather station. If
#' missing, taken from WeatherStation
#' @param elev     elevation in meters of the weather station. If missing,
#' taken from WeatherStation
#' @param ET       "ETo" for short crops, similar to clipped, cool-season
#' grass; or "ETr" (default) for tall crops, similar to 0.5 m tall full-cover
#' alfalfa.
#' @param long.z   longitude of the centre of the time zone of the weather
#' station clock (decimal degrees, negative west). By default it is computed
#' from the time zone of the records; see \code{\link{hourlyET}}
#' @param MTL      path to the Landsat metadata (*_MTL.txt) file. Optional:
#' the day of the overpass in the time zone of the station selects the day
#' (see date)
#' @param date     "auto" (default) uses the acquisition date in the MTL file
#' or, without MTL, the date of the conditions at the overpass of
#' WeatherStation (at.sat); otherwise a date (e.g. "2016-02-09") present in the
#' weather station data, taken in the time zone of the station records.
#' Without MTL, date or at.sat, dailyET() stops
#' @details
#' As in ASCE-EWRI (2005), the day is the 24 hourly periods ending at 01:00 to
#' 24:00 (00:00 of the next day); 23 or 25 periods on days when daylight saving
#' time starts or ends. The cloudiness function of the net long-wave radiation
#' is computed while the sun is higher than 0.3 rad, and its last value is
#' carried over dusk, night and dawn (Eq. 46); before the first sunny period
#' the previous evening is used, or the last sunny period of the day if the
#' records do not include the previous day.
#' @return Numeric. Daily reference ET (mm/day): the sum of the hourly values
#' of the day. If the records do not cover the whole day, NA with a warning
#' (only the last period, ending at 24:00, may be missing: it is then left
#' out, with a message).
#' @family ET functions
#' @author Guillermo Federico Olmedo
#' @export
#' @references
#' ASCE-EWRI (2005). The ASCE Standardized Reference Evapotranspiration
#' Equation. Report of the ASCE-EWRI Task Committee on Standardization of
#' Reference Evapotranspiration.
#' @examples
#' csvfile <- system.file("extdata", "apples.csv", package="water")
#'
#' WeatherStation <- read.WSdata(WSdata = csvfile, date.format = "%d/%m/%Y",
#' lat=-35.42222, long= -71.38639, elev=201, height= 2.2, cf=c(1,0.2777778,1),
#' tz = "America/Santiago")
#'
#' dailyET(WeatherStation = WeatherStation, ET="ETo", date = "2013-02-15")
#' 
dailyET <- function(WeatherStation, DOY, height, lat, long, elev, ET="ETr", 
                    long.z=NULL, date = "auto", MTL){
  checkChoice(ET, c("ETr", "ETo"), "ET")
  if(inherits(WeatherStation, "waterWeatherStation")){
    if(missing(height)){height <- WeatherStation$location$height}
    if(missing(lat)){lat <- WeatherStation$location$lat}
    if(missing(long)){long <- WeatherStation$location$long}
    if(missing(elev)){elev <- WeatherStation$location$elev}
    
    ## Day of the image: date, the acquisition date in the MTL file or the
    ## date of the conditions at the overpass (at.sat)
    MTL <- findMTL(if(!missing(MTL)) MTL, required = FALSE)
    hourly <- WeatherStation$hourly
    if(!identical(date, "auto")){
      day <- as.character(as.Date(date))  # a Date or a character date
    } else {
      if(!is.null(MTL)){
        # local day of the overpass: DATE_ACQUIRED and SCENE_CENTER_TIME are
        # UTC (east of about UTC+10 the local day is the next one)
        MTL <- readLines(MTL, warn=FALSE)
        date.line <- grep("DATE_ACQUIRED",MTL,value=TRUE)
        time.line <- grep("SCENE_CENTER_TIME",MTL,value=TRUE)
        sat.date <-regmatches(date.line,regexec(text=date.line,
                                                pattern="([0-9]{4})(-)([0-9]{2})(-)([0-9]{2})"))[[1]][1]
        sat.time <- regmatches(time.line, regexec(text = time.line,
                                                  pattern = "([0-9]{2}):([0-9]{2}):([0-9]{2})"))[[1]][1]
        sat.datetime <- as.POSIXct(paste(sat.date, sat.time), format = "%Y-%m-%d %H:%M:%S",
                                   tz = "UTC")
        day <- format(sat.datetime, "%Y-%m-%d",
                      tz = attr(as.POSIXlt(hourly$datetime[1]), "tzone")[1])
      } else if(!is.null(WeatherStation$at.sat)){
        day <- format(WeatherStation$at.sat$datetime, "%Y-%m-%d")
      } else {
        stop("dailyET() needs the day of the image: pass MTL or date (e.g. ",
             "date = \"2016-02-09\"), or read the weather station with MTL",
             call. = FALSE)
      }
    }
    # the day in the time zone of the station records: the hourly periods
    # ending at 01:00 to 24:00 (ASCE-EWRI 2005); 23 or 25 with daylight saving
    tz <- attr(as.POSIXlt(hourly$datetime[1]), "tzone")[1]
    if(is.null(tz)){tz <- ""}
    start <- as.POSIXct(paste(day, "00:00:00"), tz = tz)
    end <- as.POSIXct(paste(as.Date(day) + 1, "00:00:00"), tz = tz)
    expected <- round(as.numeric(difftime(end, start, units = "hours")))
    t <- hourly$datetime
    target <- t > start & t <= end
    if(!any(target)){
      stop("No weather station records on ", day, call. = FALSE)
    }
    # with the previous day, if available, for the cloudiness at night
    rows <- t > start - 86400 & t <= end
    d <- hourly[rows, ]
    d$lat <- lat; d$long <- long; d$elev <- elev; d$height <- height
    lt <- as.POSIXlt(d$datetime)
    hours <- lt$hour + lt$min/60 + lt$sec/3600
    DOY <- lt$yday + 1
    lz <- if(is.null(long.z)) tzLongitude(d$datetime) else long.z
    first <- asceHourly(d, hours = hours, DOY = DOY, long.z = lz, ET = ET)
    # cloudiness function: from Rs/Rso while the sun is higher than 0.3 rad,
    # then carried over dusk, night and dawn (Eq. 46). Before the first such
    # period, the previous evening or else the last one of the day
    tgt <- target[rows]
    sunny <- first$beta > 0.3
    fcd <- first$fcd
    current <- if(any(sunny & tgt)) utils::tail(fcd[sunny & tgt], 1) else fcd[1]
    for(i in seq_along(fcd)){
      if(sunny[i]){current <- fcd[i]} else {fcd[i] <- current}
    }
    ET.hourly <- asceHourly(d, hours = hours, DOY = DOY, long.z = lz, ET = ET,
                            fcd = fcd)$ET[tgt]
    n <- sum(tgt)
    if(n < expected){
      if(n == expected - 1 && !any(t == end)){
        message("The record of 24:00 on ", day, " (00:00 of the next day) is ",
                "missing: the last hour (at night) is not included")
      } else {
        warning("Only ", n, " of ", expected, " hours on ", day, ": daily ", ET,
                " is NA. Provide the whole day, or the daily reference ET from ",
                "another source", call. = FALSE)
        return(NA_real_)
      }
    }
  } else {
    stop("WeatherStation must be a waterWeatherStation object (see read.WSdata())",
         call. = FALSE)
  }
  return(sum(ET.hourly))
}
