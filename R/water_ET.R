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
  tempK <- WeatherStation$temp + 273.16
  Rs <- WeatherStation$radiation * 3600 / 1e6
  P <- 101.3*((293-0.0065*WeatherStation$elev)/293)^5.26
  psi <- 0.000665*P
  Delta <- 2503 * exp((17.27*WeatherStation$temp)/
                  (WeatherStation$temp+237.3))/((WeatherStation$temp+237.3)^2)
  ea.sat <- 0.61078*exp((17.269*WeatherStation$temp)/
                          (WeatherStation$temp+237.3))
  ea <- (WeatherStation$RH/100)*ea.sat
  DPV <- ea.sat - ea
  if(is.null(long.z)){
    long.z <- WeatherStation$long.z
    if(is.null(long.z) && inherits(WeatherStation$datetime, "POSIXt")){
      long.z <- tzLongitude(WeatherStation$datetime)
    }
    if(is.null(long.z)){long.z <- WeatherStation$long}  # time zone unknown
  }
  if(ET.instantaneous==TRUE){hours <- hours +0.5}
  Rso <- clearSkyRso(DOY = DOY, hours = hours, lat = WeatherStation$lat,
                     long = WeatherStation$long, long.z = long.z,
                     elev = WeatherStation$elev)
  Rs.Rso <- ifelse(Rs/Rso<=0.3, 0, ifelse(Rs/Rso>=1, 1, Rs/Rso))
  fcd <- ifelse(1.35*Rs.Rso-0.35<=0.05, 0.05, 
                ifelse(1.35*Rs.Rso-0.35<1, 1.35*Rs.Rso-0.35,1))
  Rn.a <- ((1-0.23)*Rs) - (2.042e-10*fcd*(0.34-0.14*(ea^0.5))*tempK^4)
  Ghr <- c(0.1, 0.5)
  if(ET=="ETr"){Ghr <- c(0.04, 0.2)}
  G.day <- Rn.a * ifelse(Rn.a > 0, Ghr[1], Ghr[2])
  wind.2 <- WeatherStation$wind *(4.87/(log(67.8*WeatherStation$height-5.42)))
  windfactors <- c(0.24, 0.96)
  if(ET=="ETr"){windfactors <-  c(0.25, 1.7)}
  windf <- ifelse(WeatherStation$radiation > 0, windfactors[1], windfactors[2])
  ET.hourly <- ((0.408*Delta*(Rn.a-G.day))+(psi*(37/tempK)*wind.2*(DPV)))/
    (Delta+(psi*(1+(windf*wind.2))))
  if(ET=="ETr"){
    ET.hourly <- ((0.408*Delta*(Rn.a-G.day))+(psi*(66/tempK)*wind.2*(DPV)))/
      (Delta+(psi*(1+(windf*wind.2))))
  }
  return(ET.hourly)
}

# Clear-sky solar radiation Rso (MJ m-2 h-1, ASCE-EWRI 2005) for the hour
# ending at hours (decimal hours of the station clock); long.z is the
# longitude of the centre of its time zone
clearSkyRso <- function(DOY, hours, lat, long, long.z, elev){
  dr <- 1 + 0.033*(cos(2*pi*DOY/365))
  delta <- 0.409*sin((2*pi*DOY/365)-1.39)
  phi <- lat*(pi/180)
  b <- 2*pi*(DOY-81)/364
  Sc <- 0.1645*sin(2*b)-0.1255*cos(b)-0.025*sin(b)
  hour.angle <- (pi/12)*((hours-0.5+0.06667*(long-long.z)+Sc)-12)
  w1 <- hour.angle-((pi)/24)
  w2 <- hour.angle+((pi)/24)
  hour.angle.s <- acos(-tan(phi)*tan(delta))
  w1c <- ifelse(w1< -hour.angle.s, -hour.angle.s,
                ifelse(w1>hour.angle.s, hour.angle.s, ifelse(w1>w2, w2, w1)))
  w2c <- ifelse(w2< -hour.angle.s, -hour.angle.s,
                ifelse(w2>hour.angle.s, hour.angle.s, w2))
  Beta <- asin((sin(phi)*sin(delta)+cos(phi)*cos(delta)*cos(hour.angle)))
  Ra <- ifelse(Beta <= 0, 1e-45, ((12/pi)*4.92*dr)*
                 (((w2c-w1c)*sin(phi)*sin(delta))+(cos(phi)*cos(delta)*
                                                     (sin(w2)-sin(w1)))))
  (0.75+2e-5*elev)*Ra
}

#' Calculates daily ET from a surface energy balance and Weather Station
#' @description
#' Calculates the daily (24 h) actual evapotranspiration: the instantaneous ET
#' (from LE = Rn - G - H) divided by the hourly alfalfa reference ET (ETr, see
#' hourlyET()) at the satellite overpass gives the reference ET fraction, which
#' is multiplied by the daily reference ET.
#' @param Rn             Net radiation (W/m2). See \code{\link{netRadiation}}
#' @param G              Soil Heat Flux (W/m2). See \code{\link{soilHeatFlux}}
#' @param H              Sensible Heat Flux (W/m2). See \code{\link{calcH}}
#' @param Ts             Land surface temperature (K). See \code{\link{surfaceTemperature}}
#' @param WeatherStation WeatherStation data at the satellite overpass.
#' Can be a waterWeatherStation object calculated using read.WSdata and MTL file
#' @param ETr.daily      daily reference ET (mm/day) for the day of the image,
#' e.g. from dailyET() (whose default is the alfalfa reference, ETr)
#' @param C.rad          correction term used in sloping terrain to correct for 
#' variation in 24 h versus instantaneous energy availability. See Allen (2007)
#' @return SpatRaster with daily (24 h) actual evapotranspiration (mm/day);
#' negative values are set to 0. The result is also plotted.
#' @family ET functions
#' @author Guillermo Federico Olmedo
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @importFrom grDevices colorRampPalette
#' @export
ET24h <- function(Rn, G, H, Ts, WeatherStation, ETr.daily, C.rad=1){
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
  saveLoadClean(imagestack = ET.24, 
                file = "ET24", overwrite=TRUE)
}


#' Calculates daily reference ET as the sum of hourly Penman-Monteith values
#' @param WeatherStation waterWeatherStation object, see \code{\link{read.WSdata}}. Other
#' inputs are not implemented
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
#' @param MTL      Optional path to the Landsat metadata (MTL) file. Its
#' acquisition date selects the day (see date)
#' @param date     "auto" (default) uses the acquisition date in the MTL file
#' or, without MTL, the date of the conditions at the overpass of
#' WeatherStation (at.sat); otherwise a date (e.g. "2016-02-09") present in the
#' weather station data. Without MTL, date or at.sat it is an error
#' @return Numeric. Daily reference ET (mm/day): the sum of the 24 hourly values
#' computed with hourlyET().
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
    if(date != "auto"){
      day <- as.character(as.Date(date))
    } else {
      if(!is.null(MTL)){
        MTL <- readLines(MTL, warn=FALSE)
        date.line <- grep("DATE_ACQUIRED",MTL,value=TRUE)
        sat.date <-regmatches(date.line,regexec(text=date.line,
                                                pattern="([0-9]{4})(-)([0-9]{2})(-)([0-9]{2})"))[[1]][1]
        sat.date <- strptime(sat.date,
                             format = "%Y-%m-%d", tz="GMT")
        day <- as.character(as.Date(sat.date))
      } else if(!is.null(WeatherStation$at.sat)){
        day <- format(WeatherStation$at.sat$datetime, "%Y-%m-%d")
      } else {
        stop("dailyET() needs the day of the image: pass MTL or date (e.g. ",
             "date = \"2016-02-09\"), or read the weather station with MTL",
             call. = FALSE)
      }
    }
    # day in the time zone of the station records
    data <- hourly[as.character(trunc(hourly$datetime, "days")) == day,]
    if(nrow(data) == 0){
      stop("No weather station records on ", day, call. = FALSE)
    }
    ET.daily <- vector()
    for(i in 1:24){
      ET.daily <- c(ET.daily, hourlyET(data[i,], lat=lat,
                                       long = long, elev=elev, ET=ET,
                                       height = height, long.z = long.z))
    }
  } else {
    stop("WeatherStation must be a waterWeatherStation object (see read.WSdata())",
         call. = FALSE)
  }
  return(sum(ET.daily, na.rm=T))
}
