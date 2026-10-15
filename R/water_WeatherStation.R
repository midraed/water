#' Prepares weather station data
#' @param WSdata             csv file with weather station data or data.frame
#' @param ...                additional parameter to pass to read.csv()
#' @param height             weather station sensors height in meters
#' @param lat                latitude of weather station in decimal degrees. 
#' Negative values for south latitude
#' @param long               longitude of weather station in decimal degrees. 
#' Negative values for west longitude
#' @param elev               elevation of weather station in meters
#' @param columns            vector with the column numbers in WSdata for the 
#'                           date, time, radiation, wind, RH, temperature and 
#'                           rain. If date and time are in the same column, the 
#'                           column number has to be the same. Names in this
#'                           vector are ignored and are presented on Usage and 
#'                           examples only as a reference.
#' @param date.format        date format. See strptime format argument.
#' @param time.format        time format. See strptime format argument.
#' @param datetime.format    datetime format. See strptime format argument.
#' @param tz                 timezone of the weather station dates. If not 
#' present assumes the same timezone as the computer running the code. See 
#' strptime for details. 
#' @param cf                 conversion factor to convert radiation, wind, and 
#' temperature to W/m2; m/s and Celsius. See Details.
#' @param MTL                Metadata file. If not provided will look for one on
#' working directory. If provided or present will calculate weather conditions
#' on satellite overpass.
#' @details 
#' For cf, if your data is in W/m2, km/h and Celsius (radiation, wind, 
#' temperature), cf should be: cf = c(1,0.2777778,1)
#' @examples 
#' csvfile <- system.file("extdata", "apples.csv", package="water")
#' MTLfile <- system.file("extdata", "L7.MTL.txt", package="water")
#' WS <- read.WSdata(WSdata = csvfile, date.format = "%d/%m/%Y", 
#'                   lat=-35.42222, long= -71.38639, elev=201, height= 2.2,
#'                   columns=c("date" = 1, "time" = 2, "radiation" = 3,
#'                   "wind" = 4, "RH" = 6, "temp" = 7, "rain" = 8), 
#'                   MTL = MTLfile, tz = "America/Santiago")
#' print(WS)
#' plot(WS, alldata=FALSE)
#' plot(WS, alldata=TRUE)
#' @author Guillermo Federico Olmedo
#' @return waterWeatherStation object, with data.frames with all data, hourly
#' data and conditions at satellite flyby.
#' @family Weather station related functions
#' @seealso \code{\link{read.WSdata2}} for the equivalent using read.csv2()
#' @references 
#' Landsat 7 Metadata example file available from the U.S. Geological Survey.
#' Weather Station example file courtesy of CITRA, Universidad de Talca, Chile
#' @export
read.WSdata <- function(WSdata, ..., height = 2.2, lat, long, elev,
                        columns=c("date" = 1, "time" = 2, "radiation" = 3,
                                  "wind" = 4, "RH" = 6, "temp" = 7, "rain" = 8),
                        date.format = "%Y-%m-%d", time.format = "%H:%M:%S", 
                        datetime.format = "%Y-%m-%d %H:%M:%S", tz = "",
                        cf = c(1, 1, 1), MTL){
  if(tz==""){warning(paste('As tz = "", assuming the weather station time zone is', Sys.timezone()))
  } else if(!tz %in% OlsonNames()){
    stop("Unknown time zone: \"", tz, "\". Use a name from OlsonNames(), e.g. ",
         "\"America/Santiago\", or \"Etc/GMT+6\" for 6 hours behind UTC")
  }
  if("temp" %in% columns){warning("The parameter columns has changed and the old
  sintaxis is now deprecated. In future versions of water package ONLY the 
  new (numeric and simple!) sintaxis will work.")}
  if("pp" %in% columns){columns[columns == "pp"] <- "rain"}  ## TODO: deprecated "pp"
  if(is.character(WSdata)){
    if(!file.exists(WSdata)){stop("Weather station file not found: ", WSdata)}
    WSdata <- utils::read.csv(WSdata, ...)
  } else if(is.data.frame(WSdata)){
    WSdata <- as.data.frame(WSdata)  # e.g. a tibble
  } else stop("WSdata should be a character string with the name of the csv file
  or a data.frame with the WS data")
  num <- function(x, name){
    if(!is.numeric(x)){
      stop("The weather station column for ", name, " is not numeric (first ",
           "value: \"", x[1], "\"). If the file uses decimal commas, use ",
           "read.WSdata2() or pass dec = \",\"")
    }
    x
  }
  result <- list()
  result$location <- data.frame(lat=lat, long=long, elev=elev, height=height)
  
  ### Old method for columns : the character method!
  if(is.character(columns)){
    if("date" %in% columns & "time" %in% columns){
      raw  <- paste(WSdata[, which(columns == "date")],
                    WSdata[, which(columns == "time")])
      fmt <- paste(date.format, time.format)
    } else {
      if("datetime" %in% columns){
        raw  <- WSdata[, which(columns == "datetime")]
        fmt <- datetime.format
      } else {stop("date and time or datetime are needed columns")}
    }
    if("rain" %in% columns){rain = num(WSdata[, which(columns == "rain")], "rain")}
    else rain = NA
    radiation = num(WSdata[, which(columns == "radiation")], "radiation") * cf[1]
    wind =  num(WSdata[, which(columns == "wind")], "wind") * cf[2]
    RH =  num(WSdata[, which(columns == "RH")], "RH")
    temp =  num(WSdata[, which(columns == "temp")], "temp") * cf[3]
  }
  
  ### New method for columns : the numeric method!
  # vector with the column numbers in WSdata for the 
  # date, time, radiation, wind, RH, temp and rain. If 
  # date and time are in the same column, the column 
  # number has to be the same.
  if(is.numeric(columns)){
    if(columns[1] != columns[2]){
      raw  <- paste(WSdata[, columns[1]],
                    WSdata[, columns[2]])
      fmt <- paste(date.format, time.format)
    } else {
      raw  <- WSdata[, columns[1]]
      fmt <- datetime.format
    }
    if(!is.na(columns[7])){rain = num(WSdata[, columns[7]], "rain")}
    else rain = NA
    radiation = num(WSdata[, columns[3]], "radiation") * cf[1]
    wind =  num(WSdata[, columns[4]], "wind") * cf[2]
    RH =  num(WSdata[, columns[5]], "RH")
    temp =  num(WSdata[, columns[6]], "temp") * cf[3]

  }
  datetime <- strptime(raw, format = fmt, tz = tz)
  bad <- is.na(datetime)
  if(all(bad)){
    stop("No date could be read with the format \"", fmt, "\" (first value: \"",
         raw[1], "\"). Check date.format, time.format or datetime.format ",
         "(see strptime)")
  }
  if(any(bad)){
    warning(sum(bad), " records with dates that could not be read with the ",
            "format \"", fmt, "\" were removed (first: \"", raw[bad][1], "\")")
    datetime <- datetime[!bad]
    radiation <- radiation[!bad]; wind <- wind[!bad]; RH <- RH[!bad]
    temp <- temp[!bad]
    if(length(rain) > 1){rain <- rain[!bad]}
  }
  
  ea = (RH/100)*0.6108*exp((17.27*temp)/(temp+237.3))
  WSdata <- data.frame(datetime=datetime, radiation=radiation, wind=wind,
                       RH=RH, ea=ea, temp=temp, rain=rain)
  result$alldata <- WSdata
  ## Daily
  WSdata$date <- as.Date(WSdata$datetime, tz = tz)
  result$daily <- data.frame(date=unique(WSdata$date), 
                             radiation_sum=tapply(WSdata$radiation, WSdata$date, sum),
                             wind_mean=tapply(WSdata$wind, WSdata$date, mean),
                             RH_mean=tapply(WSdata$RH, WSdata$date, mean),
                             temp_mean=tapply(WSdata$temp, WSdata$date, mean),
                             temp_max=tapply(WSdata$temp, WSdata$date, max),
                             temp_min=tapply(WSdata$temp, WSdata$date, min),
                             ea_mean=tapply(WSdata$ea, WSdata$date, mean),
                             rain_sum=tapply(WSdata$rain, WSdata$date, sum))
  ## Hourly

  result$hourly <- list()
  datetime <- as.POSIXlt(WSdata$datetime)
  if(datetime[1]$min==0 & datetime[1]$sec==0){
    first = datetime[1]
  } else {
    first <- datetime[1] + 3600 - datetime$min[1] * 60 - datetime$sec[1]
  }
    
  sequence <- seq.POSIXt(from=first, to = tail(datetime,1),by= "1 hour")
  # Time interpolation
  radiationFun <- stats::approxfun(result$alldata$datetime, result$alldata$radiation)
  windFun <- stats::approxfun(result$alldata$datetime, result$alldata$wind)
  RHFun <- stats::approxfun(result$alldata$datetime, result$alldata$RH)
  eaFun <- stats::approxfun(result$alldata$datetime, result$alldata$ea)
  tempFun <- stats::approxfun(result$alldata$datetime, result$alldata$temp)
  if(!all(is.na(rain))){
    rainFun <- stats::approxfun(result$alldata$datetime, result$alldata$rain)
  } else { rainFun <- function(seq){seq <- NA}}
  result$hourly <- data.frame(datetime = sequence,
                               radiation = radiationFun(sequence),
                               wind = windFun(sequence),
                               RH = RHFun(sequence),
                               ea = eaFun(sequence),
                               temp = tempFun(sequence),
                               rain = rainFun(sequence))
    ## Join with satellite data
  MTL <- findMTL(if(!missing(MTL)) MTL, required = FALSE)
  if(!is.null(MTL)){
    MTL <- readLines(MTL, warn=FALSE)
    time.line <- grep("SCENE_CENTER_TIME",MTL,value=TRUE)
    date.line <- grep("DATE_ACQUIRED",MTL,value=TRUE)
    sat.time <-regmatches(time.line,regexec(text=time.line,
           pattern="([0-9]{2})(:)([0-9]{2})(:)([0-9]{2})(.)([0-9]{2})"))[[1]][1]
    sat.date <-regmatches(date.line,regexec(text=date.line,
                        pattern="([0-9]{4})(-)([0-9]{2})(-)([0-9]{2})"))[[1]][1]
    sat.datetime <- strptime(paste(sat.date, sat.time),
                             format = "%Y-%m-%d %H:%M:%S", tz="UTC")
    records <- range(WSdata$datetime)
    if(!(sat.datetime > records[1] && sat.datetime < records[2])){
      warning("The satellite overpass (", format(sat.datetime, usetz = TRUE),
              ") is outside the weather station records (",
              format(records[1], usetz = TRUE), " to ",
              format(records[2], usetz = TRUE), "), so the conditions at the ",
              "overpass (at.sat) were not computed")
      class(result) <- "waterWeatherStation"
      return(result)
    }
    WS.prev<-suppressWarnings(WSdata[WSdata$datetime == tail(datetime[datetime < 
                                                       sat.datetime],1),])
    WS.after <- suppressWarnings(WSdata[WSdata$datetime == datetime[datetime > 
                                                     sat.datetime][1],])
    delta1 <- as.numeric(difftime(WS.after$datetime, 
                                  WS.prev$datetime, units="secs"))
    delta2 <- as.numeric(difftime(sat.datetime, WS.prev$datetime, units="secs"))
    sat.data <- WS.prev + (WS.after - WS.prev)/delta1 * delta2
    sat.data[,2:6] <- round(sat.data[,2:6],2)
    result$at.sat <- sat.data
  }
  class(result) <- "waterWeatherStation"
  return(result)
}


#' Prepares weather station data 2
#' @param WSdata             csv file with weather station data
#' @param ...                additional parameter to pass to read.csv()
#' @param height             weather station sensors height in meters
#' @param lat                latitude of weather station in decimal degrees. 
#' Negative values for south latitude
#' @param long               longitude of weather station in decimal degrees. 
#' Negative values for west longitude
#' @param elev               elevation of weather station in meters
#' @param columns            columns order of needed data. Vector containing 
#' "date", "time", "radiation", "wind", "RH", "temp" and "rain". Other values are 
#' ignored. If you have a column with date and time in the same column, you can
#' include "datetime" and "date" and "time" are no longer needed.
#' @param date.format        date format. See strptime format argument.
#' @param time.format        time format. See strptime format argument.
#' @param datetime.format    datetime format. See strptime format argument.
#' @param tz                 timezone of the weather station dates. If not 
#' present assumes the same timezone as the computer running the code. See 
#' strptime for details. 
#' @param cf                 conversion factor to convert radiation, wind, and 
#' temperature to W/m2; m/s and Celsius. See Details.
#' @param MTL                Metadata file. If not provided will look for one on
#' working directory. If provided or present will calculate weather conditions
#' on satellite flyby.
#' @details 
#' For cf, if your data is in W/m2, km/h and Celsius (radiation, wind, 
#' temperature), cf should be: cf = c(1,0.2777778,1)
#' @family Weather station related functions
#' @author Guillermo Federico Olmedo
#' @export
read.WSdata2 <- function(WSdata, ..., height = 2.2, lat, long, elev,
                         columns=c("date" = 1, "time" = 2, "radiation" = 3,
                                   "wind" = 4, "RH" = 6, "temp" = 7, "rain" = 8),
                         date.format = "%d/%m/%Y", time.format = "%H:%M:%S", 
                         datetime.format = "%Y-%m-%d %H:%M:%S", tz = "",
                         cf = c(1, 1, 1), MTL){
  if(is.character(WSdata)){
    if(!file.exists(WSdata)){stop("Weather station file not found: ", WSdata)}
    WSdata <- utils::read.csv2(WSdata, ...)
  }
  read.WSdata(WSdata, height = height, lat = lat, long = long, elev = elev,
              columns = columns, date.format = date.format,
              time.format = time.format, datetime.format = datetime.format,
              tz = tz, cf = cf, MTL = MTL)
}


#' Plot method for waterWeatherStation S3 class
#' @param x       waterWeatherStation object. See read.WSdata()
#' @param hourly  If TRUE will plot only hourly data, instead of all data 
#' from the waterWeatherStation object
#' @param sat     If TRUE, and if the waterWeatherStation object was created
#' using a Landsat Metadata File, will plot only data from day of the satellite
#' overpass 
#' @param date    When sat = FALSE, the date to plot in a multi-day object
#' @param ...      additional parameters to pass to plot()
#' @author Guillermo Federico Olmedo
#' @importFrom utils read.csv
#' @importFrom graphics abline axis axis.POSIXct mtext par
#' @export
#' @family Weather station related functions
#' @method plot waterWeatherStation
plot.waterWeatherStation <- function(x, hourly=FALSE, sat=TRUE, date, ...){
  # Based on http://evolvingspaces.blogspot.cl/2011/05/multiple-y-axis-in-r-plot.html
  WSp <- x$hourly
  if(hourly == FALSE) {WSp <- x$alldata}
  long=FALSE
  daily_mean = ""
  if(sat == TRUE & "datetime" %in% names(x$at.sat)){
    atsat  <- as.POSIXlt(x$at.sat$datetime)
    WSp <- WSp[as.POSIXlt(WSp$datetime)$yday == atsat$yday 
               & as.POSIXlt(WSp$datetime)$year == atsat$year,]
  } else {
    if(nrow(x$daily)>5){WSp <- x$daily
    names(WSp) <- c("datetime", "radiation", "wind", "RH", "temp", 
                    "temp_max", "temp_min", "ea", "rain")
    long=TRUE
    daily_mean = "(mean daily value)"}
  }
  
  if(!missing(date)){
    date <- as.POSIXlt(date)
    WSp <- WSp[as.POSIXlt(WSp$datetime)$yday == date$yday 
               & as.POSIXlt(WSp$datetime)$year == date$year,]
  }
  WSp <- WSp[order(WSp$datetime),]
  time <- WSp$datetime
  graphics::par(mar=c(5, 7, 4, 9.5) + 0.1)
  plot(time, WSp$radiation, axes=F, ylim=c(0,max(WSp$radiation)), xlab="", 
       ylab="",type="l",col="red", xlim=range(time), ...)
  if(sat == TRUE & "datetime" %in% names(x$at.sat)){
    graphics::abline(v=as.POSIXct(atsat), lwd=5, col="gray")
    graphics::text(atsat, max(WSp$radiation)*0.85, "satellite overpass", cex=0.7, 
                   adj=c(NA, -0.5), srt=90, col="gray")
  }
  graphics::points(time,WSp$radiation,pch=20,col="red")
  graphics::axis(2, ylim=c(0,max(WSp$radiation)),col="red",lwd=1, cex.axis=0.5, tcl=-0.25)
  graphics::mtext(2,text=paste("Solar radiation (W.m-2)", daily_mean),line=1.7, cex=0.7)
  # Rain
  #WSp[WSp$rain == 0] <- NA
  if(!all(is.na(WSp$rain))){ 
    graphics::par(new=T)
    plot(time, WSp$rain, axes=F, ylim=c(0,max(WSp$rain)*1.7), xlab="", ylab="", 
         type="h",col="light blue",lty=1, main="",xlim=range(time),lwd=3)
    graphics::axis(4, ylim=c(0,max(WSp$rain)*1.7),col="light blue", line = 6.7, lwd=1, cex.axis=0.5, tcl=-0.25)
    graphics::mtext(4,text=paste("Rain (mm)", daily_mean), cex=0.7, line = 8.4)
    }
  #
  graphics::par(new=T)
  plot(time, WSp$wind, axes=F, ylim=c(0,max(WSp$wind)), xlab="", ylab="", 
       type="l",col="green",lty=2, main="",xlim=range(time),lwd=2)
  graphics::axis(2, ylim=c(0,max(WSp$wind)),col="green",lwd=1,line=3.5, cex.axis=0.5, tcl=-0.25)
  graphics::points(time, WSp$wind,pch=20,col="green")
  graphics::mtext(2,text=paste("Wind speed (m.s-1)", daily_mean),line=5.5, cex=0.7)
  graphics::par(new=T)
  plot(time, WSp$temp, axes=F, ylim=c(0,max(WSp$temp)), xlab="", ylab="", 
       type="l",col="black",lty=2, main="",xlim=range(time),lwd=2)
  graphics::axis(4, ylim=c(0,max(WSp$temp)),col="black",lwd=1, cex.axis=0.5, tcl=-0.25)
  graphics::points(time, WSp$temp,pch=20,col="black")
  graphics::mtext(4,text=paste("Temperature (C)", daily_mean),line=1.7, cex=0.7)
  graphics::par(new=T)
  plot(time, WSp$ea, axes=F, ylim=c(0,max(WSp$ea)), xlab="", ylab="", 
       type="l",col="blue",lty=2, main="",xlim=range(time),lwd=2)
  graphics::axis(4, ylim=c(0,max(WSp$ea)),col="blue",lwd=1,line=3.5, cex.axis=0.5, tcl=-0.25)
  graphics::points(time, WSp$ea,pch=20,col="blue")
  graphics::mtext(4,text=paste("vapor pressure (kPa)", daily_mean),line=5.5, cex=0.7)
  if(long==FALSE){r <- as.POSIXct(round(range(time), "hours"))
  graphics::axis.POSIXct(1, at = seq(r[1], r[2], by = "hour"), format = "%H:%M")}
  if(long==TRUE){r <- range(time)
  graphics::axis.Date(1, at = seq(r[1], r[2], by = "day"))}
  graphics::mtext("Time",side=1,col="black",line=2)
}

#' Print method for waterWeatherStation S3 class
#' @param x        waterWeatherStation object. See read.WSdata()
#' @param ...      additional parameters to pass to print()    
#' @author Guillermo Federico Olmedo
#' @export
#' @family Weather station related functions
#' @method print waterWeatherStation
print.waterWeatherStation <- function(x, ...){
  cat("Weather Station @ lat:", round(x$location$lat, 2), "long:", 
      round(x$location$long, 2), "elev:", round(x$location$elev, 2), "\n")
  if(!is.null(x$at.sat)){
  cat("Summary:\n")
  print(summary(x$alldata[format(x$alldata$datetime, format = "%d-%m-%Y") == 
                            format(x$at.sat$datetime, format = "%d-%m-%Y"), 2:7]))
  cat("\n Conditions at satellite overpass:\n")
  print(x$at.sat)} else {
    cat("Summary:\n")
    print(summary(x$alldata[,1:7]))
  }
}


#' Export data.frame from waterWeatherStation Object
#' @description 
#' Export weather conditions at satellite flyby from waterWeatherStation object
#' @param WeatherStation    waterWeatherStation object. See read.WSdata()
#' @author Guillermo Federico Olmedo
#' @export
getDataWS <- function(WeatherStation){
  if(is.null(WeatherStation$at.sat)){
    stop("The weather station data has no conditions at the satellite overpass ",
         "(at.sat): use read.WSdata(..., MTL = <Landsat MTL file>) with records ",
         "that include the overpass")
  }
  date <- as.POSIXlt(WeatherStation$at.sat$datetime, format="%Y-%m-%d %H:%M:%S")
  WeatherStation2 <- data.frame(wind=WeatherStation$at.sat$wind, 
                               RH=WeatherStation$at.sat$RH, 
                               temp=WeatherStation$at.sat$temp, 
                               radiation=WeatherStation$at.sat$radiation, 
                               height=WeatherStation$location$height,
                               lat=WeatherStation$location$lat, 
                               long=WeatherStation$location$long, 
                               elev=WeatherStation$location$elev, 
                               DOY=date$yday+1, hours=date$hour + date$min/60 + 
                                 date$sec/3600)
  return(WeatherStation2)
}


