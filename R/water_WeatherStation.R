#' Prepares weather station data
#' @examples 
#' csvfile <- system.file("extdata", "apples.csv", package="water")
#' MTLfile <- system.file("extdata", "L7.MTL.txt", package="water")
#' WS <- read.WSdata(WSdata = csvfile, date.format = "%d/%m/%Y", 
#'                   lat=-35.42222, long= -71.38639, elev=201, height= 2.2,
#'                   MTL = MTLfile)
#' print(WS)
#' plot(WS, alldata=FALSE)
#' plot(WS, alldata=TRUE)
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' Landsat 7 Metadata example file available from the U.S. Geological Survey.
#' Weather Station example file courtesy of CITRA, Universidad de Talca, Chile
#' @export
read.WSdata <- function(WSdata, ..., height = 2.2, lat, long, elev,
                        columns = c("date", "time", "radiation", "wind", NA,
                                    "RH", "temp", "pp"),
                        date.format = "%Y-%m-%d", time.format = "%H:%M:%S", 
                        datetime.format = "%Y-%m-%d %H:%M:%S", tz = "",
                        cf = c(1, 1, 1, 1), MTL){
  result <- list()
  result$location <- data.frame(lat=lat, long=long, elev=elev, height=height)
  WSdata <- read.csv(WSdata, ...)
  if("date" %in% columns & "time" %in% columns){
    datetime  <- paste(WSdata[, which(columns == "date")], 
                       WSdata[, which(columns == "time")])
    datetime <- strptime(datetime, format = paste(date.format, time.format), 
                         tz = tz)
  } else {
    if("datetime" %in% columns){
      datetime  <- WSdata[, which(columns == "datetime")]
      datetime <- strptime(datetime, format = datetime.format, tz = tz)
    } else {message("ERROR: date and time or datetime are needed")}
  }
  radiation = WSdata[, which(columns == "radiation")] * cf[1]
  wind =  WSdata[, which(columns == "wind")] * cf[2]
  RH =  WSdata[, which(columns == "RH")] * cf[3]
  temp =  WSdata[, which(columns == "temp")] * cf[4]
  ea = (RH/100)*0.6108*exp((17.27*temp)/(temp+237.3))
  WSdata <- data.frame(datetime=datetime, radiation=radiation, wind=wind,
                       RH=RH, ea=ea, temp=temp)
  result$alldata <- WSdata
  result$hourly <- WSdata[datetime$min==0,] 
  ## Join with satellite data
  if(missing(MTL)){MTL <- list.files(pattern = "MTL.txt", full.names = T)}
  if(length(MTL)!=0){
    MTL <- readLines(MTL, warn=FALSE)
    time.line <- grep("SCENE_CENTER_TIME",MTL,value=TRUE)
    date.line <- grep("DATE_ACQUIRED",MTL,value=TRUE)
    sat.time <-regmatches(time.line,regexec(text=time.line,
           pattern="([0-9]{2})(:)([0-9]{2})(:)([0-9]{2})(.)([0-9]{2})"))[[1]][1]
    sat.date <-regmatches(date.line,regexec(text=date.line,
                        pattern="([0-9]{4})(-)([0-9]{2})(-)([0-9]{2})"))[[1]][1]
    sat.datetime <- strptime(paste(sat.date, sat.time), 
                             format = "%Y-%m-%d %H:%M:%S", tz="GMT")
    WS.prev<-WSdata[WSdata$datetime == tail(datetime[datetime < 
                                                       sat.datetime],1),]
    WS.after <- WSdata[WSdata$datetime == datetime[datetime > 
                                                     sat.datetime][1],]
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
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @export
read.WSdata2 <- function(WSdata, ..., height = 2.2, lat, long, elev,
                         columns = c("date", "time", "radiation", "wind", NA,
                                     "RH", "temp", "pp"),
                         date.format = "%d/%m/%Y", time.format = "%H:%M:%S", 
                         datetime.format = "%Y-%m-%d %H:%M:%S", tz = "",
                         cf = c(1, 3.6, 1, 1), MTL){
  read.WSdata(WSdata, ..., height = 2.2, lat = lat, long = long, elev = elev,
              columns = columns, date.format = date.format, 
              time.format = time.format, datetime.format = datetime.format,
              tz = tz, cf = cf, MTL)
}


#' Plot method for waterWeatherStation S3 class
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @export
#' @method plot waterWeatherStation
plot.waterWeatherStation <- function(x, alldata=TRUE, ...){
  # Based on http://evolvingspaces.blogspot.cl/2011/05/multiple-y-axis-in-r-plot.html
  WSp <- x$hourly
  atsat  <- as.POSIXct(x$at.sat$datetime)
  if(alldata == TRUE) {WSp <- x$alldata}
  time <- WSp$datetime
  par(mar=c(5, 7, 4, 7) + 0.1)
  plot(time, WSp$radiation, axes=F, ylim=c(0,max(WSp$radiation)), xlab="", 
       ylab="",type="l",col="red", main="",xlim=range(time), ...)
  abline(v=atsat, lwd=5, col="gray")
  graphics::text(atsat, max(WSp$radiation)*0.85, "Satellite Flyby", cex=0.7, 
                 adj=c(NA, -0.5), srt=90, col="gray")
  points(time,WSp$radiation,pch=20,col="red")
  axis(2, ylim=c(0,max(WSp$radiation)),col="red",lwd=1, cex.axis=0.5)
  mtext(2,text="Solar radiation (W.m-2)",line=1.7, cex=0.7)
  par(new=T)
  plot(time, WSp$wind, axes=F, ylim=c(0,max(WSp$wind)), xlab="", ylab="", 
       type="l",col="green",lty=2, main="",xlim=range(time),lwd=2,...)
  axis(2, ylim=c(0,max(WSp$wind)),col="green",lwd=1,line=3.5, cex.axis=0.5)
  points(time, WSp$wind,pch=20,col="green")
  mtext(2,text="Wind speed (m.s-1)",line=5.5, cex=0.7)
  par(new=T)
  plot(time, WSp$temp, axes=F, ylim=c(0,max(WSp$temp)), xlab="", ylab="", 
       type="l",col="black",lty=2, main="",xlim=range(time),lwd=2,...)
  axis(4, ylim=c(0,max(WSp$temp)),col="black",lwd=1, cex.axis=0.5)
  points(time, WSp$temp,pch=20,col="black")
  mtext(4,text="Temperature (C)",line=1.7, cex=0.7)
  par(new=T)
  plot(time, WSp$ea, axes=F, ylim=c(0,max(WSp$ea)), xlab="", ylab="", 
       type="l",col="blue",lty=2, main="",xlim=range(time),lwd=2, ...)
  axis(4, ylim=c(0,max(WSp$ea)),col="blue",lwd=1,line=3.5, cex.axis=0.5)
  points(time, WSp$ea,pch=20,col="blue")
  mtext(4,text="vapor pressure (kPa)",line=5.5, cex=0.7)
  axis.POSIXct(1, range(time))
  mtext("Time",side=1,col="black",line=2)
}

#' Print method for waterWeatherStation S3 class
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @export
#' @method print waterWeatherStation
print.waterWeatherStation <- function(x, ...){
  cat("Weather Station at lat:", round(x$location$lat, 2), "long:", 
      round(x$location$lat, 2), "elev:", round(x$location$elev, 2), "\n")
  cat("Summary:\n")
  print(summary(x$alldata[,2:6]), ...)
  if(!is.null(x$at.sat)){
  cat("\n Conditions at satellite flyby:\n")
  print(x$at.sat)}
}


#' Export data.frame from waterWeatherStation Object
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @export
getDataWS <- function(WeatherStation){
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


