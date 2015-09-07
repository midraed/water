#' Calculates ET using Penman Monteith hourly formula
#' @param WeatherStation a data frame with all the needed fields (see example)
#' @param hours time of the day in hours in 24hs format
#' @param DOY day of year
#' @param long.z longitude for local time
#' @return ET hourly in mm.h-1
#' @author Guillermo F Olmedo
#' @examples 
#' WeatherStation  <- data.frame(wind=4.72,
#'                               RH=59, 
#'                               Ta=24.3,
#'                               Gr.Rad=675, 
#'                               height=2.2, 
#'                               lat=-35.37, 
#'                               long=71.5946, 
#'                               elev=124)
#'   hourlyET(WeatherStation, hours=10.5, DOY=363, long.z=71.635)
#' @export
#' @references 
#' Allen 2005 ASCE
hourlyET <- function(WeatherStation, hours, DOY, long.z=WeatherStation$long){
  TaK <- WeatherStation$Ta + 273.16
  Rs <- WeatherStation$Gr.Rad * 3600 / 1e6
  P <- 101.3*((293-0.0065*WeatherStation$elev)/293)^5.26
  psi <- 0.000665*P
  Delta <- 2503 * exp((17.27*WeatherStation$Ta)/
                        (WeatherStation$Ta+237.3))/((WeatherStation$Ta+237.3)^2)
  ea.sat <- 0.6108*exp((17.27*WeatherStation$Ta)/(WeatherStation$Ta+237.3))
  ea <- (WeatherStation$RH/100)*ea.sat
  DPV <- ea.sat - ea
  dr <- 1 + 0.033*(cos(2*pi*DOY/365))
  delta <- 0.409*sin((2*pi*DOY/365)-1.39)
  phi <- WeatherStation$lat*(pi/180)
  b <- 2*pi*(DOY-81)/364
  Sc <- 0.1645*sin(2*b)-0.1255*cos(b)-0.025*sin(b)
  hour.angle <- (pi/12)*((hours+0.06667*(WeatherStation$long*pi/180-long.z*pi/180)+Sc)-12)
  w1 <- hour.angle-((pi)/24)
  w2 <- hour.angle+((pi)/24)
  hour.angle.s <- acos(-tan(phi)*tan(delta))
  w1c <- ifelse(w1< -hour.angle.s, -hour.angle.s, 
                ifelse(w1>hour.angle.s, hour.angle.s, ifelse(w1>w2, w2, w1)))
  w2c <- ifelse(w2< -hour.angle.s, -hour.angle.s, 
                ifelse(w2>hour.angle.s, hour.angle.s, w2))
  Beta <- asin((sin(phi)*sin(delta)+cos(phi)*cos(delta)*cos(hour.angle)))
  Ra <- ifelse(Beta <= 0, 1e-45, ((12/pi)*4.92*dr)*
                 (((w2c-w1c)*sin(phi)*sin(delta))+(cos(phi)*cos(delta)*(sin(w2)-sin(w1)))))
  Rso <- (0.75+2e-5*WeatherStation$elev)*Ra
  Rs.Rso <- ifelse(Rs/Rso<=0.3, 0, ifelse(Rs/Rso>=1, 1, Rs/Rso))
  fcd <- ifelse(1.35*Rs.Rso-0.35<=0.05, 0.05, 
                ifelse(1.35*Rs.Rso-0.35<1, 1.35*Rs.Rso-0.35,1))
  Rn.a <- ((1-0.23)*Rs) - (2.042e-10*fcd*(0.34-0.14*(ea^0.5))*TaK^4)
  G.day <- Rn.a * 0.1
  wind.2 <- WeatherStation$wind *(4.87/(log(67.8*WeatherStation$height-5.42)))
  ETo.hourly <- ((0.408*Delta*(Rn.a-G.day))+(psi*(37/TaK)*wind.2*(DPV)))/
    (Delta+(psi*(1+(0.24*wind.2))))
  return(ETo.hourly)
}

#' Calculates ET-24hs from energy balance and Weather Station 
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @export
ET24h <- function(Rn, G, H, Ts, WeatherStation, ETr.daily, C.rad=1){
  LE = Rn - G - H
  ET.inst <- 3600*LE/((2.501 - 0.00236 * (Ts - 273.15)) * (1e6))
  ETo.hourly <- hourlyET(WeatherStation, WeatherStation$hours, WeatherStation$DOY)
  ETr.Fr <- ET.inst/ETo.hourly
  ET.24 <- ETr.Fr * ETr.daily * C.rad
  rgb.palette <- colorRampPalette(c("red3","snow2","blue"),  space = "rgb")
  print(spplot(ET.24, col.regions=rgb.palette, main= "24-Hour Evapotranspiration (mm/day)",
               colorkey=list(height=1), at=seq(0,ceiling(ETr.daily*1.5),length.out=50), maxpixels=ncell(ET.24) * 0.3))
  saveLoadClean(imagestack = ET.24, 
                file = "ET24", overwrite=TRUE)
}


#' Calculates daily ET using Penman Monteith hourly formula for every hour
#' @param WeatherStation a data frame with all the needed fields (see example)
#' @param DOY day of year
#' @param long.z longitude for local time
#' @return ET daily in mm.h-1
#' @author Guillermo F Olmedo
#' @export
#' @references 
#' Allen 2005 ASCE
dailyEToPM <- function(WeatherStation, DOY, long.z=WeatherStation$long){
  print("not yet")
  return()
}
