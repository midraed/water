aerodynamic.transport.i <- function(anchors, Ts, LAI, n=1, anchors.method= "random",
                                    wind, ETr, ETp.coef= 1.05, Z.om.ws=0.0018, height.ws=2, 
                                    mountainous=FALSE, elev.ws, DEM,
                                    Rn, G, plots=TRUE){
  start <- Sys.time ()
  ### Some values used later
  Ts.datum <- Ts - (DEM - 702) * 6.49 / 1000
  P <- 101.3*((293-0.0065 * DEM)/293)^5.26
  air.density <- 1000 * P / (1.01*(Ts)*287)
  latent.heat.vaporization <- (2.501-0.00236*(Ts-273.15))# En el paper dice por 1e6
  ### We create anchors points if they dont exist
  if(missing(anchors)){
    if(anchors.method=="random"){
      if(!max(values(Ts))>=310){
        warning(paste("Ts max value is", round(max(values(Ts)),2), 
                      "and the expected was Ts >= 310 for hot pixel"))
        hot <- sample(which(values(Ts)>quantile(Ts,0.9999)),n)  
      } else hot <- sample(which(values(Ts>310)),n) 
      if(!max(values(LAI))>=4){
        warning(paste("LAI max value is", round(max(values(LAI)),2), 
                      "and the expected was LAI >= 4 for cold pixel"))
        cold <- sample(which(values(LAI)>quantile(LAI,0.9999)),n)  
      } else cold <- sample(which(values(LAI>4)),n)  
    }}
  else{
    hot <- as.numeric(hc.data$pixel[hc.data$type =="hot"])
    cold <- as.numeric(hc.data$pixel[hc.data$type =="cold"])
  }
  ### We can plot anchor points
  if(plots==TRUE){
    plot(LAI, main="LAI and hot and cold pixels")
    points(xyFromCell(LAI, hot), col="red", pch=3)
    points(xyFromCell(LAI, cold), col="blue", pch=4)
  }
  ### We calculate the initial conditions assuming neutral stability
  friction.velocity <- 0.41 * wind / log(height.ws/Z.om.ws) # Tasumi 2003
  u200 <- wind * (log(200/Z.om)/0.41) # Tasumi 2003 / not the same as Allen 2007
  if(mountainous==TRUE){
    u200 <- u200 * (1+0.1*((DEM-elev.ws)/1000))
  }
  friction.velocity <- 0.41 * u200 / log(200/Z.om)
  r.ah <- log(2/0.1)/friction.velocity*0.41
  ### Iteractive process start here:
  delta.r.ah <- vector()
  delta.fric <- vector()
  for(i in 1:15){
    prev.r.ah <- r.ah
    print(paste("iteraction #", i))
    ### We calculate dT and H 
    dT.hot <- (Rn[hot] - G[hot])*r.ah[hot]/(air.density[hot]*1007)
    LE.cold <- ETp.coef * ETr * latent.heat.vaporization[cold]# instead of ETr better use ETr~LAI
    H.cold <- Rn[cold]-G[cold]-LE.cold
    dT.cold <- H.cold*r.ah[cold]/(air.density[cold]*1007)
    a <- mean((dT.hot-dT.cold)/(Ts.datum[hot]-Ts.datum[cold]), na.rm=T)
    b <- mean((dT.hot-a)/Ts.datum[hot], na.rm=T)
    dT <- a+b*Ts.datum
    H <- air.density*1007*(dT/r.ah)
    print(a)
    print(b)
    print(r.ah[cold])
    print(dT[cold])
    print(r.ah[hot])
    print(dT[hot])
    ### Then we calculate L and phi200, phi2, and phi0.1
    Monin.Obukhov.L <- (air.density * 1007 * friction.velocity^3 * Ts) / (0.41 * 9.807 * H)
    phi.200 <- raster(Monin.Obukhov.L) # copy raster extent and pixel size, not values!
    phi.2 <- raster(Monin.Obukhov.L)
    phi.01 <- raster(Monin.Obukhov.L)
    ## stable condition = L > 0
    phi.200[Monin.Obukhov.L > 0] <- -5*(2/Monin.Obukhov.L)[Monin.Obukhov.L > 0]
    phi.2[Monin.Obukhov.L > 0] <- -5*(2/Monin.Obukhov.L)[Monin.Obukhov.L > 0]
    phi.01[Monin.Obukhov.L > 0] <-  -5*(0.1/Monin.Obukhov.L)[Monin.Obukhov.L > 0]
    ## unstable condition = L < 0
    x.200 <- (1- 16*(200/Monin.Obukhov.L))^0.25
    x.2 <- (1- 16*(2/Monin.Obukhov.L))^0.25
    x.01 <- (1- 16*(0.1/Monin.Obukhov.L))^0.25
    phi.200[Monin.Obukhov.L < 0] <- (2 * log(1+x.200/2) + log(1 + x.200^2 /2) - 2* atan(x.200) + 0.5 * pi)[Monin.Obukhov.L < 0]
    phi.2[Monin.Obukhov.L < 0] <- (2 * log(1 + x.2^2 / 2))[Monin.Obukhov.L < 0]
    phi.01[Monin.Obukhov.L < 0] <- (2 * log(1 + x.01^2 / 2))[Monin.Obukhov.L < 0]
    ## And finally, r.ah and friction velocity
    friction.velocity <- 0.41 * u200 / (log(200/Z.om) - phi.200)
    ### SAVE
    r.ah <- (log(2/0.1) - phi.2 + phi.01) / friction.velocity * 0.41
    r.ah[r.ah > quantile(r.ah, 0.999)] <-  NA # to eliminate very high values
    delta.r.ah <- c(delta.r.ah, mean(values(prev.r.ah), na.rm=T) - mean(values(r.ah), na.rm=T))
  } # End interactive process
  print (Sys.time () - start)
  return(delta.r.ah)
}



# 
#   # Prepare and populate result list
#   result <- list()
#   result$dT <- dT
#   result$hot.and.cold <- data.frame(pixel=integer(),  X=integer(), Y=integer(), Ts=double(), 
#                                     LAI=double(), type=factor(levels = c("hot", "cold")))
#   for(i in 1:n){result$hot.and.cold[i, ] <- c(hot[i], xyFromCell(LAI, hot[i]), Ts[hot][i], round(LAI[hot][i],2), "hot")}
#   for(i in 1:n){result$hot.and.cold[i+n, ] <- c(cold[i], xyFromCell(LAI, hot[i]), Ts[cold][i], round(LAI[cold][i],2), "cold")}
#   removeTmpFiles(h=0)
# 
#   return(result)
