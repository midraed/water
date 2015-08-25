H <- function(anchors, Ts, LAI, albedo, Z.om, n=1, anchors.method= "random", 
              WeatherStation, ETo.hourly, ETp.coef= 1.05, Z.om.ws=0.0018, 
              mountainous=FALSE, DEM, Rn, G, plots=TRUE){
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
      } else cold <- sample(which(values(LAI>4)),n)}
    if(anchors.method=="random2"){
      if(!max(values(Ts))>=310){
        warning(paste("Ts max value is", round(max(values(Ts)),2), 
                      "and the expected was Ts >= 310 for hot pixel"))
        hot <- sample(which(values(Ts)>quantile(Ts,0.9999)& 
                              values(albedo>=0.13) & values(albedo<=0.15) &
                              values(Z.om>0) & values(Z.om<=0.005)),n)  # CITRA
      } else hot <- sample(which(values(Ts>310) & 
                                   values(albedo>=0.13) & values(albedo<=0.15) &
                                   values(Z.om>0) & values(Z.om<=0.005)),n) # CITRA
      if(!max(values(LAI))>=3){
        warning(paste("LAI max value is", round(max(values(LAI)),2), 
                      "and the expected was LAI >= 3 for cold pixel"))
        cold <- sample(which(values(LAI)>quantile(LAI,0.9999)& values(LAI<=6) & values(Ts>288) & 
                               values(albedo>=0.18) & values(albedo<=0.25) &
                               values(Z.om>=0.08)),n)  # CITRA
      } else cold <- sample(which(values(LAI>=3) & values(LAI<=6) & values(Ts>288) & 
                                    values(albedo>=0.18) & values(albedo<=0.25) &
                                    values(Z.om>=0.08)),n)} # CITRA 
    }
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
  u.ws <- WeatherStation$wind * 0.41 / log(WeatherStation$height/Z.om.ws)
  u200 <- u.ws / 0.41 * log(200/Z.om.ws)
  if(mountainous==TRUE){
    u200 <- u200 * (1+0.1*((DEM-WeatherStation$elev)/1000))
  }
  friction.velocity <- 0.41 * u200 / log(200/Z.om) 
  r.ah <- log(2/0.1)/(friction.velocity*0.41) #ok
  ### Iteractive process start here:
  #delta.r.ah <- vector()
  #delta.H <- vector()
  ETr.hourly <-  0.56 ## need to calculate this first - MCB uses PM.ASCE hourly eq
  LE.cold <- ETr.hourly * ETp.coef * (2.501 - 0.002361*(Ts[cold]-273.15))*(1e6)/3600 
  # here uses latent.heat.vapo
  H.cold <- Rn[cold] - G[cold] - LE.cold #ok
  result <- list()
  for(i in 1:5){
    print(paste("iteraction #", i))
    ### We calculate dT and H 
    dT.cold <- H.cold * r.ah[cold] / (air.density[cold]*1004)
    dT.hot <- (Rn[hot] - G[hot]) * r.ah[hot] / (air.density[hot]*1004)
    a <- (dT.hot - dT.cold) / (Ts.datum[hot] - Ts.datum[cold])
    b <- -a * Ts.datum[cold] + dT.cold
    print(a)
    print(b)
    dT <- as.numeric(a) * Ts.datum + as.numeric(b)   #ok
    rho <- 349.467*((((Ts-dT)-0.0065*DEM)/(Ts-dT))^5.26)/Ts  
    H <- rho * 1004 * dT / r.ah
    Monin.Obukhov.L <- (air.density * -1004 * friction.velocity^3 * Ts) / 
      (0.41 * 9.807 * H)
    ### Then we calculate L and phi200, phi2, and phi0.1 
    ## !!! This is very time consumig... maybe only for hot and cold pixels?
    phi.200 <- raster(Monin.Obukhov.L) 
    # copy raster extent and pixel size, not values!
    phi.2 <- raster(Monin.Obukhov.L)
    phi.01 <- raster(Monin.Obukhov.L)
    ## stable condition = L > 0
    phi.200[Monin.Obukhov.L > 0] <- -5*(2/Monin.Obukhov.L)[Monin.Obukhov.L > 0] #ok
    phi.2[Monin.Obukhov.L > 0] <- -5*(2/Monin.Obukhov.L)[Monin.Obukhov.L > 0]  #ok
    phi.01[Monin.Obukhov.L > 0] <-  -5*(0.1/Monin.Obukhov.L)[Monin.Obukhov.L > 0] #ok
    ## unstable condition = L < 0
    x.200 <- (1- 16*(200/Monin.Obukhov.L))^0.25 #ok
    x.2 <- (1- 16*(2/Monin.Obukhov.L))^0.25 #ok
    x.01 <- (1- 16*(0.1/Monin.Obukhov.L))^0.25 # ok
    phi.200[Monin.Obukhov.L < 0] <- (2 * log((1+x.200)/2) + log((1 + x.200^2) /2) - 
                                       2* atan(x.200) + 0.5 * pi)[Monin.Obukhov.L < 0] #ok
    phi.2[Monin.Obukhov.L < 0] <- (2 * log((1 + x.2^2) / 2))[Monin.Obukhov.L < 0]
    phi.01[Monin.Obukhov.L < 0] <- (2 * log((1 + x.01^2) / 2))[Monin.Obukhov.L < 0]
    print(r.ah[cold])
    print(dT[cold])
    print(r.ah[hot])
    print(dT[hot])
    ## And finally, r.ah and friction velocity
    friction.velocity <- 0.41 * u200 / (log(200/Z.om) - phi.200)
    r.ah <- (log(2/0.1) - phi.2 + phi.01) / (friction.velocity * 0.41) # ok ok
    #r.ah[r.ah > quantile(r.ah, 0.9999)] <-  NA # to eliminate very high alues
    #delta.r.ah <- c(delta.r.ah, mean(values(prev.r.ah), na.rm=T) - mean(values(r.ah), na.rm=T))
    #delta.H
    #prev.r.ah <- r.ah
    #prev.H <- H
  } # End interactive process
  dT <- save.load.clean(imagestack = dT, file = paste0(result.folder,"dT.tif"), overwrite=TRUE)
  H <- save.load.clean(imagestack = H, file = paste0(result.folder,"H.tif"), overwrite=TRUE)
  result$a <- a
  result$b <- b
  result$dT <- dT
  result$H <- H
  result$hot.and.cold <- data.frame(pixel=integer(),  X=integer(), Y=integer(), Ts=double(), 
                                    LAI=double(), type=factor(levels = c("hot", "cold")))
  for(i in 1:n){result$hot.and.cold[i, ] <- c(hot[i], xyFromCell(LAI, hot[i]),
                                              Ts[hot][i], round(LAI[hot][i],2), "hot")}
  for(i in 1:n){result$hot.and.cold[i+n, ] <- c(cold[i], xyFromCell(LAI, hot[i]), 
                                                Ts[cold][i], round(LAI[cold][i],2), "cold")}
  return(result)
}




