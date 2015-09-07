#' Iterative function to estimate H and R.ah - CITRA method
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @references 
#' CITRA y MCB (com pers)
#' @export
sensibleHeatFlux.CITRA <- function(anchors, image, Ts, LAI, albedo, Z.om, n=1,
                                   anchors.method= "CITRA-MCB", WeatherStation, 
                                   ETp.coef= 1.05, Z.om.ws=0.0018, sat="auto", 
                                   ESPA=F, mountainous=FALSE, DEM, Rn, G, 
                                   plots=TRUE, deltaTemp=5) {
  ### Some values used later
  if(sat=="auto"){sat <- getSat(getwd())}
  if(sat=="L8" & ESPA==T){
    removeTmpFiles(h=0)
    sr.red <- raster(list.files(path = path, 
                                pattern = "_sr_band4.tif", full.names = T))
    sr.nir <- raster(list.files(path = path, 
                                pattern = "_sr_band5.tif", full.names = T))
    sr.4.5 <- stack(sr.red, sr.nir)
    sr.4.5 <- aoiCrop(sr.4.5, aoi)
  }
  if(sat=="L7"){
    sr.4.5 <- stack(image[[3]], image[[4]])
  }
  NDVI <- (sr.4.5[[2]] - sr.4.5[[1]])/(sr.4.5[[1]] + sr.4.5[[2]])
  ETo.hourly <- hourlyET(WeatherStation, WeatherStation$hours, WeatherStation$DOY)
  Ts.datum <- Ts - (DEM - WeatherStation$elev) * 6.49 / 1000
  P <- 101.3*((293-0.0065 * DEM)/293)^5.26
  air.density <- 1000 * P / (1.01*(Ts)*287)
  latent.heat.vaporization <- (2.501-0.00236*(Ts-273.15))# En el paper dice por 1e6

  
### We create anchors points if they dont exist---------------------------------
  if(missing(anchors)){
    if(anchors.method=="CITRA-MCB"){
      minT <- quantile(Ts[LAI>=3&LAI<=6&albedo>=0.18&albedo<=0.25&Z.om>=0.03&Z.om<=0.08], 0.25)
      if(minT+deltaTemp<288){minT = 288 + deltaTemp}
      cold <- sample(which(values(LAI>=3) & values(LAI<=6) &  
                             values(albedo>=0.18) & values(albedo<=0.25) &
                             values(Z.om>=0.03) & values(Z.om<=0.08) &
                             values(Ts<(minT+deltaTemp)) & values(Ts>288)),n)
      maxT <- max(Ts[albedo>=0.13&albedo<=0.15&NDVI>=0.1&NDVI<=0.28&Z.om<=0.005])
      hot <- sample(which(values(albedo>=0.13) & values(albedo<=0.15) &
                            values(NDVI>=0.1) & values(NDVI<=0.28) &
                            values(Z.om<=0.005) & values(Ts>(maxT-deltaTemp))),n)
      }
    }
    else{
      hot <- as.numeric(hc.data$pixel[hc.data$type =="hot"])
      cold <- as.numeric(hc.data$pixel[hc.data$type =="cold"])
    }
  print("Cold pixels")
  print(data.frame(cbind("LAI"=LAI[cold], "NDVI"=NDVI[cold], 
                         "albedo"=albedo[cold], "Z.om"=Z.om[cold])))
  print("Hot pixels")
  print(data.frame(cbind("LAI"=LAI[hot], "NDVI"=NDVI[hot], 
                         "albedo"=albedo[hot], "Z.om"=Z.om[hot])))


### End anchors selection ------------------------------------------------------
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
    LE.cold <- ETo.hourly * ETp.coef * (2.501 - 0.002361*(Ts[cold]-273.15))*(1e6)/3600 
    # here uses latent.heat.vapo
    H.cold <- Rn[cold] - G[cold] - LE.cold #ok
    result <- list()
    print("starting conditions")
    print("Cold")
    print(data.frame(cbind("Ts"=Ts[cold], "Ts_datum"=Ts.datum[cold], "Rn"=Rn[cold], 
                           "G"=G[cold], "Z.om"=Z.om[cold], "u200"=u200[cold], 
                           "u*"=friction.velocity[cold])))
    print("Hot")
    print(data.frame(cbind("Ts"=Ts[hot], "Ts_datum"=Ts.datum[hot], "Rn"=Rn[hot], 
                           "G"=G[hot], "Z.om"=Z.om[hot], "u200"=u200[hot], 
                           "u*"=friction.velocity[hot])))
    
    converge <- FALSE
    last.loop <- FALSE 
    i <- 1
### Start of iterative process -------------------------------------------------    
    while(!converge){
      i <-  i + 1 
      print(paste("iteraction #", i))
      ### We calculate dT and H 
      dT.cold <- H.cold * r.ah[cold] / (air.density[cold]*1004)
      dT.hot <- (Rn[hot] - G[hot]) * r.ah[hot] / (air.density[hot]*1004)
      a <- (dT.hot - dT.cold) / (Ts.datum[hot] - Ts.datum[cold])
      b <- -a * Ts.datum[cold] + dT.cold
      print(paste("a",a))
      print(paste("b",b))
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
      print(paste("r.ah cold", r.ah[cold]))
      print(paste("r.ah hot", r.ah[hot]))
      print(paste("dT cold", dT[cold]))
      print(paste("dT hot", dT[hot]))
      print("##############")
      ## And finally, r.ah and friction velocity
      friction.velocity <- 0.41 * u200 / (log(200/Z.om) - phi.200)
      # converge condition
      r.ah.hot.previous <- r.ah[hot]
      r.ah.cold.previous <- r.ah[cold]
      ### -----------
      r.ah <- (log(2/0.1) - phi.2 + phi.01) / (friction.velocity * 0.41) # ok ok
      # Check convergence
      if(last.loop == TRUE){
        converge <- TRUE
        print (paste0("convergence reached at iteration #", i))
        }
      delta.r.ah.hot <- (r.ah[hot] - r.ah.hot.previous) / r.ah[hot] * 100
      delta.r.ah.cold <- (r.ah[cold] - r.ah.cold.previous) / r.ah[cold] * 100
      print (paste("delta rah hot", delta.r.ah.hot))
      print (paste("delta rah cold", delta.r.ah.cold))
      print ("### -------")
      if(abs(delta.r.ah.hot) < 1 & abs(delta.r.ah.cold) < 1){last.loop <-  TRUE}
    } 
    
    
### End interactive process ----------------------------------------------------
    dT <- saveLoadClean(imagestack = dT, file = "dT", overwrite=TRUE)
    H <- saveLoadClean(imagestack = H, file = "H", overwrite=TRUE)
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