aerodynamic.transport.0 <- function(Z.om, wind, height.ws=2, Z.omw = 0.0018, mountainous=FALSE, ws.elevation, surface.model){
  u200 <- wind * log(200/Z.omw)/log(height.ws/Z.omw)
  if(mountainous==TRUE){
    u200 <- u200 * (1+0.1*((surface.model$DEM-ws.elevation)/1000))
  }
  friction.velocity <- 0.41 * u200 / log(200/Z.om)
  r.ah.0 <- log(2/0.1)/friction.velocity*0.41
}

hot.and.cold <- function(method="random", n=1, Ts, ETr, ETp.coef= 1.05, Rn, G, r.ah = r.ah.0, DEM, LAI, aoi){
  Ts.datum <- Ts - (DEM - 702) * 6.49 / 1000
  P <- 101.3*((293-0.0065 * DEM)/293)^5.26
  air.density <- 1000 * P / (1.01*(Ts)*287)
  if(method=="random"){
    if(!max(values(Ts))>=310){
      warning(paste("Ts max value is", round(max(values(Ts)),2), "and the expected was Ts >= 310 for hot pixel"))
      hot <- sample(which(values(Ts)>quantile(Ts,0.9999)),n)  
    } else hot <- sample(which(values(Ts>310)),n) 
    if(!max(values(LAI))>=4){
      warning(paste("LAI max value is", round(max(values(LAI)),2), "and the expected was LAI >= 4 for cold pixel"))
      cold <- sample(which(values(LAI)>quantile(LAI,0.9999)),n)  
    } else cold <- sample(which(values(LAI>4)),n)  
  }
  dT.hot <- (Rn[hot] - G[hot])*r.ah[hot]/(air.density[hot]*1007)
  latent.heat.vaporization <- (2.501-0.00236*(Ts-273.15))  # En el paper dice por 1e6
  LE.cold <- ETp.coef * ETr * latent.heat.vaporization[cold] # instead of ETr better use ETr ~ LAI
  H.cold <- Rn[cold]-G[cold]-LE.cold
  dT.cold <- H.cold*r.ah[cold]/(air.density[cold]*1007)
  a <- mean((dT.hot-dT.cold)/(Ts.datum[hot]-Ts.datum[cold]), na.rm=T)
  b <- mean((dT.hot-a)/Ts.datum[hot], na.rm=T)
  dT <- a+b*Ts.datum
  # Prepare and populate result list
  result <- list()
  result$dT <- dT
  result$hot.and.cold <- data.frame(pixel=integer(),  X=integer(), Y=integer(), Ts=double(), 
                                    LAI=double(), type=factor(levels = c("hot", "cold")))
  for(i in 1:n){result$hot.and.cold[i, ] <- c(hot[i], xyFromCell(LAI, hot[i]), Ts[hot][i], round(LAI[hot][i],2), "hot")}
  for(i in 1:n){result$hot.and.cold[i+n, ] <- c(cold[i], xyFromCell(LAI, hot[i]), Ts[cold][i], round(LAI[cold][i],2), "cold")}
  removeTmpFiles(h=0)
  plot(LAI, main="LAI and hot and cold pixels")
  points(xyFromCell(LAI, hot), col="red", pch=3)
  points(xyFromCell(LAI, cold), col="blue", pch=4)
  return(result)
}

## Add iteractive solution for r.ah
## Add iteraction between r.ah.0 and r.ah
aerodynamic.transport <- function(Z.om, wind, height.ws=2, Z.omw = 0.0018, dT, mountainous=FALSE, ws.elevation, surface.model){
  u200 <- wind * log(200/Z.omw)/log(height.ws/Z.omw)
  if(mountainous==TRUE){
    u200 <- u200 * (1+0.1*((surface.model$DEM-ws.elevation)/1000))
  }
  friction.velocity.0 <- 0.41 * u200 / log(200/Z.om)
  r.ah.0 <- log(2/0.1)/friction.velocity.0*0.41
  P <- 101.3*((293-0.0065 * surface.model$DEM)/293)^5.26
  air.density <- 1000 * P / (1.01*(Ts)*287)
  H <- air.density*1007*(dT/r.ah.0)
  Monin.Obukhov.L <- (air.density * 1007 * friction.velocity.0^3 * Ts) / (0.41 * 9.807 * H)
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
  ##
  friction.velocity <- 0.41 * u200 / (log(200/Z.om) - phi.200)
  r.ah <- (log(2/0.1) - phi.2 + phi.01) / friction.velocity * 0.41
  r.ah[r.ah > quantile(r.ah, 0.999)] <-  NA # to eliminate very high values
  return(r.ah)
}

update.dT <- function(hc.data, Ts, DEM, Rn, G, r.ah, ETr, ETp.coef= 1.05){
  hot <- as.numeric(hc.data$pixel[hc.data$type =="hot"])
  cold <- as.numeric(hc.data$pixel[hc.data$type =="cold"])
  Ts.datum <- Ts - (DEM - 702) * 6.49 / 1000
  P <- 101.3*((293-0.0065 * DEM)/293)^5.26
  air.density <- 1000 * P / (1.01*(Ts)*287)
  dT.hot <- (Rn[hot] - G[hot])*r.ah[hot]/(air.density[hot]*1007)
  latent.heat.vaporization <- (2.501-0.00236*(Ts-273.15))  # En el paper dice por 1e6
  LE.cold <- ETp.coef * ETr * latent.heat.vaporization[cold] # instead of ETr better use ETr ~ LAI
  H.cold <- Rn[cold]-G[cold]-LE.cold
  dT.cold <- H.cold*r.ah[cold]/(air.density[cold]*1007)
  a <- mean((dT.hot-dT.cold)/(Ts.datum[hot]-Ts.datum[cold]), na.rm=T)
  b <- mean((dT.hot-a)/Ts.datum[hot], na.rm=T)
  dT <- a+b*Ts.datum
  return(dT)
}

## Calculates again the air.density... maybe I should make a insolate function for air.density
sensible.heat.flux <- function(Ts, dem, dT, r.ah){
  P <- 101.3*((293-0.0065 * dem)/293)^5.26
  air.density <- 1000 * P / (1.01*(Ts)*287)
  air.density*1007*dT/r.ah
}