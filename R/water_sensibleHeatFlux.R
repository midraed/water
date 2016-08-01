#' Calculates Momentum Roughness Length
#' @description           this function estimates Momentum Roughness Length (Zom) from the average vegetation height around the weather station.
#' @param method          method selected to calculate momentum roughness length. Use 
#' "short.crops" for short crops methods from Allen et al (2007); "custom" for custom
#' method also in Allen et al (2007); Or "Perrier" to use Perrier equation as in 
#' Santos et al (2012) and Pocas et al (2014).
#' @param LAI             rasterLayer with Leaf Area Index. See LAI(). Only needed for method = "short.crops"
#' @param NDVI            rasterLayer with Normalized Difference Vegetation Index. Only needed for method = "custom"
#' @param albedo          broadband surface albedo. See albedo()
#' @param a               "a" coefficients for Allen (2007) custom function to estimate Momentum roughness length. Only needed for method = "custom"
#' @param b               "b" coefficients for Allen (2007) custom function to estimate Momentum roughness length. Only needed for method = "custom" 
#' @param fLAI.Perrier    proportion of LAI lying above h/2. Only needed for method = "Perrier"
#' @param h.Perrier       crop height in meters. Only needed for method = "Perrier"
#' @param mountainous      empirical adjustment for effects of general terrain roughness on momentum and heat transfer. See Allen (2007)
#' @param surface.model   surface model with a RasterLayer called "Slope" needed is mountainous = TRUE. See surface.model()
#' @details According Allen et al,. 2010 Zom is a measure of the form drag and skin friction for the layer of air that interacts with the surface.
#' @author Guillermo Federico Olmedo
#' @author de la Fuente-Saiz, Daniel
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#' 
#' Pocas, I., Paco, T.A., Cunha, M., Andrade, J.A., Silvestre, J., Sousa, A., Santos, F.L., Pereira, L.S., Allen, R.G., 2014. Satellite-based evapotranspiration of a super-intensive olive orchard: Application of METRIC algorithms. Biosystems Engineering 128, 69-81. doi:10.1016/j.biosystemseng.2014.06.019 \cr
#'
#' Santos, C., Lorite, I.J., Allen, R.G., Tasumi, M., 2012. Aerodynamic Parameterization of the Satellite-Based Energy Balance (METRIC) Model for ET Estimation in Rainfed Olive Orchards of Andalusia, Spain. Water Resour Manage 26, 3267-3283. doi:10.1007/s11269-012-0071-8 \cr
#' @export
## Create a function to estimate a and b coefficients or the function between Z.om and NDVI
## using some points and tabulated z.om for their covers.
## Perrier by Santos 2012 and Pocas 2014.
momentumRoughnessLength <- function(method="short.crops", LAI, NDVI, 
                                    albedo, a, b, fLAI.Perrier, h.Perrier, 
                                    mountainous=FALSE, surface.model){
  if(method=="short.crops"){
    Z.om <- (0.018*LAI)
  }
  if(method=="custom"){
    Z.om <- exp((a*NDVI/albedo)+b)
  }
  if(method=="Perrier"){
    if(fLAI.Perrier >=0.5){ a <- (2*(1-fLAI.Perrier))^-1 }
    if(fLAI.Perrier <0.5){ a <- 2*fLAI.Perrier }
    Z.om <- ((1-exp(-a*LAI/2))*exp(-a*LAI/2))^h.Perrier
  }
  if(mountainous==TRUE){
    Z.om <- Z.om * (1 + (180/pi*surface.model$Slope - 5)/20)
  }
  Z.om <- saveLoadClean(imagestack = Z.om, 
                        file = "Z.om", overwrite=TRUE)
  return(Z.om)
}

#' Select anchors pixels for H function 
#' @description            automatically search end members within the satellite scene (extreme wet and dry conditions).
#' @param image            top-of-atmosphere landsat reflectance image
#' @param Ts               land surface temperature in K. See surfaceTemperature()
#' @param LAI              rasterLayer with Leaf Area Index. See LAI()
#' @param albedo           broandband surface albedo. See albedo()
#' @param Z.om             momentum roughness lenght. See momentumRoughnessLength()
#' @param n                number of pair of anchors pixels to calculate
#' @param aoi              area of interest to limit the search. If 
#' waterOptions(autoAOI) == TRUE, It'll use aoi object from .GlobalEnv
#' @param anchors.method   method for the selection of anchor pixels. "CITRA-MCBr" for
#' random selection of hot and cold candidates according to CITRA-MCB method, or 
#' "CITRA-MCBbc" for selecting the best candidates
#' @param WeatherStation Optional. WeatherStation data at the satellite overpass. 
#' Should be a waterWeatherStation object calculated using read.WSdata and MTL file.
#' If you provide a WeatherStation object it will restrict the location of anchors
#' pixels to less than 30 km away from it.
#' @param plots            Logical. If TRUE will plot position of anchors points
#' selected. Points in red are selected hot pixels, blue are the cold ones and the 
#' black represents the position of the Weather Station
#' @param deltaTemp        deltaTemp for method "CITRA-MCBs" or "CITRA-MCBr"
#' @param buffer           minimun distance allowed for two anchor pixels of the same kind
#' @param verbose          Logical. If TRUE will print aditional data to console
#' @author Guillermo Federico Olmedo
#' @author de la Fuente-Saiz, Daniel
#' @references 
#' CITRA y MCB (com pers)
#' @export
calcAnchors  <- function(image, Ts, LAI, albedo, Z.om, n=1, aoi,
                         anchors.method= "CITRA-MCBbc", WeatherStation,
                         plots=TRUE, deltaTemp=5, buffer = 500, verbose=FALSE) {
  ### Some values used later
  NDVI <- (image$NIR - image$R) / (image$NIR + image$R)
  if(!missing(WeatherStation)){
    WSloc <- WeatherStation$location
    coordinates(WSloc) <- ~ long + lat
    WSloc@proj4string <- sp::CRS("+init=epsg:4326")
    WSloc <- sp::spTransform(WSloc, Ts@crs)
    ## Longer but avoids to use rgeos
    WScell <- extract(Ts, WSloc, cellnumbers=T)[1]
    WSbuffer <- raster(Ts)
    values(WSbuffer)[WScell] <- 1
    WSbuffer <- buffer(WSbuffer, width = 30000)
  } else {
    WSbuffer <- raster(Ts)
    values(WSbuffer)<- 1
  }
  if(anchors.method=="CITRA-MCBr"){
    minT <- quantile(Ts[LAI>=3&LAI<=6&albedo>=0.18&albedo<=0.25&Z.om>=0.03&
                          Z.om<=0.08], 0.05, na.rm=TRUE)
    if(minT+deltaTemp<288){minT = 288 + deltaTemp}
    ## NDVI used in cold isn't the same as CITRA!
    maxT <- max(Ts[albedo>=0.13&albedo<=0.15&NDVI>=0.1&NDVI<=0.28&
                     Z.om<=0.005], na.rm=TRUE)
    cold.candidates <- values(LAI>=3) & values(LAI<=6) &  
      values(albedo>=0.18) & values(albedo<=0.25) &
      values(NDVI>=max(values(NDVI), na.rm=T)-0.15) &
      values(Z.om>=0.03) & values(Z.om<=0.08) &
      values(Ts<(minT+deltaTemp)) & values(WSbuffer == 1)
    hot.candidates <- values(albedo>=0.13) & values(albedo<=0.15) &
      values(NDVI>=0.1) & values(NDVI<=0.28) &
      values(Z.om<=0.005) & values(Ts>(maxT-deltaTemp)) & values(WSbuffer == 1)
    # First cold sample
    try(cold <- sample(which(cold.candidates),1), silent=TRUE)
    if(!exists("cold")){
      stop("There are no pixels that meet the conditions for cold pixels")
    }
    if(n>1){  ## Next samples...
      for(nsample in 1:(n-1)){
        distbuffer <- raster(Ts)
        values(distbuffer)[cold] <- 1
        distbuffer <- buffer(distbuffer, width = buffer) ### 500m buffer
        distbuffer <- is.na(distbuffer)
        newAnchor <- NA
        cold.candidates <- values(LAI>=3) & values(LAI<=6) &  
          values(albedo>=0.18) & values(albedo<=0.25) &
          values(NDVI>=max(values(NDVI), na.rm=T)-0.15) &
          values(Z.om>=0.03) & values(Z.om<=0.08) &
          values(Ts<(minT+deltaTemp)) & values(distbuffer==1) & values(WSbuffer == 1)
        if(length(which(cold.candidates))<2){
          warning(paste("I can only find ", nsample, " anchors with cold pixel conditions"))
          break
        }
        try(newAnchor <- sample(which(cold.candidates),1), silent = FALSE)
        if(!is.na(newAnchor)){cold <- c(cold, newAnchor)} 
      }}
    
    # First hot sample
    try(hot <- sample(which(hot.candidates & values(Ts>quantile(Ts[hot.candidates], 0.75))),1), silent=TRUE)
    if(!exists("hot")){
      stop("There are no pixels that meet the conditions for hot pixels")
    }
    if(n>1){  ## Next samples...
      for(nsample in 1:(n-1)){
        distbuffer <- raster(Ts)
        values(distbuffer)[hot] <- 1
        distbuffer <- buffer(distbuffer, width = buffer) ### 500m buffer
        distbuffer <- is.na(distbuffer)
        newAnchor <- NA
        hot.candidates <- values(albedo>=0.13) & values(albedo<=0.15) &
          values(NDVI>=0.1) & values(NDVI<=0.28) & values(distbuffer==1) &
          values(Z.om<=0.005) & values(Ts>(maxT-deltaTemp)) & values(WSbuffer == 1)
        if(length(which(hot.candidates))<2){
          warning(paste("I can only find ", nsample, " anchors with hot pixel conditions"))
          break
        }
        try(newAnchor <- sample(which(hot.candidates),1), silent = FALSE)
        if(!is.na(newAnchor)){hot <- c(hot, newAnchor)} 
      }}
  }
  if(anchors.method=="CITRA-MCB"  | anchors.method=="CITRA-MCBbc"){
    minT <- quantile(Ts[LAI>=3&LAI<=6&albedo>=0.18&albedo<=0.25&Z.om>=0.03&
                          Z.om<=0.08], 0.05, na.rm=TRUE)
    if(minT+deltaTemp<288){minT = 288 + deltaTemp}
    ## NDVI used in cold isn't the same as CITRA!
    maxT <- max(Ts[albedo>=0.13&albedo<=0.15&NDVI>=0.1&NDVI<=0.28&
                     Z.om<=0.005], na.rm=TRUE)
    cold.candidates <- values(LAI>=3) & values(LAI<=6) &  
      values(albedo>=0.18) & values(albedo<=0.25) &
      values(NDVI>=max(values(NDVI), na.rm=T)-0.15) &
      values(Z.om>=0.03) & values(Z.om<=0.08) &
      values(Ts<(minT+deltaTemp)) & values(WSbuffer == 1)
    hot.candidates <- values(albedo>=0.13) & values(albedo<=0.15) &
      values(NDVI>=0.1) & values(NDVI<=0.28) &
      values(Z.om<=0.005) & values(Ts>(maxT-deltaTemp)) & values(WSbuffer == 1)
    # Cold samples
    Ts.cold <- Ts
    values(Ts.cold)[!cold.candidates] <- NA
    cold <- raster::which.min(Ts.cold)[1]
    if(n>1){  ## Next samples...
      for(nsample in 1:(n-1)){
        distbuffer <- raster(Ts)
        values(distbuffer)[cold] <- 1
        distbuffer <- buffer(distbuffer, width = buffer) ### 500m buffer
        distbuffer <- is.na(distbuffer)
        newAnchor <- NA
        cold.candidates <- values(LAI>=3) & values(LAI<=6) &  
          values(albedo>=0.18) & values(albedo<=0.25) &
          values(NDVI>=max(values(NDVI), na.rm=T)-0.15) &
          values(Z.om>=0.03) & values(Z.om<=0.08) &
          values(Ts<(minT+deltaTemp)) & values(distbuffer==1) & values(WSbuffer == 1)
        values(Ts.cold)[!cold.candidates] <- NA
        if(length(which(cold.candidates))<2){
          warning(paste("I can only find ", nsample, " anchors with cold pixel conditions"))
          break
        }
        try(newAnchor <- raster::which.min(Ts.cold)[1], silent = FALSE)
        if(!is.na(newAnchor)){cold <- c(cold, newAnchor)} 
      }}
    
    # hot samples
    Ts.hot <- Ts
    values(Ts.hot)[!hot.candidates] <- NA
    hot <- raster::which.max(Ts.hot)
    if(n>1){  ## Next samples...
      for(nsample in 1:(n-1)){
        distbuffer <- raster(Ts)
        values(distbuffer)[hot] <- 1
        distbuffer <- buffer(distbuffer, width = buffer) ### 500m buffer
        distbuffer <- is.na(distbuffer)
        newAnchor <- NA
        hot.candidates <- values(albedo>=0.13) & values(albedo<=0.15) &
          values(NDVI>=0.1) & values(NDVI<=0.28) & values(distbuffer==1) &
          values(Z.om<=0.005) & values(Ts>(maxT-deltaTemp)) & values(WSbuffer == 1)
        values(Ts.hot)[!hot.candidates] <- NA
        if(length(which(hot.candidates))<2){
          warning(paste("I can only find ", nsample, " anchors with hot pixel conditions"))
          break
        }
        try(newAnchor <- raster::which.max(Ts.hot)[1], silent = FALSE)
        if(!is.na(newAnchor)){hot <- c(hot, newAnchor)} 
      }}
    
  }
  if(verbose==TRUE){
    print("Cold pixels")
    print(data.frame(cbind(pixel=cold, "LAI"=LAI[cold], "NDVI"=NDVI[cold], 
                           "albedo"=albedo[cold], "Z.om"=Z.om[cold]), Ts=Ts[cold]))
    print("Hot pixels")
    print(data.frame(cbind(pixel=hot, "LAI"=LAI[hot], "NDVI"=NDVI[hot], 
                           "albedo"=albedo[hot], "Z.om"=Z.om[hot], Ts=Ts[hot])))    
  }
  ### End anchors selection ------------------------------------------------------
  ### We can plot anchor points
  if(plots==TRUE){
    plot(LAI, main="LAI and hot and cold pixels")
    graphics::points(xyFromCell(LAI, hot), col="red", pch=3)
    graphics::points(xyFromCell(LAI, cold), col="blue", pch=4)
    graphics::points(WSloc, pch=13)
  }
  hot.and.cold <- data.frame(pixel=integer(),  X=integer(), 
                             Y=integer(), Ts=double(), LAI=double(), 
                             type=factor(levels = c("hot", "cold")))
  for(i in 1:length(hot)){hot.and.cold[i, ] <- c(hot[i], xyFromCell(LAI, hot[i]),
                                                 Ts[hot][i], round(LAI[hot][i],2), "hot")}
  for(i in 1:length(cold)){hot.and.cold[i+length(hot), ] <- c(cold[i], xyFromCell(LAI, cold[i]), 
                                                              Ts[cold][i], round(LAI[cold][i],2), "cold")}
  for(i in 1:5){
    hot.and.cold[,i] <- as.numeric(hot.and.cold[,i])
  }
  return(hot.and.cold)
}



#' Iterative function to estimate H and R.ah
#' @description          generates an iterative solution to estimate r.ah and H because both are unknown at each pixel.
#' @param anchors        anchors points. Can be the result from calcAnchors() or
#' a spatialPointDataframe o Dataframe with X, Y, and type. type should be 
#' "cold" or "hot"
#' @param method         Method when using more than 1 pair of anchors pixels. 
#' method = "mean" will use the mean value for the cold pixels vs the mean value
#' for the hot pixels.
#' @param Ts             Land surface temperature in K. See surfaceTemperature()
#' @param Z.om           momentum roughness lenght. See momentumRoughnessLength()
#' @param WeatherStation WeatherStation data at the satellite overpass. 
#' Can be a waterWeatherStation object calculated using read.WSdata and MTL file
#' @param ETp.coef       ETp coefficient usually 1.05 or 1.2 for alfalfa
#' @param Z.om.ws        momentum roughness lenght for WeatherStation. Usually
#' a value of 0.03 might be reasonable for a typical agricultural weather station 
#' sited over vegetation that is about 0.3 m tall.  For clipped grass, use 0.015 m
#' @param mountainous    Logical. If TRUE heat transfer equation will be 
#' adjusted for mountainous terrain
#' @param DEM            Digital Elevation Model in meters.
#' @param Rn             Net radiation. See netRadiation()
#' @param G              Soil Heat Flux. See soilHeatFlux()
#' @param verbose        Logical. If TRUE will print information about every 
#' iteration to console
#' @param maxit          Maximun number of iteration. Default 20.
#' @details Sensible heat flux is the rate of heat loss to the air by convection 
#' and conduction, due to a temperature difference.This parameter is computed using the following one-dimensional,
#'aerodynamic,temperature gradient based equation for heat transport, this method is difficult to solve because 
#'there are two unknowns, rah and dT. To facilitate this computation, METRIC utilize the two "anchor" 
#'pixels  and solve for dT that satisfies eq. given the aerodynamic roughness and wind speed at a given height.
#'Aerodynamic resistance, and heat transfer is impacted by buoyancy of heated, light air at the surface, 
#'especially when H is large. Therefore, correction to rah is needed to account for buoyancy effects. However, 
#'H is needed to make this correction. An iterative solution for both H and rah is used.
#' @author Guillermo Federico Olmedo
#' @author de la Fuente-Saiz, Daniel
#' @author Fernando Fuentes Peñailillo
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#'
#' Allen, R., Irmak, A., Trezza, R., Hendrickx, J.M.H., Bastiaanssen, W., Kjaersgaard, J., 2011. Satellite-based ET estimation in agriculture using SEBAL and METRIC. Hydrol. Process. 25, 4011-4027. doi:10.1002/hyp.8408 \cr
#' @export
calcH  <- function(anchors, method = "mean", Ts, Z.om, WeatherStation, ETp.coef= 1.05, 
                   Z.om.ws=0.03, mountainous=FALSE, 
                   DEM, Rn, G, verbose=FALSE, maxit = 20) {
  if(class(WeatherStation)== "waterWeatherStation"){
    WeatherStation <- getDataWS(WeatherStation)
  }
  if(class(anchors) != "SpatialPointsDataFrame"){
    coordinates(anchors) <- ~ X + Y  
  }
  hot <- as.numeric(extract(Ts, anchors[anchors@data$type=="hot",], 
                            cellnumbers=T)[,1])
  cold <- as.numeric(extract(Ts, anchors[anchors@data$type=="cold",], 
                             cellnumbers=T)[,1])
  ###
  ETo.hourly <- hourlyET(WeatherStation, hours = WeatherStation$hours, 
                         DOY = WeatherStation$DOY, ET.instantaneous = TRUE,
                         ET= "ETor")
  Ts.datum <- Ts - (DEM - WeatherStation$elev) * 6.49 / 1000
  P <- 101.3*((293-0.0065 * DEM)/293)^5.26
  air.density <- 1000 * P / (1.01*(Ts)*287)
  latent.heat.vaporization <- (2.501-0.00236*(Ts-273.15))# En el paper dice por 1e6
  ### We calculate the initial conditions assuming neutral stability
  u.ws <- WeatherStation$wind * 0.41 / log(WeatherStation$height/Z.om.ws)
  u200.v <- u.ws / 0.41 * log(200/Z.om.ws)
  if(u200.v < 1){warning(paste0("u200 less than threshold value = ", 
                                round(u200.v,4), "m/s. using u200 = 4m/s"))
    u200.v <- 4
  }
  u200 <- raster(DEM)
  values(u200) <- u200.v
  if(mountainous==TRUE){
    u200 <- u200 * (1+0.1*((DEM-WeatherStation$elev)/1000))
  }
  friction.velocity <- 0.41 * u200 / log(200/Z.om) 
  r.ah <- log(2/0.1)/(friction.velocity*0.41) #ok
  
  LE.cold <- ETo.hourly * ETp.coef * (2.501 - 0.002361*(mean(Ts[cold])-273.15))*
    (1e6)/3600 
  # here uses latent.heat.vapo
  H.cold <- mean(Rn[cold]) - mean(G[cold]) - LE.cold #ok
  result <- list()
  if(verbose==TRUE){
    print("starting conditions")
    print("Cold")
    print(data.frame(cbind("Ts"=mean(Ts[cold]), "Ts_datum"=mean(Ts.datum[cold]), 
                           "Rn"=mean(Rn[cold]), "G"=mean(G[cold]), "Z.om"=mean(Z.om[cold]), 
                           "u200"=u200[cold], "u*"=mean(friction.velocity[cold]) )))
    print("Hot")
    print(data.frame(cbind("Ts"=mean(Ts[hot]), "Ts_datum"=mean(Ts.datum[hot]), "Rn"=mean(Rn[hot]), 
                           "G"=mean(G[hot]), "Z.om"=mean(Z.om[hot]), "u200"=u200[hot], 
                           "u*"=mean(friction.velocity[hot]))))
  }
  plot(1, mean(r.ah[hot]), xlim=c(0,15), ylim=c(0, mean(r.ah[hot])), 
       col="red", ylab="aerodynamic resistance s m-1", xlab="iteration", pch=20)
  graphics::points(1, mean(r.ah[cold]), col="blue", pch=20)
  converge <- FALSE
  last.loop <- FALSE 
  i <- 1
  if(method == "mean"){
    ### Start of iterative process -------------------------------------------------    
    while(!converge){
      #     ## For meta functions like METRIC.EB
      #     if(exists(x = "on.meta", envir=METRIC.EB)){
      #       if(i == 3){setTxtProgressBar(pb, 52)}
      #       if(i == 5){setTxtProgressBar(pb, 65)}
      #       if(i == 9){setTxtProgressBar(pb, 85)}
      #     }
      i <-  i + 1 
      if(verbose==TRUE){
        print(paste("iteraction #", i))
      }
      ### We calculate dT and H 
      dT.cold <- H.cold * mean(r.ah[cold]) / (mean(air.density[cold])*1004)
      dT.hot <- (mean(Rn[hot]) - mean(G[hot])) * mean(r.ah[hot]) / (mean(air.density[hot])*1004)
      a <- (dT.hot - dT.cold) / (mean(Ts.datum[hot]) - mean(Ts.datum[cold]))
      b <- -a * mean(Ts.datum[cold]) + dT.cold
      if(verbose==TRUE){
        print(paste("a",a))
        print(paste("b",b))
      }
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
      phi.200[Monin.Obukhov.L < 0] <- (2 * log((1+x.200)/2) + 
                                         log((1 + x.200^2) /2) - 
                                         2* atan(x.200) + 0.5 * pi)[Monin.Obukhov.L < 0] #ok
      phi.2[Monin.Obukhov.L < 0] <- (2 * log((1 + x.2^2) / 2))[Monin.Obukhov.L < 0]
      phi.01[Monin.Obukhov.L < 0] <- (2 * log((1 + x.01^2) / 2))[Monin.Obukhov.L < 0]
      if(verbose==TRUE){
        print(paste("r.ah cold", mean(r.ah[cold])))
        print(paste("r.ah hot", mean(r.ah[hot])))
        print(paste("dT cold", mean(dT[cold])))
        print(paste("dT hot", mean(dT[hot])))
        print("##############")
      }
      ## And finally, r.ah and friction velocity
      friction.velocity <- 0.41 * u200 / (log(200/Z.om) - phi.200)
      # converge condition
      r.ah.hot.previous <- mean(r.ah[hot])
      r.ah.cold.previous <- mean(r.ah[cold])
      ### -----------
      r.ah <- (log(2/0.1) - phi.2 + phi.01) / (friction.velocity * 0.41) # ok ok
      ## Update plot
      graphics::points(i, mean(r.ah[hot]), col="red", pch=20)
      graphics::points(i, mean(r.ah[cold]), col="blue", pch=20)
      lines(c(i, i-1), c(mean(r.ah[hot]), r.ah.hot.previous), col="red")
      lines(c(i, i-1), c(mean(r.ah[cold]), r.ah.cold.previous), col="blue")
      # Check convergence
      if(last.loop == TRUE){
        converge <- TRUE
        if(verbose==TRUE){print (paste0("convergence reached at iteration #", i))}
      }
      delta.r.ah.hot <- (mean(r.ah[hot]) - r.ah.hot.previous) / mean(r.ah[hot]) * 100
      delta.r.ah.cold <- (mean(r.ah[cold]) - r.ah.cold.previous) / mean(r.ah[cold]) * 100
      if(verbose==TRUE){
        print (paste("delta rah hot", delta.r.ah.hot))
        print (paste("delta rah cold", delta.r.ah.cold))
        print ("### -------")
      }
      if(abs(delta.r.ah.hot) < 1 & abs(delta.r.ah.cold) < 1){last.loop <-  TRUE}
      if(i == maxit){warning(paste0("maxit reached. Not solution found at iterarion #", i))
        break}
    } 
    ### End interactive process --------------------------------------------------
  }
  dT <- saveLoadClean(imagestack = dT, file = "dT", overwrite=TRUE)
  H <- saveLoadClean(imagestack = H, file = "H", overwrite=TRUE)
  result$a <- a
  result$b <- b
  result$dT <- dT
  result$H <- H
  return(result)
}
