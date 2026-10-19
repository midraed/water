#' Calculates Momentum Roughness Length
#' @description           Estimates the momentum roughness length (Zom, m) for
#' each pixel from LAI ("short.crops"), from NDVI and albedo ("custom") or from
#' LAI and crop height ("Perrier").
#' @param method          method selected to calculate momentum roughness length. Use
#' "short.crops" for short crops methods from Allen et al (2007); "custom" for custom
#' method also in Allen et al (2007); Or "Perrier" to use Perrier equation as in
#' Santos et al (2012) and Pocas et al (2014).
#' @param LAI             SpatRaster with leaf area index, from LAI(). Needed for
#' method = "short.crops" and "Perrier"
#' @param NDVI            SpatRaster with the Normalized Difference Vegetation
#' Index. Only needed for method = "custom"
#' @param albedo          broadband surface albedo, from albedo(). Only needed
#' for method = "custom"
#' @param a               "a" coefficient for Allen (2007) custom function to estimate Momentum roughness length. Only needed for method = "custom"
#' @param b               "b" coefficient for Allen (2007) custom function to estimate Momentum roughness length. Only needed for method = "custom"
#' @param fLAI            proportion of LAI lying above h/2. Only needed for method = "Perrier"
#' @param h               crop height in meters. Only needed for method = "Perrier"
#' @param mountainous     Logical. If TRUE applies the empirical adjustment for
#' the effects of general terrain roughness on momentum and heat transfer
#' (Allen et al., 2007); needs surface.model
#' @param surface.model   SpatRaster from METRICtopo() with a layer named Slope.
#' Only needed if mountainous = TRUE
#' @details According to Allen et al. (2010), Zom is a measure of the form drag
#' and skin friction for the layer of air that interacts with the surface.
#' @return SpatRaster with the momentum roughness length (m).
#' @author Guillermo Federico Olmedo
#' @author Daniel de la Fuente-Sáiz
#' @family sensible heat flux functions
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
                                    albedo, a, b, fLAI, h,
                                    mountainous=FALSE, surface.model){
  checkChoice(method, waterChoices$Zom.method, "method")
  number <- function(x, name){
    if(!is.numeric(x) || length(x) != 1 || !is.finite(x)){
      stop(name, " must be a finite number", call. = FALSE)
    }
  }
  if(method=="short.crops"){
    Z.om <- (0.018*LAI)
  }
  if(method=="custom"){
    if(missing(a) || missing(b)){stop("method = \"custom\" needs a and b", call. = FALSE)}
    number(a, "a"); number(b, "b")
    Z.om <- exp((a*NDVI/albedo)+b)
  }
  if(method=="Perrier"){
    if(missing(fLAI) || missing(h)){stop("method = \"Perrier\" needs fLAI and h", call. = FALSE)}
    if(!is.numeric(fLAI) || length(fLAI) != 1 || !isTRUE(fLAI >= 0 && fLAI <= 1)){
      stop("fLAI must be a number between 0 and 1", call. = FALSE)
    }
    if(!is.numeric(h) || length(h) != 1 || !isTRUE(h > 0)){
      stop("h (crop height, m) must be a number above 0", call. = FALSE)
    }
    if(fLAI <0.5){ a <- (2*(1-fLAI))^-1 }
    if(fLAI >=0.5){ a <- 2*fLAI }
    Z.om <- ((1-exp(-a*LAI/2))*exp(-a*LAI/2))*h
  }
  if(mountainous==TRUE){
    Z.om <- Z.om * (1 + (180/pi*surface.model$Slope - 5)/20)
  }
  Z.om <- saveLoadClean(imagestack = Z.om,
                        file = "Z.om", overwrite=TRUE)
  return(Z.om)
}

#' Select anchor pixels for the H function
#' @description            Automatically searches for end members (anchor
#' pixels) within the satellite scene (extreme wet and dry conditions): cold
#' pixels with well-watered full-cover vegetation and hot pixels with dry bare
#' soil.
#' @param image            reflectance SpatRaster (top-of-atmosphere
#' reflectance, or surface reflectance when no DN image is available) with
#' layers named R and NIR, used to compute NDVI
#' @param Ts               land surface temperature in K. See \code{\link{surfaceTemperature}}
#' @param LAI              SpatRaster with leaf area index, from LAI()
#' @param albedo           broadband surface albedo, from albedo()
#' @param Z.om             momentum roughness length, from momentumRoughnessLength()
#' @param n                number of pairs of anchor pixels to select
#' @param aoi              Not used; kept for backward compatibility
#' @param anchors.method   method for the selection of anchor pixels: "flexible"
#' (default), "best" or "random". See Details.
#' @param WeatherStation   Optional. waterWeatherStation object (see
#' read.WSdata()). If supplied, the anchors are searched within WSbuffer of the
#' station
#' @param plots            Logical. If TRUE will plot position of anchors points
#' selected. Points in red are selected hot pixels, blue are the cold ones and the
#' black represents the position of the Weather Station (if supplied)
#' @param deltaTemp        temperature margin (K, default 5): cold candidates must
#' be colder than the 5\% quantile of Ts over well-vegetated pixels plus
#' deltaTemp, and hot candidates warmer than the maximum Ts over bare-soil
#' pixels minus deltaTemp
#' @param minDist          minimum distance (m) between two anchor pixels of the
#' same type
#' @param WSbuffer         maximum distance (m) to the weather station. If no
#' pixel is within WSbuffer, the whole image is used, with a warning
#' @param verbose          Logical. If TRUE will print additional data to console
#' @details
#' Cold candidates have LAI 3 to 6, albedo 0.18 to 0.25, Z.om 0.03 to 0.08, NDVI
#' within 0.15 of the maximum NDVI of the image and Ts below the threshold given
#' by deltaTemp. Hot candidates have albedo 0.13 to 0.15, NDVI 0.1 to 0.28,
#' Z.om up to 0.005 and Ts above the threshold given by deltaTemp.
#'
#' anchors.method = "best" selects the coldest and the hottest candidates;
#' "random" selects them at random among the candidates (CITRA-MCB method);
#' "flexible" (default) works like "best" but, when there are not enough
#' candidates, relaxes the limits of the criteria (LAI, albedo and Z.om for cold
#' pixels; albedo, NDVI, Z.om and Ts for hot pixels) in 1\% steps up to 100\%,
#' first near the weather station (within WSbuffer) and then in the whole
#' image, with a warning. The old names "CITRA-MCB" and "CITRA-MCBbc" (now
#' "best") and "CITRA-MCBr" (now "random") are deprecated.
#' @return A data.frame with one row per anchor pixel and columns pixel (cell
#' number), X and Y (cell coordinates in the CRS of Ts), Ts (K), LAI, NDVI and
#' type ("hot" or "cold"). For anchors.method = "flexible" a column flex is
#' added: the relaxation applied to the cold-pixel criteria (0 = none, 1 =
#' 100\%). If plots = TRUE the anchors are also plotted over LAI.
#' @author Guillermo Federico Olmedo
#' @author Daniel de la Fuente-Sáiz
#' @family sensible heat flux functions
#' @references
#' CITRA and MCB (personal communication)
#' @export
calcAnchors  <- function(image, Ts, LAI, albedo, Z.om, n=1, aoi,
                         anchors.method= "flexible", WeatherStation,
                         plots=TRUE, deltaTemp=5, minDist = 500, WSbuffer = 30000,
                         verbose=FALSE) {
  sanitizeCellIndex <- function(index) {
    as.integer(stats::na.omit(unlist(index, use.names = FALSE)))
  }
  # old method names (CITRA-MCBx) are mapped, with a warning
  anchors.method <- checkAnchorsMethod(anchors.method)
  checkKelvin(Ts)
  ### Some values used later
  rn <- landsatBands(image, c("R", "NIR"), pos = 3:4)
  NDVI <- (rn[[2]] - rn[[1]]) / (rn[[2]] + rn[[1]])
  NDVI[NDVI < -1]  <-  NA
  NDVI[NDVI > 1]  <-  NA
  if(!missing(WeatherStation)){
    WSloc <- WeatherStation$location
    WSloc <- sf::st_as_sf(WSloc, coords = c("long", "lat"), crs = 4326)
    WSloc <- sf::st_transform(WSloc, crs = terra::crs(Ts))
    # Cells within WSbuffer of the centre of the station cell or, when the
    # station lies outside the image, of the station itself
    WScell <- terra::cellFromXY(Ts, sf::st_coordinates(WSloc))
    WSpoint <- if(is.na(WScell)){terra::vect(WSloc)} else {
      terra::vect(terra::xyFromCell(Ts, WScell), crs = terra::crs(Ts))
    }
    WS.buffer <- terra::distance(Ts, WSpoint) <= WSbuffer
    if(!isTRUE(global(WS.buffer, "max", na.rm = TRUE)[1, 1] == 1)){
      warning("No pixel is within WSbuffer (", WSbuffer, " m) of the weather ",
              "station; searching anchors in the whole image")
      WS.buffer <- rast(Ts, vals = 1)
    }
  } else {
    WS.buffer <- rast(Ts, vals = 1)
  }
  if(anchors.method=="random"){
    minT <- quantile(Ts[LAI>=3&LAI<=6&albedo>=0.18&albedo<=0.25&Z.om>=0.03&
                          Z.om<=0.08], 0.05, na.rm=TRUE)
    if(minT+deltaTemp<288 | is.na(minT)){minT = 288 + deltaTemp}
    ## NDVI used in cold isn't the same as CITRA!
    maxT <- max(Ts[albedo>=0.13&albedo<=0.15&NDVI>=0.1&NDVI<=0.28&
                     Z.om<=0.005], na.rm=TRUE)
    cold.candidates <- values(LAI>=3) & values(LAI<=6) &
      values(albedo>=0.18) & values(albedo<=0.25) &
      values(NDVI>=max(values(NDVI), na.rm=T)-0.15) &
      values(Z.om>=0.03) & values(Z.om<=0.08) &
      values(Ts<(minT+deltaTemp)) & values(WS.buffer == 1)
    hot.candidates <- values(albedo>=0.13) & values(albedo<=0.15) &
      values(NDVI>=0.1) & values(NDVI<=0.28) &
      values(Z.om<=0.005) & values(Ts>(maxT-deltaTemp)) & values(WS.buffer == 1)
    ### Test # anchors
    cold.n <- sum(as.numeric(cold.candidates), na.rm = T)
    hot.n <- sum(as.numeric(hot.candidates), na.rm = T)
    if(cold.n < 1 | hot.n < 1){
      stop("Not enough pixels with the conditions for anchor pixels: found ",
           cold.n, " cold and ", hot.n, " hot candidates. ", anchorHint(albedo))
    }
    # First cold sample
    try(cold <- sampleOne(which(cold.candidates)), silent=TRUE)
    if(n>1){  ## Next samples...
      for(nsample in 1:(n-1)){
        distbuffer <- rast(Ts, vals = NA)
        values(distbuffer)[sanitizeCellIndex(cold)] <- 1

        # BUG in raster package
        # workaround https://gis.stackexchange.com/questions/264133/raster-buffer-error-with-package-updates/264154
        distbuffer <- terra::buffer(distbuffer, width = minDist) ### 500m buffer
        
        
        distbuffer <- !distbuffer  # terra::buffer() is TRUE inside, FALSE outside
        newAnchor <- NA
        cold.candidates <- values(LAI>=3) & values(LAI<=6) &
          values(albedo>=0.18) & values(albedo<=0.25) &
          values(NDVI>=max(values(NDVI), na.rm=T)-0.15) &
          values(Z.om>=0.03) & values(Z.om<=0.08) &
          values(Ts<(minT+deltaTemp)) & values(distbuffer==1) & values(WS.buffer == 1)
        if(length(which(cold.candidates))<2){
          warning(paste("I can only find ", nsample, " anchors with cold pixel conditions"))
          break
        }
        try(newAnchor <- sampleOne(which(cold.candidates)), silent = FALSE)
        if(!is.na(newAnchor)){cold <- c(cold, newAnchor)}
      }}

    # First hot sample
    hot.pool <- which(hot.candidates & values(Ts > quantile(Ts[hot.candidates], 0.75, na.rm = TRUE)))
    if(length(hot.pool) == 0){
      hot.pool <- which(hot.candidates)
    }
    try(hot <- sampleOne(hot.pool), silent=TRUE)
    if(n>1){  ## Next samples...
      for(nsample in 1:(n-1)){
        distbuffer <- rast(Ts, vals = NA)
        values(distbuffer)[sanitizeCellIndex(hot)] <- 1

        # BUG in raster package
        # workaround https://gis.stackexchange.com/questions/264133/raster-buffer-error-with-package-updates/264154
        distbuffer <- terra::buffer(distbuffer, width = minDist) ### 500m buffer
        
        distbuffer <- !distbuffer  # terra::buffer() is TRUE inside, FALSE outside
        newAnchor <- NA
        hot.candidates <- values(albedo>=0.13) & values(albedo<=0.15) &
          values(NDVI>=0.1) & values(NDVI<=0.28) & values(distbuffer==1) &
          values(Z.om<=0.005) & values(Ts>(maxT-deltaTemp)) & values(WS.buffer == 1)
        if(length(which(hot.candidates))<2){
          warning(paste("I can only find ", nsample, " anchors with hot pixel conditions"))
          break
        }
        try(newAnchor <- sampleOne(which(hot.candidates)), silent = FALSE)
        if(!is.na(newAnchor)){hot <- c(hot, newAnchor)}
      }}
  }
  if(anchors.method=="best"){
    minT <- quantile(Ts[LAI>=2.8&LAI<=6&albedo>=0.15&albedo<=0.25&Z.om>=0.03&
                          Z.om<=0.08], 0.05, na.rm=TRUE)
    if(minT+deltaTemp<288 | is.na(minT)){minT = 288 + deltaTemp}
    ## NDVI used in cold isn't the same as CITRA!
    maxT <- max(Ts[albedo>=0.13&albedo<=0.15&NDVI>=0.1&NDVI<=0.28&
                     Z.om<=0.005], na.rm=TRUE)
    cold.candidates <- values(LAI>=3) & values(LAI<=6) &
      values(albedo>=0.18) & values(albedo<=0.25) &
      values(NDVI>=max(values(NDVI), na.rm=T)-0.15) &
      values(Z.om>=0.03) & values(Z.om<=0.08) &
      values(Ts<(minT+deltaTemp)) & values(WS.buffer == 1)
    hot.candidates <- values(albedo>=0.13) & values(albedo<=0.15) &
      values(NDVI>=0.1) & values(NDVI<=0.28) &
      values(Z.om<=0.005) & values(Ts>(maxT-deltaTemp)) & values(WS.buffer == 1)
    ### Test # anchors
    cold.n <- sum(as.numeric(cold.candidates), na.rm = T)
    hot.n <- sum(as.numeric(hot.candidates), na.rm = T)
    if(cold.n < 1 | hot.n < 1){
      stop("Not enough pixels with the conditions for anchor pixels: found ",
           cold.n, " cold and ", hot.n, " hot candidates. ", anchorHint(albedo))
    }
    # Cold samples
    Ts.cold <- Ts
    values(Ts.cold)[!(cold.candidates %in% TRUE)] <- NA
    cold <- which.min(values(Ts.cold))[1]
    if(n>1){  ## Next samples...
      for(nsample in 1:(n-1)){
        distbuffer <- rast(Ts, vals = NA)
        values(distbuffer)[sanitizeCellIndex(cold)] <- 1

        # BUG in raster package
        # workaround https://gis.stackexchange.com/questions/264133/raster-buffer-error-with-package-updates/264154
        distbuffer <- terra::buffer(distbuffer, width = minDist) ### 500m buffer
        
        distbuffer <- !distbuffer  # terra::buffer() is TRUE inside, FALSE outside
        newAnchor <- NA
        cold.candidates <- values(LAI>=3) & values(LAI<=6) &
          values(albedo>=0.18) & values(albedo<=0.25) &
          values(NDVI>=max(values(NDVI), na.rm=T)-0.15) &
          values(Z.om>=0.03) & values(Z.om<=0.08) &
          values(Ts<(minT+deltaTemp)) & values(distbuffer==1) & values(WS.buffer == 1)
        values(Ts.cold)[!(cold.candidates %in% TRUE)] <- NA
        if(length(which(cold.candidates))<2){
          warning(paste("I can only find ", nsample, " anchors with cold pixel conditions"))
          break
        }
        try(newAnchor <- which.min(values(Ts.cold))[1], silent = FALSE)
        if(!is.na(newAnchor)){cold <- c(cold, newAnchor)} 
      }}

    # hot samples
    Ts.hot <- Ts
    values(Ts.hot)[!(hot.candidates %in% TRUE)] <- NA
    hot <- terra::which.max(values(Ts.hot))
    if(n>1){  ## Next samples...
      for(nsample in 1:(n-1)){
        distbuffer <- rast(Ts, vals = NA)
        values(distbuffer)[sanitizeCellIndex(hot)] <- 1

        # BUG in raster package
        # workaround https://gis.stackexchange.com/questions/264133/raster-buffer-error-with-package-updates/264154
        distbuffer <- terra::buffer(distbuffer, width = minDist) ### 500m buffer
        
        distbuffer <- !distbuffer  # terra::buffer() is TRUE inside, FALSE outside
        newAnchor <- NA
        hot.candidates <- values(albedo>=0.13) & values(albedo<=0.15) &
          values(NDVI>=0.1) & values(NDVI<=0.28) & values(distbuffer==1) &
          values(Z.om<=0.005) & values(Ts>(maxT-deltaTemp)) & values(WS.buffer == 1)
        values(Ts.hot)[!(hot.candidates %in% TRUE)] <- NA
        if(length(which(hot.candidates))<2){
          warning(paste("I can only find ", nsample, " anchors with hot pixel conditions"))
          break
        }
        try(newAnchor <- which.max(values(Ts.hot))[1], silent = FALSE)
        if(!is.na(newAnchor)){hot <- c(hot, newAnchor)} 
      }}

  }

  if(anchors.method=="flexible"){    ### method = "flexible" ####
    ### Find minT and maxT or fall back to default values
    minT <- quantile(Ts[LAI>=2.8&LAI<=6&albedo>=0.15&albedo<=0.25&Z.om>=0.03&
                          Z.om<=0.08], 0.05, na.rm=TRUE)
    if(minT+deltaTemp<288 | is.na(minT)){minT = 288 + deltaTemp}
    maxT <- suppressWarnings(max(Ts[albedo>=0.13&albedo<=0.15&NDVI>=0.1&NDVI<=0.28&
                     Z.om<=0.005], na.rm=TRUE))
    if(!is.finite(maxT)){maxT <- stats::quantile(values(Ts), 0.95, na.rm = TRUE)}
    ### create data.frames with optimal values for anchors
    optValCold <- data.frame(LAI = c(3,6), albedo = c(0.18, 0.25),
                                Z.om = c(0.03, 0.08))
    optValHot <- data.frame(albedo = c(0.13, 0.15), NDVI = c(0.1, 0.28),
                               Z.om = c(NA, 0.005), Ts = c(maxT-deltaTemp, NA))
    NDVI.cold <- max(values(NDVI), na.rm=T)-0.15
    coldCandidates <- function(opt, inBuffer){
      values(LAI>=opt$LAI[1]) & values(LAI<=opt$LAI[2]) &
        values(albedo>=opt$albedo[1]) & values(albedo<=opt$albedo[2]) &
        values(NDVI>=NDVI.cold) &
        values(Z.om>=opt$Z.om[1]) & values(Z.om<=opt$Z.om[2]) &
        values(Ts<(minT+deltaTemp)) & inBuffer
    }
    hotCandidates <- function(opt, inBuffer){
      values(albedo>=opt$albedo[1]) & values(albedo<=opt$albedo[2]) &
        values(NDVI>=opt$NDVI[1]) & values(NDVI<=opt$NDVI[2]) &
        values(Z.om<=opt$Z.om[2]) & values(Ts>(opt$Ts[1])) & inBuffer
    }
    ### Relax the criteria in 1% steps (up to 100%) until more than n candidates
    ### are found: first near the weather station, then in the whole image
    relaxCriteria <- function(candidates, opt, type){
      for(useBuffer in c(TRUE, FALSE)){
        inBuffer <- if(useBuffer) values(WS.buffer == 1) else TRUE
        optRelaxed <- opt
        flex <- 0
        found <- candidates(optRelaxed, inBuffer)
        while(sum(found, na.rm = TRUE) <= n && flex < 1){
          flex <- flex + 0.01
          print(paste0("relaxing criteria for ", type, " pixels: ", round(flex*100), "%"))
          optRelaxed[1,] <- opt[1,] * (1-flex)
          optRelaxed[2,] <- opt[2,] * (1+flex)
          found <- candidates(optRelaxed, inBuffer)
        }
        if(sum(found, na.rm = TRUE) > n){
          if(!useBuffer){
            warning("Not enough ", type, " anchor candidates within WSbuffer (",
                    WSbuffer, " m) of the weather station; searched the whole image")
          }
          if(flex != 0 || !useBuffer){
            warning(paste0("Criteria used for ", type, " pixels: ",
                           paste(names(optRelaxed), sapply(optRelaxed, function(x)
                             paste(round(x, 4), collapse = " to ")), collapse = "; "),
                           if(!useBuffer) " (whole image)"))
          }
          return(list(candidates = found, flex = flex, opt = optRelaxed,
                      inBuffer = inBuffer))
        }
      }
      stop("Automatic selection of ", type, " anchors FAILED even after relaxing ",
           "the criteria by 100% in the whole image. ", anchorHint(albedo))
    }
    coldSearch <- relaxCriteria(coldCandidates, optValCold, "cold")
    cold.candidates <- coldSearch$candidates
    flex.cold <- coldSearch$flex
    optValCold <- coldSearch$opt  # criteria used, also for the next anchors
    hotSearch <- relaxCriteria(hotCandidates, optValHot, "hot")
    hot.candidates <- hotSearch$candidates
    flex.hot <- hotSearch$flex
    optValHot <- hotSearch$opt
    # Cold samples
    Ts.cold <- Ts
    values(Ts.cold)[!(cold.candidates %in% TRUE)] <- NA
    cold <- terra::which.min(values(Ts.cold))[1]
    if(n>1){  ## Next samples...
      for(nsample in 1:(n-1)){
        distbuffer <- rast(Ts, vals = NA)
        values(distbuffer)[sanitizeCellIndex(cold)] <- 1

        # BUG in raster package
        # workaround https://gis.stackexchange.com/questions/264133/raster-buffer-error-with-package-updates/264154
        distbuffer <- terra::buffer(distbuffer, width = minDist) ### 500m buffer
        
        distbuffer <- !distbuffer  # terra::buffer() is TRUE inside, FALSE outside
        newAnchor <- NA
        cold.candidates <- values(LAI>=optValCold$LAI[1]) & values(LAI<=optValCold$LAI[2]) &
          values(albedo>=optValCold$albedo[1]) & values(albedo<=optValCold$albedo[2]) &
          values(NDVI>=max(values(NDVI), na.rm=T)-0.15) &
          values(Z.om>=optValCold$Z.om[1]) & values(Z.om<=optValCold$Z.om[2]) &
          values(Ts<(minT+deltaTemp)) & coldSearch$inBuffer &
          values(distbuffer==1)

        values(Ts.cold)[!(cold.candidates %in% TRUE)] <- NA
        if(length(which(cold.candidates))<=2){
          warning(paste("I can only find ", nsample, " anchors with cold pixel conditions"))
          break
        }
        try(newAnchor <- terra::which.min(values(Ts.cold))[1], silent = FALSE)
        if(!is.na(newAnchor)){cold <- c(cold, newAnchor)} 
      }}

    # hot samples
    Ts.hot <- Ts
    values(Ts.hot)[!(hot.candidates %in% TRUE)] <- NA
    hot <- which.max(values(Ts.hot))
    if(n>1){  ## Next samples...
      for(nsample in 1:(n-1)){
        distbuffer <- rast(Ts, vals = NA)
        values(distbuffer)[sanitizeCellIndex(hot)] <- 1

        # BUG in raster package
        # workaround https://gis.stackexchange.com/questions/264133/raster-buffer-error-with-package-updates/264154
        distbuffer <- terra::buffer(distbuffer, width = minDist) ### 500m buffer
        
        distbuffer <- !distbuffer  # terra::buffer() is TRUE inside, FALSE outside
        newAnchor <- NA
        hot.candidates <- values(albedo>=optValHot$albedo[1]) & values(albedo<=optValHot$albedo[2]) &
          values(NDVI>=optValHot$NDVI[1]) & values(NDVI<=optValHot$NDVI[2]) &
          values(Z.om<=optValHot$Z.om[2]) & values(Ts>(optValHot$Ts[1])) & hotSearch$inBuffer &
          values(distbuffer==1)
        values(Ts.hot)[!(hot.candidates %in% TRUE)] <- NA
        if(length(which(hot.candidates))<2){
          warning(paste("I can only find ", nsample, " anchors with hot pixel conditions"))
          break
        }
        try(newAnchor <- which.max(values(Ts.hot))[1], silent = FALSE)
        if(!is.na(newAnchor)){hot <- c(hot, newAnchor)} 
      }}

    }


  hot <- sanitizeCellIndex(hot)
  cold <- sanitizeCellIndex(cold)
  if(length(hot) == 0 || length(cold) == 0){
    stop("Automatic anchor selection failed to return valid hot/cold pixels. Try reducing n, increasing deltaTemp, or providing anchors manually.")
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
    graphics::points(terra::xyFromCell(LAI, sanitizeCellIndex(hot)), col="red", pch=3)
    graphics::points(terra::xyFromCell(LAI, sanitizeCellIndex(cold)), col="blue", pch=4)
    if(!missing(WeatherStation)){graphics::points(WSloc, pch=13)}
  }
  hot.and.cold <- data.frame(pixel=integer(),  X=integer(),
                             Y=integer(), Ts=double(), LAI=double(), NDVI=double(),
                             type=factor(levels = c("hot", "cold")))
  for(i in 1:length(hot)){hot.and.cold[i, ] <- c(hot[i], terra::xyFromCell(LAI, sanitizeCellIndex(hot[i])),
                                                 Ts[hot[i]][1,1], round(LAI[hot[i]][1,1],2), round(NDVI[hot[i]][1,1],2), "hot")}
  for(i in 1:length(cold)){hot.and.cold[i+length(hot), ] <- c(cold[i], terra::xyFromCell(LAI, sanitizeCellIndex(cold[i])),
                                                              Ts[cold[i]][1,1], round(LAI[cold[i]][1,1],2), round(NDVI[cold[i]][1,1],2), "cold")}
  for(i in 1:5){
    hot.and.cold[,i] <- as.numeric(hot.and.cold[,i])
  }
  if(anchors.method=="flexible"){
    hot.and.cold$flex <- flex.cold
    #hot.and.cold$flex[hot.and.cold$type == "hot"] <- flex.hot
  }
  return(hot.and.cold)
}



#' Iterative function to estimate H and R.ah
#' @description          Generates an iterative solution to estimate r.ah and H because both are unknown at each pixel.
#' @param anchors        anchor pixels: the data.frame returned by calcAnchors(),
#' another data.frame with X, Y (coordinates in the CRS of Ts) and type
#' columns, or an sf POINT object with a type column (transformed to the CRS
#' of Ts). type must be "hot" or "cold" (case is ignored). The anchors are
#' checked before iterating: they must be inside the image, Ts, Rn, G, Z.om
#' and DEM must not be NA there, and the hot anchors must be at least 5 K
#' warmer than the cold ones (Ts corrected to the datum elevation; for
#' method = "lm", in every pair). A smaller difference means that the area
#' lacks the wet and dry extremes that METRIC needs
#' @param method         Method when using more than 1 pair of anchor pixels:
#' "mean" (default) uses the mean values of the cold pixels and of the hot
#' pixels; "lm" iterates each hot/cold pair in turn (maxit applies to each
#' pair) and returns a and b of the last pair
#' @param Ts             Land surface temperature in K. See \code{\link{surfaceTemperature}}
#' @param Z.om           momentum roughness length, from momentumRoughnessLength()
#' @param WeatherStation WeatherStation data at the satellite overpass.
#' Can be a waterWeatherStation object calculated using read.WSdata and MTL file
#' @param ETp.coef       ratio of the ET of the cold pixels to the hourly alfalfa
#' reference ET (ETr) at the satellite overpass (see Details). Allen et al.
#' (2007) use 1.05, the default
#' @param Z.om.ws        momentum roughness length (m) of the weather station
#' surface, used to extrapolate the wind speed to 200 m (Allen et al., 2007,
#' Eq. 32). It is about 0.123 times the height of the vegetation around the
#' station (Allen et al., 1998): 0.015 m (the default) for a standard station
#' over clipped grass 0.12 m tall, 0.03 m for vegetation about 0.25 m tall
#' @param mountainous    Logical. If TRUE heat transfer equation will be
#' adjusted for mountainous terrain
#' @param DEM            Digital Elevation Model in meters.
#' @param Rn             Net radiation. See \code{\link{netRadiation}}
#' @param G              Soil Heat Flux. See \code{\link{soilHeatFlux}}
#' @param verbose        Logical. If TRUE will print information about every
#' iteration to console
#' @param maxit          Maximum number of iterations. Default 20.
#' @details Sensible heat flux is the rate of heat loss to the air by convection
#' and conduction, due to a temperature difference. This parameter is computed
#' using a one-dimensional, aerodynamic, temperature-gradient-based equation for
#' heat transport. This method is difficult to solve because there are two
#' unknowns, rah and dT. To facilitate this computation, METRIC utilizes the two
#' "anchor" pixels and solves for the dT that satisfies the equation given the
#' aerodynamic roughness and wind speed at a given height. Aerodynamic
#' resistance and heat transfer are impacted by buoyancy of heated, light air at
#' the surface, especially when H is large. Therefore, correction to rah is
#' needed to account for buoyancy effects. However, H is needed to make this
#' correction. An iterative solution for both H and rah is used.
#'
#' The ET of the cold pixels is computed as ETp.coef times the hourly alfalfa
#' reference ET (ETr, ASCE standardized equation; see \code{\link{hourlyET}}) at the
#' satellite overpass, as in Allen et al. (2007). ET24h() and dailyET() use the
#' same reference.
#' @return A list with a and b (coefficients of dT = a * Ts_datum + b), dT
#' (SpatRaster, near-surface temperature difference, K) and H (SpatRaster,
#' sensible heat flux, W/m2). The aerodynamic resistance of the anchor pixels at
#' each iteration is plotted.
#' @author Guillermo Federico Olmedo
#' @author Daniel de la Fuente-Sáiz
#' @author Fernando Fuentes-Peñailillo
#' @family sensible heat flux functions
#' @references
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007 \cr
#'
#' Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration: guidelines for computing crop water requirements. FAO Irrigation and Drainage Paper 56. FAO, Rome.
#'
#' Allen, R., Irmak, A., Trezza, R., Hendrickx, J.M.H., Bastiaanssen, W., Kjaersgaard, J., 2011. Satellite-based ET estimation in agriculture using SEBAL and METRIC. Hydrol. Process. 25, 4011-4027. doi:10.1002/hyp.8408 \cr
#' @export
calcH  <- function(anchors, method = "mean", Ts, Z.om, WeatherStation, ETp.coef= 1.05,
                   Z.om.ws=0.015, mountainous=FALSE,
                   DEM, Rn, G, verbose=FALSE, maxit = 20) {
  checkChoice(method, c("mean", "lm"), "method")

  # ---- helper: safe numeric values by cell (Option A) ----
  cell_vals <- function(r, cells) {
    x <- r[cells]
    
    # Case 1: terra returns a data.frame (like your r.ah[hot] output)
    if (is.data.frame(x)) {
      # common: one column with layer name (e.g., "B")
      return(as.numeric(x[[1]]))
    }
    
    # Case 2: terra returns a SpatRaster
    if (inherits(x, "SpatRaster")) {
      return(as.numeric(terra::values(x, mat = FALSE)))
    }
    
    stop("Unexpected return type from r[cells]: ", paste(class(x), collapse = ", "))
  }
  
  cell_mean <- function(r, cells, na.rm = TRUE) {
    mean(cell_vals(r, cells), na.rm = na.rm)
  }
  
  if(inherits(WeatherStation, "waterWeatherStation")){
    WeatherStation <- getDataWS(WeatherStation)
  }
  checkKelvin(Ts)

  # alfalfa reference ET, as in Allen et al. (2007) and in ET24h()
  ETr.hourly <- hourlyET(WeatherStation, hours = WeatherStation$hours,
                         DOY = WeatherStation$DOY, ET.instantaneous = TRUE,
                         ET= "ETr")

  Ts.datum <- Ts + (DEM - WeatherStation$elev) * 6.49 / 1000
  # hot and cold cells, checked before iterating
  anchor.cells <- anchorCells(anchors, Ts = Ts, Ts.datum = Ts.datum, Rn = Rn,
                              G = G, Z.om = Z.om, DEM = DEM, method = method)
  hot <- anchor.cells$hot
  cold <- anchor.cells$cold
  P <- 101.3*((293-0.0065 * DEM)/293)^5.26
  air.density <- 1000 * P / (1.01*(Ts)*287)
  latent.heat.vaporization <- (2.501-0.00236*(Ts-273.15))
  
  # initial conditions assuming neutral stability
  u.ws <- WeatherStation$wind * 0.41 / log(WeatherStation$height/Z.om.ws)
  u200.v <- u.ws / 0.41 * log(200/Z.om.ws)
  if(u200.v < 1){
    warning(paste0("u200 less than threshold value = ",
                   round(u200.v,4), "m/s. using u200 = 4m/s"))
    u200.v <- 4
  }
  
  u200 <- rast(DEM)
  terra::values(u200) <- u200.v
  if(mountainous==TRUE){
    u200 <- u200 * (1+0.1*((DEM-WeatherStation$elev)/1000))
  }
  
  friction.velocity <- 0.41 * u200 / log(200/Z.om)
  friction.velocity[friction.velocity==0] <- 0.1
  
  r.ah <- log(2/0.1)/(friction.velocity*0.41)
  
  LE.cold <- ETr.hourly * ETp.coef *
    (2.501 - 0.002361*(cell_mean(Ts, cold)-273.15)) * (1e6)/3600
  
  H.cold <- cell_mean(Rn, cold) - cell_mean(G, cold) - LE.cold
  
  result <- list()
  
  if(verbose==TRUE){
    print("starting conditions")
    print("Cold")
    print(data.frame(cbind(
      "Ts"       = cell_mean(Ts, cold),
      "Ts_datum" = cell_mean(Ts.datum, cold),
      "Rn"       = cell_mean(Rn, cold),
      "G"        = cell_mean(G, cold),
      "Z.om"     = cell_mean(Z.om, cold),
      "u200"     = cell_vals(u200, cold)[1],
      "u*"       = cell_mean(friction.velocity, cold)
    )))
    print("Hot")
    print(data.frame(cbind(
      "Ts"       = cell_mean(Ts, hot),
      "Ts_datum" = cell_mean(Ts.datum, hot),
      "Rn"       = cell_mean(Rn, hot),
      "G"        = cell_mean(G, hot),
      "Z.om"     = cell_mean(Z.om, hot),
      "u200"     = cell_vals(u200, hot)[1],
      "u*"       = cell_mean(friction.velocity, hot)
    )))
  }
  
  plot(1, cell_mean(r.ah, hot), xlim=c(0,15), ylim=c(0, cell_mean(r.ah, hot)),
       col="red", ylab="aerodynamic resistance s m-1", xlab="iteration", pch=20)
  graphics::points(1, cell_mean(r.ah, cold), col="blue", pch=20)
  
  converge <- FALSE
  last.loop <- FALSE
  i <- 1
  
  if(method == "mean"){
    
    while(!converge){
      
      if(verbose==TRUE){
        print(paste("iteration #", i))
      }
      i <-  i + 1
      
      dT.cold <- H.cold * cell_mean(r.ah, cold) / (cell_mean(air.density, cold)*1004)
      
      dT.hot <- (cell_mean(Rn, hot) - cell_mean(G, hot)) * cell_mean(r.ah, hot) /
        (cell_mean(air.density, hot)*1004)
      
      a <- (dT.hot - dT.cold) / (cell_mean(Ts.datum, hot) - cell_mean(Ts.datum, cold))
      b <- -a * cell_mean(Ts.datum, cold) + dT.cold
      
      if(verbose==TRUE){
        print(paste("a",a))
        print(paste("b",b))
      }
      
      dT <- as.numeric(a) * Ts.datum + as.numeric(b)
      rho <- 349.467*((((Ts-dT)-0.0065*DEM)/(Ts-dT))^5.26)/Ts
      H <- rho * 1004 * dT / r.ah
      
      Monin.Obukhov.L <- (air.density * -1004 * friction.velocity^3 * Ts) /
        (0.41 * 9.807 * H)
      
      phi.200 <- Monin.Obukhov.L * NA  # empty rast() breaks [<- when no cell matches
      phi.2   <- Monin.Obukhov.L * NA
      phi.01  <- Monin.Obukhov.L * NA
      
      # stable condition = L > 0
      phi.200[Monin.Obukhov.L > 0] <- (-5*(2/Monin.Obukhov.L))[Monin.Obukhov.L > 0]
      phi.2  [Monin.Obukhov.L > 0] <- (-5*(2/Monin.Obukhov.L))[Monin.Obukhov.L > 0]
      phi.01 [Monin.Obukhov.L > 0] <- (-5*(0.1/Monin.Obukhov.L))[Monin.Obukhov.L > 0]
      
      # unstable condition = L < 0
      x.200 <- (1- 16*(200/Monin.Obukhov.L))^0.25
      x.2   <- (1- 16*(2/Monin.Obukhov.L))^0.25
      x.01  <- (1- 16*(0.1/Monin.Obukhov.L))^0.25
      
      phi.200[Monin.Obukhov.L < 0] <- (2 * log((1+x.200)/2) +
                                         log((1 + x.200^2) /2) -
                                         2* atan(x.200) + 0.5 * pi)[Monin.Obukhov.L < 0]
      phi.2  [Monin.Obukhov.L < 0] <- (2 * log((1 + x.2^2) / 2))[Monin.Obukhov.L < 0]
      phi.01 [Monin.Obukhov.L < 0] <- (2 * log((1 + x.01^2) / 2))[Monin.Obukhov.L < 0]
      
      if(verbose==TRUE){
        print(paste("r.ah cold", cell_mean(r.ah, cold)))
        print(paste("r.ah hot",  cell_mean(r.ah, hot)))
        print(paste("dT cold",   cell_mean(dT, cold)))
        print(paste("dT hot",    cell_mean(dT, hot)))
        print("##############")
      }
      
      friction.velocity <- 0.41 * u200 / (log(200/Z.om) - phi.200)
      friction.velocity[friction.velocity==0] <- 0.1
      
      r.ah.hot.previous  <- cell_mean(r.ah, hot)
      r.ah.cold.previous <- cell_mean(r.ah, cold)
      
      r.ah <- (log(2/0.1) - phi.2 + phi.01) / (friction.velocity * 0.41)
      
      graphics::points(i, cell_mean(r.ah, hot),  col="red",  pch=20)
      graphics::points(i, cell_mean(r.ah, cold), col="blue", pch=20)
      lines(c(i, i-1), c(cell_mean(r.ah, hot),  r.ah.hot.previous),  col="red")
      lines(c(i, i-1), c(cell_mean(r.ah, cold), r.ah.cold.previous), col="blue")
      
      if(last.loop == TRUE){
        converge <- TRUE
        if(verbose==TRUE){print (paste0("convergence reached at iteration #", i))}
      }
      
      delta.r.ah.hot  <- (cell_mean(r.ah, hot)  - r.ah.hot.previous)  / cell_mean(r.ah, hot)  * 100
      delta.r.ah.cold <- (cell_mean(r.ah, cold) - r.ah.cold.previous) / cell_mean(r.ah, cold) * 100
      
      if(verbose==TRUE){
        print (paste("delta rah hot",  delta.r.ah.hot))
        print (paste("delta rah cold", delta.r.ah.cold))
        print ("### -------")
      }
      
      if(abs(delta.r.ah.hot) < 1 & abs(delta.r.ah.cold) < 1){last.loop <-  TRUE}
      if(i == maxit){
        warning(paste0("No convergence after ", i, " iterations: try different anchor values?"))
        break
      }
    }
    
  } else if(method=="lm"){
    
    converge <- FALSE
    last.loop <- FALSE
    i <- 1
    
    npairs <- min(c(length(cold), length(hot)))
    for(pair in 1:npairs){

      converge <- FALSE
      last.loop <- FALSE
      i <- 1  # maxit applies to each pair
      
      while(!converge){
        
        i <-  i + 1
        if(verbose==TRUE){
          print(paste("iteration #", i))
        }
        
        dT.cold <- H.cold * cell_vals(r.ah, cold[pair]) /
          (cell_vals(air.density, cold[pair]) * 1004)
        
        dT.hot <- (cell_vals(Rn, hot[pair]) - cell_vals(G, hot[pair])) *
          cell_vals(r.ah, hot[pair]) / (cell_vals(air.density, hot[pair]) * 1004)
        
        a <- (dT.hot - dT.cold) / (cell_vals(Ts.datum, hot[pair]) - cell_vals(Ts.datum, cold[pair]))
        b <- -a * cell_vals(Ts.datum, cold[pair]) + dT.cold
        
        if(verbose==TRUE){
          print(paste("a",a))
          print(paste("b",b))
        }
        
        dT <- as.numeric(a) * Ts.datum + as.numeric(b)
        rho <- 349.467*((((Ts-dT)-0.0065*DEM)/(Ts-dT))^5.26)/Ts
        H <- rho * 1004 * dT / r.ah
        
        Monin.Obukhov.L <- (air.density * -1004 * friction.velocity^3 * Ts) /
          (0.41 * 9.807 * H)
        
        phi.200 <- Monin.Obukhov.L * NA  # empty rast() breaks [<- when no cell matches
        phi.2   <- Monin.Obukhov.L * NA
        phi.01  <- Monin.Obukhov.L * NA
        
        phi.200[Monin.Obukhov.L > 0] <- (-5*(2/Monin.Obukhov.L))[Monin.Obukhov.L > 0]
        phi.2  [Monin.Obukhov.L > 0] <- (-5*(2/Monin.Obukhov.L))[Monin.Obukhov.L > 0]
        phi.01 [Monin.Obukhov.L > 0] <- (-5*(0.1/Monin.Obukhov.L))[Monin.Obukhov.L > 0]
        
        x.200 <- (1- 16*(200/Monin.Obukhov.L))^0.25
        x.2   <- (1- 16*(2/Monin.Obukhov.L))^0.25
        x.01  <- (1- 16*(0.1/Monin.Obukhov.L))^0.25
        
        phi.200[Monin.Obukhov.L < 0] <- (2 * log((1+x.200)/2) +
                                           log((1 + x.200^2) /2) -
                                           2* atan(x.200) + 0.5 * pi)[Monin.Obukhov.L < 0]
        phi.2  [Monin.Obukhov.L < 0] <- (2 * log((1 + x.2^2) / 2))[Monin.Obukhov.L < 0]
        phi.01 [Monin.Obukhov.L < 0] <- (2 * log((1 + x.01^2) / 2))[Monin.Obukhov.L < 0]
        
        if(verbose==TRUE){
          print(paste("r.ah cold", cell_vals(r.ah, cold[pair])))
          print(paste("r.ah hot",  cell_vals(r.ah, hot[pair])))
          print(paste("dT cold",   cell_vals(dT, cold[pair])))
          print(paste("dT hot",    cell_vals(dT, hot[pair])))
          print("##############")
        }
        
        friction.velocity <- 0.41 * u200 / (log(200/Z.om) - phi.200)
        friction.velocity[friction.velocity==0] <- 0.1
        
        r.ah.hot.previous  <- cell_vals(r.ah, hot[pair])
        r.ah.cold.previous <- cell_vals(r.ah, cold[pair])
        
        r.ah <- (log(2/0.1) - phi.2 + phi.01) / (friction.velocity * 0.41)
        
        graphics::points(i, cell_vals(r.ah, hot[pair]),  col="red",  pch=20)
        graphics::points(i, cell_vals(r.ah, cold[pair]), col="blue", pch=20)
        lines(c(i, i-1), c(cell_vals(r.ah, hot[pair]),  r.ah.hot.previous),  col="red")
        lines(c(i, i-1), c(cell_vals(r.ah, cold[pair]), r.ah.cold.previous), col="blue")
        
        if(last.loop == TRUE){
          converge <- TRUE
          if(verbose==TRUE){print (paste0("convergence reached at iteration #", i))}
        }
        
        delta.r.ah.hot  <- (cell_vals(r.ah, hot[pair])  - r.ah.hot.previous)  / cell_vals(r.ah, hot[pair])  * 100
        delta.r.ah.cold <- (cell_vals(r.ah, cold[pair]) - r.ah.cold.previous) / cell_vals(r.ah, cold[pair]) * 100
        
        if(verbose==TRUE){
          print (paste("delta rah hot",  delta.r.ah.hot))
          print (paste("delta rah cold", delta.r.ah.cold))
          print ("### -------")
        }
        
        if(abs(delta.r.ah.hot) < 1 & abs(delta.r.ah.cold) < 1){last.loop <-  TRUE}
        if(i >= maxit){
          warning(paste0("No convergence after ", i, " iterations for anchor pair ",
                         pair, ": try different anchor values?"))
          converge <- TRUE
        }
      }
    }
  }
  
  dT <- saveLoadClean(imagestack = dT, file = "dT", overwrite=TRUE)
  H  <- saveLoadClean(imagestack = H,  file = "H",  overwrite=TRUE)
  
  result$a  <- a
  result$b  <- b
  result$dT <- dT
  result$H  <- H
  return(result)
}
