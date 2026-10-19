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
#' @param fLAI            proportion of LAI lying above h/2: a number or a SpatRaster
#' on the grid of LAI. Only needed for method = "Perrier"
#' @param h               crop height in meters: a number or a SpatRaster on the
#' grid of LAI (e.g. rasterized from a field map; cells with 0 or NA, such as
#' bare soil, get the minimum Z.om of 0.005 m). Only needed for method = "Perrier"
#' @param mountainous     Logical. If TRUE applies the empirical adjustment for
#' the effects of general terrain roughness on momentum and heat transfer
#' (Allen et al., 2007); needs surface.model
#' @param surface.model   SpatRaster from METRICtopo() with a layer named Slope.
#' Only needed if mountainous = TRUE
#' @param slope.min1      logical. If TRUE, the terrain factor 1 + (slope - 5)
#' / 20 is at least 1 (METRIC manual), so slopes below 5 degrees do not lower
#' Z.om. See \code{\link{calcH.control}}
#' @details According to Allen et al. (2010), Zom is a measure of the form drag
#' and skin friction for the layer of air that interacts with the surface.
#' Values below 0.005 m (bare agricultural soil, Allen et al., 2007) are set
#' to 0.005 m, for every method.
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
#' @examples
#' # a small window of the Landsat 7 example scene, to run quickly
#' MTL <- system.file("extdata", "L7.MTL.txt", package = "water")
#' WeatherStation <- read.WSdata(WSdata = system.file("extdata", "apples.csv", package = "water"),
#'   date.format = "%d/%m/%Y", lat = -35.42222, long = -71.38639, elev = 201, height = 2.2,
#'   cf = c(1, 0.2777778, 1), tz = "America/Santiago", MTL = MTL)
#' image.DN <- terra::rast(L7_Talca)[1:60, 1:80, drop = FALSE]
#' DEM <- terra::rast(DEM_Talca)[1:60, 1:80, drop = FALSE]
#' surface.model <- METRICtopo(DEM)
#' solar.angles <- solarAngles(surface.model = surface.model,
#'                             WeatherStation = WeatherStation, MTL = MTL)
#' image.TOAr <- calcTOAr(image.DN, sat = "L7", MTL = MTL,
#'                        incidence.rel = solar.angles$incidence.rel)
#' LAI <- LAI(method = "metric2010", image = image.TOAr, L = 0.1)
#' Z.om <- momentumRoughnessLength(LAI = LAI, mountainous = TRUE, method = "short.crops",
#'                                 surface.model = surface.model)
#' plot(Z.om)
## Create a function to estimate a and b coefficients or the function between Z.om and NDVI
## using some points and tabulated z.om for their covers.
## Perrier by Santos 2012 and Pocas 2014.
momentumRoughnessLength <- function(method="short.crops", LAI, NDVI,
                                    albedo, a, b, fLAI, h,
                                    mountainous=FALSE, surface.model,
                                    slope.min1 = FALSE){
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
    # fLAI and h: one number, or a SpatRaster on the grid of LAI
    ras <- function(x) inherits(x, "SpatRaster")
    if(ras(fLAI) || ras(h)){
      for(x in list(fLAI, h)){
        if(ras(x) && !terra::compareGeom(x, LAI, stopOnError = FALSE)){
          stop("fLAI and h as rasters must have the grid of LAI", call. = FALSE)
        }
      }
    }
    if(ras(fLAI)){
      rg <- unlist(terra::global(fLAI, c("min", "max"), na.rm = TRUE))
      if(!all(is.finite(rg)) || rg[1] < 0 || rg[2] > 1){
        stop("fLAI must be between 0 and 1", call. = FALSE)
      }
    } else if(!is.numeric(fLAI) || length(fLAI) != 1 || !isTRUE(fLAI >= 0 && fLAI <= 1)){
      stop("fLAI must be a number (or a raster) between 0 and 1", call. = FALSE)
    }
    if(ras(h)){
      # 0 (or NA) cells, e.g. bare soil, end at the minimum roughness of 0.005 m
      if(!isTRUE(terra::global(h, "min", na.rm = TRUE)[1, 1] >= 0)){
        stop("h (crop height, m) cannot be negative", call. = FALSE)
      }
    } else if(!is.numeric(h) || length(h) != 1 || !isTRUE(h > 0)){
      stop("h (crop height, m) must be a number (or a raster) above 0", call. = FALSE)
    }
    a <- if(ras(fLAI)) terra::ifel(fLAI < 0.5, (2*(1-fLAI))^-1, 2*fLAI) else
      if(fLAI < 0.5) (2*(1-fLAI))^-1 else 2*fLAI
    Z.om <- ((1-exp(-a*LAI/2))*exp(-a*LAI/2))*h
  }
  if(mountainous==TRUE){
    terrain <- 1 + (180/pi*surface.model$Slope - 5)/20  # Allen et al. (2007), Eq. 35
    if(slope.min1){terrain <- terra::ifel(terrain < 1, 1, terrain)}
    Z.om <- Z.om * terrain
  }
  # minimum 0.005 m, bare agricultural soil (Allen et al., 2007, after Eq. 33)
  Z.om <- terra::ifel(Z.om < 0.005, 0.005, Z.om)
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
#' layers named R and NIR, used to compute NDVI.
#' Layers are selected by name; a stack without the standard names is used by
#' position, with a warning, and a stack with these names that lacks a needed
#' layer is an error.
#' @param Ts               land surface temperature in K. See \code{\link{surfaceTemperature}}.
#' Must be in Kelvin: values that look like Celsius, or more than 1\% of the
#' cells below 200 K (unmasked fill or clouds), stop; fewer cells below 200 K
#' give a warning.
#' @param LAI              SpatRaster with leaf area index, from LAI()
#' @param albedo           broadband surface albedo, from albedo()
#' @param Z.om             momentum roughness length, from momentumRoughnessLength()
#' @param n                number of pairs of anchor pixels to select (not used
#' by "percentile")
#' @param aoi              Not used (a warning is given when it is passed); kept for
#' backward compatibility
#' @param anchors.method   method for the selection of anchor pixels:
#' "percentile" (default), "flexible", "best" or "random". See Details.
#' @param WeatherStation   Optional. waterWeatherStation object (see
#' read.WSdata()). If supplied, the anchors are searched within WSbuffer of the
#' station
#' @param plots            Logical. If TRUE (default: only in interactive
#' sessions) the anchors are plotted over LAI: hot in red, cold in blue and the
#' weather station (if supplied) in black
#' @param deltaTemp        (not used by "percentile") temperature margin (K,
#' default 5): cold candidates must
#' be colder than the 5\% quantile of Ts over well-vegetated pixels plus
#' deltaTemp, and hot candidates warmer than the maximum Ts over bare-soil
#' pixels minus deltaTemp
#' @param minDist          (not used by "percentile") minimum distance (m)
#' between two anchor pixels of the same type
#' @param WSbuffer         maximum distance (m) to the weather station. If no
#' pixel is within WSbuffer, the whole image is used, with a warning
#' @param verbose          Logical. If TRUE, additional data are shown (as messages)
#' @param min.candidates   minimum number of hot and of cold candidates (default
#' 20) below which a warning is given: the area may be too small or too
#' homogeneous for the calibration
#' @param control          a list from \code{\link{anchors.control}}; only its
#' percentile parameters are read (cold.NDVI, cold.LAI, cold.fraction,
#' hot.NDVI, hot.fraction and homogeneity): deltaTemp, minDist, WSbuffer and
#' min.candidates are arguments of calcAnchors(). METRIC.EB() fills both from
#' its anchors.control
#' @details
#' anchors.method = "percentile" (default) uses criteria relative to the
#' scene, in the spirit of the automated calibration of Allen et al. (2013)
#' and Morton et al. (2013): cold candidates are full-cover vegetation, with
#' NDVI at or above its 95th percentile in the area and LAI of at least 4;
#' hot candidates are bare soil, with NDVI at or below its 10th percentile,
#' NDVI above 0 (no water) and albedo 0.10 to 0.35 (no urban surfaces or
#' bright rock). Candidates must have a homogeneous surface temperature (sd of
#' Ts in the 5 x 5 pixels around them up to 1.5 K). The coldest 20\% of the
#' cold candidates and the hottest 20\% of the hot candidates are the anchors,
#' and calcH() uses their means. The candidates are searched within WSbuffer
#' of the weather station, or in the whole image if there are none. The NDVI
#' quantiles are taken over land (NDVI > 0), so water in the area does not
#' move them. If no candidate is homogeneous, the homogeneity test is dropped;
#' if there are no cold candidates, cold.LAI is lowered to 3 and then
#' cold.NDVI to the 0.85 quantile; if there are no hot candidates, hot.NDVI is
#' raised to the 0.3 quantile. Attribute relaxation records what was relaxed.
#' The parameters are in \code{\link{anchors.control}}.
#'
#' For the other methods, cold candidates have LAI 3 to 6, albedo 0.18 to 0.25, Z.om 0.03 to 0.08, NDVI
#' within 0.15 of the maximum NDVI of the image and Ts below the threshold given
#' by deltaTemp. Hot candidates have albedo 0.13 to 0.15, NDVI 0.1 to 0.28,
#' Z.om up to 0.005 and Ts above the threshold given by deltaTemp.
#'
#' anchors.method = "best" selects the coldest and the hottest candidates;
#' "random" selects them at random among the candidates (CITRA-MCB method);
#' "flexible" works like "best" but, when there are not enough
#' candidates, relaxes the limits of the criteria (LAI, albedo and Z.om for cold
#' pixels; albedo, NDVI, Z.om and Ts for hot pixels) in 1\% steps up to 100\%,
#' first near the weather station (within WSbuffer) and then in the whole
#' image, with a warning. The old names "CITRA-MCB" and "CITRA-MCBbc" (now
#' "best") and "CITRA-MCBr" (now "random") are deprecated.
#' @return A data.frame with one row per anchor pixel and columns pixel (cell
#' number), X and Y (cell coordinates in the CRS of Ts), Ts (K), LAI, NDVI and
#' type ("hot" or "cold"). For anchors.method = "flexible" a column flex is
#' added: the relaxation applied to the criteria of each type (0 = none, 1 =
#' 100\%). Attribute candidates gives the number of cold and hot candidates
#' and, for "flexible" and "percentile", attribute relaxation the criteria
#' used. If plots =
#' TRUE the anchors are also plotted over LAI.
#' @author Guillermo Federico Olmedo
#' @author Daniel de la Fuente-Sáiz
#' @family sensible heat flux functions
#' @references
#' CITRA and MCB (personal communication)
#'
#' Allen, R. G., Burnett, B., Kramber, W., Huntington, J., Kjaersgaard, J.,
#' Kilic, A., Kelly, C., & Trezza, R. (2013). Automated calibration of the
#' METRIC-Landsat evapotranspiration process. Journal of the American Water
#' Resources Association, 49(3), 563-576. https://doi.org/10.1111/jawr.12056
#'
#' Morton, C. G., Huntington, J. L., Pohll, G. M., Allen, R. G., McGwire, K. C.,
#' & Bassett, S. D. (2013). Assessing calibration uncertainty and automation
#' for estimating evapotranspiration from agricultural areas using METRIC.
#' Journal of the American Water Resources Association, 49(3), 549-562.
#' https://doi.org/10.1111/jawr.12054
#' @export
#' @examples
#' # a small window of the Landsat 7 example scene, to run quickly
#' MTL <- system.file("extdata", "L7.MTL.txt", package = "water")
#' WeatherStation <- read.WSdata(WSdata = system.file("extdata", "apples.csv", package = "water"),
#'   date.format = "%d/%m/%Y", lat = -35.42222, long = -71.38639, elev = 201, height = 2.2,
#'   cf = c(1, 0.2777778, 1), tz = "America/Santiago", MTL = MTL)
#' image.DN <- terra::rast(L7_Talca)[1:60, 1:80, drop = FALSE]
#' DEM <- terra::rast(DEM_Talca)[1:60, 1:80, drop = FALSE]
#' surface.model <- METRICtopo(DEM)
#' solar.angles <- solarAngles(surface.model = surface.model,
#'                             WeatherStation = WeatherStation, MTL = MTL)
#' Rs.inc <- incSWradiation(surface.model = surface.model, solar.angles = solar.angles,
#'                          WeatherStation = WeatherStation)
#' image.TOAr <- calcTOAr(image.DN, sat = "L7", MTL = MTL,
#'                        incidence.rel = solar.angles$incidence.rel)
#' image.SR <- suppressWarnings(calcSR(image.TOAr, sat = "L7", surface.model = surface.model,
#'                    incidence.hor = solar.angles$incidence.hor, WeatherStation = WeatherStation))
#' albedo <- albedo(image.SR, coeff = "Tasumi", sat = "L7")
#' LAI <- LAI(method = "metric2010", image = image.TOAr, L = 0.1)
#' Ts <- suppressWarnings(surfaceTemperature(image.DN = image.DN, LAI = LAI, sat = "L7",
#'                                           WeatherStation = WeatherStation))
#' Rl.out <- outLWradiation(LAI = LAI, Ts = Ts)
#' Rl.inc <- incLWradiation(WeatherStation, DEM = surface.model$DEM,
#'                          solar.angles = solar.angles, Ts = Ts)
#' Rn <- netRadiation(LAI, albedo, Rs.inc, Rl.inc, Rl.out)
#' G <- soilHeatFlux(image = image.SR, Ts = Ts, albedo = albedo, Rn = Rn, LAI = LAI)
#' Z.om <- momentumRoughnessLength(LAI = LAI, mountainous = TRUE, method = "short.crops",
#'                                 surface.model = surface.model)
#' anchors <- calcAnchors(image = image.TOAr, Ts = Ts, LAI = LAI, albedo = albedo, Z.om = Z.om,
#'                       WeatherStation = WeatherStation, plots = FALSE, verbose = FALSE)
#' head(anchors)
calcAnchors  <- function(image, Ts, LAI, albedo, Z.om, n=1,
                         anchors.method= "percentile", WeatherStation,
                         plots=interactive(), deltaTemp=5, minDist = 500, WSbuffer = 30000,
                         verbose=FALSE, min.candidates = 20, control = NULL, aoi) {
  unusedAoi(missing(aoi))
  control <- controlList(control, "anchors.control")
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
  if(anchors.method=="percentile"){
    selected <- anchorsPercentile(NDVI = NDVI, LAI = LAI, Ts = Ts, albedo = albedo,
                                  WS.buffer = WS.buffer, control = control)
    hot <- selected$hot; cold <- selected$cold
    cold.n <- selected$cold.n; hot.n <- selected$hot.n
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
      values(NDVI>=max(values(NDVI), na.rm=TRUE)-0.15) &
      values(Z.om>=0.03) & values(Z.om<=0.08) &
      values(Ts<(minT+deltaTemp)) & values(WS.buffer == 1)
    hot.candidates <- values(albedo>=0.13) & values(albedo<=0.15) &
      values(NDVI>=0.1) & values(NDVI<=0.28) &
      values(Z.om<=0.005) & values(Ts>(maxT-deltaTemp)) & values(WS.buffer == 1)
    ### Test # anchors
    cold.n <- sum(as.numeric(cold.candidates), na.rm = TRUE)
    hot.n <- sum(as.numeric(hot.candidates), na.rm = TRUE)
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
          values(NDVI>=max(values(NDVI), na.rm=TRUE)-0.15) &
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
      values(NDVI>=max(values(NDVI), na.rm=TRUE)-0.15) &
      values(Z.om>=0.03) & values(Z.om<=0.08) &
      values(Ts<(minT+deltaTemp)) & values(WS.buffer == 1)
    hot.candidates <- values(albedo>=0.13) & values(albedo<=0.15) &
      values(NDVI>=0.1) & values(NDVI<=0.28) &
      values(Z.om<=0.005) & values(Ts>(maxT-deltaTemp)) & values(WS.buffer == 1)
    ### Test # anchors
    cold.n <- sum(as.numeric(cold.candidates), na.rm = TRUE)
    hot.n <- sum(as.numeric(hot.candidates), na.rm = TRUE)
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
          values(NDVI>=max(values(NDVI), na.rm=TRUE)-0.15) &
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
    NDVI.cold <- max(values(NDVI), na.rm=TRUE)-0.15
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
          if(verbose){message("Relaxing the criteria for ", type, " pixels: ", round(flex*100), "%")}
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
    cold.n <- sum(cold.candidates, na.rm = TRUE)
    hot.n <- sum(hot.candidates, na.rm = TRUE)
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
          values(NDVI>=max(values(NDVI), na.rm=TRUE)-0.15) &
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
  cells <- c(hot, cold)
  xy <- terra::xyFromCell(LAI, cells)
  at <- function(r){as.numeric(r[cells][[1]])}
  hot.and.cold <- data.frame(pixel = as.numeric(cells), X = xy[, 1], Y = xy[, 2],
                             Ts = at(Ts), LAI = at(LAI), NDVI = at(NDVI),
                             type = factor(rep(c("hot", "cold"), c(length(hot), length(cold))),
                                           levels = c("hot", "cold")))
  if(anchors.method=="flexible"){
    hot.and.cold$flex <- flex.cold
    hot.and.cold$flex[hot.and.cold$type == "hot"] <- flex.hot
    attr(hot.and.cold, "relaxation") <- list(
      cold = list(flex = flex.cold, criteria = optValCold,
                  wholeImage = length(coldSearch$inBuffer) == 1),
      hot = list(flex = flex.hot, criteria = optValHot,
                 wholeImage = length(hotSearch$inBuffer) == 1))
  }
  if(anchors.method=="percentile"){
    attr(hot.and.cold, "relaxation") <- selected$relaxation
  }
  attr(hot.and.cold, "candidates") <- c(cold = cold.n, hot = hot.n)
  if(min(cold.n, hot.n) < min.candidates){
    warning("Only ", cold.n, " cold and ", hot.n, " hot anchor candidates ",
            "(minimum ", min.candidates, "): the area may be too small or too ",
            "homogeneous; calibrate on an area of tens of km", call. = FALSE)
  }
  return(hot.and.cold)
}



# Anchors of anchors.method = "percentile", relative to the scene (in the
# spirit of the automated calibration of Allen et al., 2013, and Morton et
# al., 2013): cold candidates have NDVI >= its cold.NDVI quantile and LAI >=
# cold.LAI; hot candidates NDVI <= its hot.NDVI quantile, NDVI > 0 (no water)
# and albedo 0.10-0.35 (no urban surfaces or bright rock); both with a
# homogeneous surface temperature (sd of Ts in 5 x 5 pixels, as the 100 m
# thermal pixel spans about 3 x 3 Landsat pixels). The coldest cold.fraction
# of the cold candidates and the hottest hot.fraction of the hot candidates
# are the anchors. Near the weather station when possible; criteria relaxed
# when a set is empty
anchorsPercentile <- function(NDVI, LAI, Ts, albedo, WS.buffer, control){
  v <- data.frame(NDVI = values(NDVI)[, 1], LAI = values(LAI)[, 1],
                  Ts = values(Ts)[, 1], albedo = values(albedo)[, 1],
                  buffer = values(WS.buffer)[, 1] %in% 1)
  sdTs <- values(terra::focal(Ts, w = 5, fun = "sd", na.rm = TRUE))[, 1]
  valid <- stats::complete.cases(v[, c("NDVI", "LAI", "Ts", "albedo")])
  homogeneous <- valid & sdTs <= control$homogeneity & !is.na(sdTs)
  relax <- list(cold = list(), hot = list())
  # quantiles over land (NDVI > 0), so that water in the area does not move them
  land <- function(use){use & v$NDVI > 0}
  hotSet <- function(use, qhot){
    q <- stats::quantile(v$NDVI[land(use)], qhot, names = FALSE)
    land(use) & v$NDVI <= q & v$albedo >= 0.10 & v$albedo <= 0.35
  }
  pick <- function(use){
    q <- stats::quantile(v$NDVI[land(use)], control$cold.NDVI, names = FALSE)
    cold.ok <- land(use) & v$NDVI >= q & v$LAI >= control$cold.LAI
    list(cold = cold.ok, hot = hotSet(use, control$hot.NDVI))
  }
  for(inBuffer in c(TRUE, FALSE)){
    use <- valid & (if(inBuffer) v$buffer else TRUE)
    if(sum(land(use)) < 25){next}
    sets <- pick(use)
    # homogeneous candidates; if there are none, all candidates (with a note)
    for(type in c("cold", "hot")){
      h <- sets[[type]] & homogeneous
      if(any(h)){sets[[type]] <- h} else {relax[[type]]$homogeneity <- FALSE}
    }
    # cold: lower the LAI and NDVI limits when there are no candidates
    LAImin <- control$cold.LAI; qcold <- control$cold.NDVI
    while(!any(sets$cold) && (LAImin > 3 || qcold > 0.85)){
      if(LAImin > 3){LAImin <- LAImin - 0.25} else {qcold <- qcold - 0.01}
      q <- stats::quantile(v$NDVI[land(use)], qcold, names = FALSE)
      sets$cold <- land(use) & v$NDVI >= q & v$LAI >= LAImin
      relax$cold$cold.LAI <- LAImin; relax$cold$cold.NDVI <- qcold
    }
    # hot: raise the NDVI quantile when there are no candidates
    qhot <- control$hot.NDVI
    while(!any(sets$hot) && qhot < 0.3){
      qhot <- qhot + 0.01
      sets$hot <- hotSet(use, qhot)
      relax$hot$hot.NDVI <- qhot
    }
    if(any(sets$cold) && any(sets$hot)){
      relax$cold$wholeImage <- relax$hot$wholeImage <- !inBuffer
      if(!inBuffer && any(v$buffer)){
        warning("Not enough anchor candidates within WSbuffer of the weather ",
                "station; searched the whole image", call. = FALSE)
      }
      cold.T <- stats::quantile(v$Ts[sets$cold], control$cold.fraction, names = FALSE)
      hot.T <- stats::quantile(v$Ts[sets$hot], 1 - control$hot.fraction, names = FALSE)
      return(list(cold = which(sets$cold & v$Ts <= cold.T),
                  hot = which(sets$hot & v$Ts >= hot.T),
                  cold.n = sum(sets$cold), hot.n = sum(sets$hot), relaxation = relax))
    }
  }
  stop("Automatic selection of anchors FAILED: no cold (full-cover vegetation) or ",
       "no hot (dry bare soil) candidates. ", anchorHint(albedo), call. = FALSE)
}


# Air density (kg m-3) at the temperature of the air, Ts - dT, and elevation
# DEM (Allen et al., 2007, Eq. 37): 1000 P / (1.01 (Ts - dT) 287), with the
# pressure P (kPa) of the standard atmosphere for that temperature
airDensity <- function(Ts, dT, DEM){
  Ta <- Ts - dT
  349.467*(((Ta-0.0065*DEM)/Ta)^5.26)/Ta
}

# One iteration of the stability correction of calcH() for the pixels in p
# (vectors Ts, Ts.datum, DEM, Z.om, u200, air.density, r.ah and fv): dT and
# H with the current aerodynamic resistance r.ah and friction velocity fv,
# then the new fv and r.ah (Allen et al., 2007)
rahStep <- function(p, a, b, damping = 0){
  dT <- a * p$Ts.datum + b
  rho <- airDensity(p$Ts, dT, p$DEM)
  H <- rho * 1004 * dT / p$r.ah
  L <- (p$air.density * -1004 * p$fv^3 * p$Ts) / (0.41 * 9.807 * H)
  phi.200 <- phi.2 <- phi.01 <- rep(NA_real_, length(L))
  stable <- !is.na(L) & L > 0
  phi.200[stable] <- (-5*(2/L))[stable]
  phi.2[stable] <- (-5*(2/L))[stable]
  phi.01[stable] <- (-5*(0.1/L))[stable]
  x.200 <- (1- 16*(200/L))^0.25
  x.2   <- (1- 16*(2/L))^0.25
  x.01  <- (1- 16*(0.1/L))^0.25
  unstable <- !is.na(L) & L < 0
  phi.200[unstable] <- (2 * log((1+x.200)/2) + log((1 + x.200^2) /2) -
                          2* atan(x.200) + 0.5 * pi)[unstable]
  phi.2[unstable] <- (2 * log((1 + x.2^2) / 2))[unstable]
  phi.01[unstable] <- (2 * log((1 + x.01^2) / 2))[unstable]
  fv <- 0.41 * p$u200 / (log(200/p$Z.om) - phi.200)
  fv[!is.na(fv) & fv == 0] <- 0.1
  fv[!is.na(fv) & fv < 0] <- NA  # no valid stability correction
  r.new <- (log(2/0.1) - phi.2 + phi.01) / (fv * 0.41)
  r.ah <- if(damping > 0) (1 - damping) * r.new + damping * p$r.ah else r.new
  list(dT = dT, H = H, fv = fv, r.ah = r.ah)
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
#' @param Ts             Land surface temperature in K. See \code{\link{surfaceTemperature}}.
#' Must be in Kelvin: values that look like Celsius, or more than 1\% of the
#' cells below 200 K (unmasked fill or clouds), stop; fewer cells below 200 K
#' give a warning.
#' @param Z.om           momentum roughness length, from momentumRoughnessLength()
#' (which keeps it at least 0.005 m; calcH() uses it as given)
#' @param WeatherStation WeatherStation data at the satellite overpass.
#' Can be a waterWeatherStation object calculated using read.WSdata and MTL file
#' (with read.WSdata(..., overpass.radiation = "clear-sky"), the clear-sky
#' radiation is used)
#' @param ETp.coef       ratio of the ET of the cold pixels to the instantaneous
#' alfalfa reference ET (ETr) at the satellite overpass (see Details). Allen et al.
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
#' @param verbose        Logical. If TRUE, information about every iteration is
#' shown (as messages)
#' @param maxit          Maximum number of iterations. Default 20.
#' @param u200.min       minimum wind speed at 200 m (m/s, default 2). See
#' \code{\link{calcH.control}}
#' @param damping        weight of the previous aerodynamic resistance in each
#' iteration (default 0, no damping). See \code{\link{calcH.control}}
#' @param ETrF.hot       reference ET fraction of the hot anchors (default 0).
#' See \code{\link{calcH.control}}
#' @param plots          Logical. If TRUE (default: only in interactive sessions)
#' the aerodynamic resistance of the anchors is plotted at each iteration
#' @param lapse          lapse rate (K/km) used to refer Ts to the elevation of
#' the weather station (default 6.49). See \code{\link{calcH.control}}
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
#' The ET of the cold pixels is computed as ETp.coef times the instantaneous alfalfa
#' reference ET (ETr, ASCE standardized equation; see \code{\link{hourlyET}}) at the
#' satellite overpass, as in Allen et al. (2007). ET24h() and dailyET() use the
#' same reference.
#' @return A list with a and b (coefficients of dT = a * Ts_datum + b), dT
#' (SpatRaster, near-surface temperature difference, K), H (SpatRaster,
#' sensible heat flux, W/m2), converged (logical), iterations (number of
#' iterations; one per pair for method = "lm") and u200 (wind speed at 200 m
#' used, m/s). The aerodynamic resistance of the anchor pixels at each
#' iteration is plotted with plots = TRUE. The iteration stops with an error when
#' it diverges.
#' @author Guillermo Federico Olmedo
#' @author Daniel de la Fuente-Sáiz
#' @author Fernando Fuentes-Peñailillo
#' @family sensible heat flux functions
#' @references
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#'
#' Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration: guidelines for computing crop water requirements. FAO Irrigation and Drainage Paper 56. FAO, Rome.
#'
#' Allen, R., Irmak, A., Trezza, R., Hendrickx, J.M.H., Bastiaanssen, W., Kjaersgaard, J., 2011. Satellite-based ET estimation in agriculture using SEBAL and METRIC. Hydrol. Process. 25, 4011-4027. doi:10.1002/hyp.8408 \cr
#' @export
#' @examples
#' # a small window of the Landsat 7 example scene, to run quickly
#' MTL <- system.file("extdata", "L7.MTL.txt", package = "water")
#' WeatherStation <- read.WSdata(WSdata = system.file("extdata", "apples.csv", package = "water"),
#'   date.format = "%d/%m/%Y", lat = -35.42222, long = -71.38639, elev = 201, height = 2.2,
#'   cf = c(1, 0.2777778, 1), tz = "America/Santiago", MTL = MTL)
#' image.DN <- terra::rast(L7_Talca)[1:60, 1:80, drop = FALSE]
#' DEM <- terra::rast(DEM_Talca)[1:60, 1:80, drop = FALSE]
#' surface.model <- METRICtopo(DEM)
#' solar.angles <- solarAngles(surface.model = surface.model,
#'                             WeatherStation = WeatherStation, MTL = MTL)
#' Rs.inc <- incSWradiation(surface.model = surface.model, solar.angles = solar.angles,
#'                          WeatherStation = WeatherStation)
#' image.TOAr <- calcTOAr(image.DN, sat = "L7", MTL = MTL,
#'                        incidence.rel = solar.angles$incidence.rel)
#' image.SR <- suppressWarnings(calcSR(image.TOAr, sat = "L7", surface.model = surface.model,
#'                    incidence.hor = solar.angles$incidence.hor, WeatherStation = WeatherStation))
#' albedo <- albedo(image.SR, coeff = "Tasumi", sat = "L7")
#' LAI <- LAI(method = "metric2010", image = image.TOAr, L = 0.1)
#' Ts <- suppressWarnings(surfaceTemperature(image.DN = image.DN, LAI = LAI, sat = "L7",
#'                                           WeatherStation = WeatherStation))
#' Rl.out <- outLWradiation(LAI = LAI, Ts = Ts)
#' Rl.inc <- incLWradiation(WeatherStation, DEM = surface.model$DEM,
#'                          solar.angles = solar.angles, Ts = Ts)
#' Rn <- netRadiation(LAI, albedo, Rs.inc, Rl.inc, Rl.out)
#' G <- soilHeatFlux(image = image.SR, Ts = Ts, albedo = albedo, Rn = Rn, LAI = LAI)
#' Z.om <- momentumRoughnessLength(LAI = LAI, mountainous = TRUE, method = "short.crops",
#'                                 surface.model = surface.model)
#' anchors <- calcAnchors(image = image.TOAr, Ts = Ts, LAI = LAI, albedo = albedo, Z.om = Z.om,
#'                       WeatherStation = WeatherStation, plots = FALSE, verbose = FALSE)
#' H <- calcH(anchors = anchors, Ts = Ts, Z.om = Z.om, WeatherStation = WeatherStation,
#'            DEM = DEM, Rn = Rn, G = G, plots = FALSE, verbose = FALSE)
#' plot(H$H)
calcH  <- function(anchors, Ts, Z.om, WeatherStation, ETp.coef= 1.05,
                   Z.om.ws=0.015, mountainous=FALSE,
                   DEM, Rn, G, verbose=FALSE, maxit = 20, u200.min = 2,
                   damping = 0, ETrF.hot = 0, lapse = 6.49, plots = interactive(),
                   method = "mean") {
  # the method was the second argument: calcH(anchors, "mean", Ts, ...)
  if(!missing(Ts) && is.character(Ts)){
    stop("calcH(): method is now after the data arguments; name the arguments ",
         "(e.g. method = \"", Ts, "\")", call. = FALSE)
  }
  checkChoice(method, c("mean", "lm"), "method")
  control <- calcH.control(maxit = maxit, u200.min = u200.min, damping = damping,
                           ETrF.hot = ETrF.hot)

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

  Ts.datum <- Ts + (DEM - WeatherStation$elev) * lapse / 1000
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
  if(u200.v < u200.min){
    warning("Wind at the station at the overpass is ", round(WeatherStation$wind, 2),
            " m/s, so the wind extrapolated to 200 m (", round(u200.v, 2), " m/s) ",
            "is below the minimum of ", u200.min, " m/s, which is used. Check the ",
            "wind units (cf in read.WSdata()) and the sensor: results are ",
            "sensitive to this value", call. = FALSE)
    u200.v <- u200.min
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
  # ET of the hot pixels, 0 for dry bare soil (Allen et al., 2007, Eq. 47)
  LE.hot <- function(cells){
    ETr.hourly * ETrF.hot * (2.501 - 0.002361*(cell_vals(Ts, cells)-273.15)) * (1e6)/3600
  }
  
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
  
  if(plots){
    plot(1, cell_mean(r.ah, hot), xlim=c(0, maxit + 1), ylim=c(0, cell_mean(r.ah, hot)),
         col="red", ylab="aerodynamic resistance s m-1", xlab="iteration", pch=20)
    graphics::points(1, cell_mean(r.ah, cold), col="blue", pch=20)
  }
  
  # The iteration of each pixel depends only on its own values and on a and b,
  # which depend only on the anchors: iterate on the anchor pixels, then
  # replay the same sequence of a and b once for all the pixels
  anchor <- unique(c(hot, cold))
  pix <- function(r){cell_vals(r, anchor)}
  state <- list(Ts = pix(Ts), Ts.datum = pix(Ts.datum), DEM = pix(DEM), Z.om = pix(Z.om),
                dT = rep(0, length(anchor)),
                u200 = pix(u200), air.density = pix(air.density),
                r.ah = pix(r.ah), fv = pix(friction.velocity))
  at <- function(x, cells){x[match(cells, anchor)]}
  Rn.a <- cell_vals(Rn, anchor); G.a <- cell_vals(G, anchor)
  ab <- list()
  iterations <- integer()
  converged <- TRUE
  pairs <- if(method == "mean") list(list(hot = hot, cold = cold)) else {
    npairs <- min(c(length(cold), length(hot)))
    if(npairs > 10){
      warning("method = \"lm\" iterates each of the ", npairs, " pairs of anchors in ",
              "turn and uses the last one: slow and no better than \"mean\" with many ",
              "anchors (e.g. from anchors.method = \"percentile\")", call. = FALSE)
    }
    lapply(1:npairs, function(k) list(hot = hot[k], cold = cold[k]))
  }
  step <- 1
  for(pair in seq_along(pairs)){
    h <- pairs[[pair]]$hot; k <- pairs[[pair]]$cold
    converge <- FALSE
    last.loop <- FALSE
    i <- 1  # maxit applies to each pair
    while(!converge){
      if(verbose==TRUE){
        print(paste("iteration #", i))
      }
      i <-  i + 1
      step <- step + 1
      # density at the temperature of the air, with dT of the previous iteration
      rhoA <- airDensity(state$Ts, state$dT, state$DEM)
      dT.cold <- H.cold * mean(at(state$r.ah, k), na.rm = TRUE) /
        (mean(at(rhoA, k), na.rm = TRUE)*1004)
      dT.hot <- (mean(at(Rn.a, h), na.rm = TRUE) - mean(at(G.a, h), na.rm = TRUE) -
                   mean(LE.hot(h))) *
        mean(at(state$r.ah, h), na.rm = TRUE) / (mean(at(rhoA, h), na.rm = TRUE)*1004)
      a <- (dT.hot - dT.cold) / (mean(at(state$Ts.datum, h), na.rm = TRUE) -
                                   mean(at(state$Ts.datum, k), na.rm = TRUE))
      b <- -a * mean(at(state$Ts.datum, k), na.rm = TRUE) + dT.cold
      checkIteration(a, b)
      if(verbose==TRUE){
        print(paste("a",a))
        print(paste("b",b))
      }
      ab[[length(ab) + 1]] <- c(a, b)
      r.ah.hot.previous  <- mean(at(state$r.ah, h), na.rm = TRUE)
      r.ah.cold.previous <- mean(at(state$r.ah, k), na.rm = TRUE)
      new <- rahStep(state, a, b, damping)
      if(verbose==TRUE){
        print(paste("r.ah cold", r.ah.cold.previous))
        print(paste("r.ah hot",  r.ah.hot.previous))
        print(paste("dT cold",   mean(at(new$dT, k), na.rm = TRUE)))
        print(paste("dT hot",    mean(at(new$dT, h), na.rm = TRUE)))
        print("##############")
      }
      checkIteration(u = at(new$fv, c(h, k)))
      state$r.ah <- new$r.ah; state$fv <- new$fv; state$dT <- new$dT
      r.ah.hot  <- mean(at(state$r.ah, h), na.rm = TRUE)
      r.ah.cold <- mean(at(state$r.ah, k), na.rm = TRUE)
      if(plots){
        graphics::points(step, r.ah.hot,  col="red",  pch=20)
        graphics::points(step, r.ah.cold, col="blue", pch=20)
        lines(c(step, step-1), c(r.ah.hot,  r.ah.hot.previous),  col="red")
        lines(c(step, step-1), c(r.ah.cold, r.ah.cold.previous), col="blue")
      }
      if(last.loop == TRUE){
        converge <- TRUE
        if(verbose==TRUE){print (paste0("convergence reached at iteration #", i))}
      }
      delta.r.ah.hot  <- (r.ah.hot  - r.ah.hot.previous)  / r.ah.hot  * 100
      delta.r.ah.cold <- (r.ah.cold - r.ah.cold.previous) / r.ah.cold * 100
      if(verbose==TRUE){
        print (paste("delta rah hot",  delta.r.ah.hot))
        print (paste("delta rah cold", delta.r.ah.cold))
        print ("### -------")
      }
      if(abs(delta.r.ah.hot) < 1 & abs(delta.r.ah.cold) < 1){last.loop <-  TRUE}
      if(!converge && i - 1 >= maxit){
        warning("No convergence after ", i - 1, " iterations",
                if(method == "lm") paste0(" for anchor pair ", pair),
                ": try calcH.control(damping = 0.5), more iterations or other anchors",
                call. = FALSE)
        converged <- FALSE
        break
      }
    }
    iterations <- c(iterations, i - 1)
  }
  # all the pixels: the same iterations, in one pass by blocks
  replay <- function(Ts, Ts.datum, DEM, Z.om, u200, air.density){
    fv <- 0.41 * u200 / log(200/Z.om)
    fv[!is.na(fv) & fv == 0] <- 0.1
    p <- list(Ts = Ts, Ts.datum = Ts.datum, DEM = DEM, Z.om = Z.om, u200 = u200,
              air.density = air.density, r.ah = log(2/0.1)/(fv*0.41), fv = fv)
    for(j in seq_along(ab)){
      new <- rahStep(p, ab[[j]][1], ab[[j]][2], damping)
      p$r.ah <- new$r.ah; p$fv <- new$fv
    }
    cbind(dT = new$dT, H = new$H)
  }
  inputs <- c(Ts, Ts.datum, DEM, Z.om, u200, air.density)
  names(inputs) <- c("Ts", "Ts.datum", "DEM", "Z.om", "u200", "air.density")
  out <- terra::lapp(inputs, replay, usenames = TRUE)
  dT <- out[[1]]
  H <- out[[2]]

  dT <- saveLoadClean(imagestack = dT, file = "dT", overwrite=TRUE)
  H  <- saveLoadClean(imagestack = H,  file = "H",  overwrite=TRUE)
  
  result$a  <- a
  result$b  <- b
  result$dT <- dT
  result$H  <- H
  result$converged <- converged
  result$iterations <- iterations
  result$u200 <- u200.v
  return(result)
}
