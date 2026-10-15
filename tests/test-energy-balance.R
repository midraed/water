# Regression tests for the energy balance: topography, solar angles, TOA
# reflectance, LAI, anchor pixels and the METRIC wrappers. Plain R, run by
# R CMD check.
library(water)
pdf(NULL)

expect_warning <- function(expr, pattern) {
  w <- NULL
  value <- withCallingHandlers(expr, warning = function(c) {
    if (grepl(pattern, conditionMessage(c))) w <<- conditionMessage(c)
    invokeRestart("muffleWarning")
  })
  if (is.null(w)) stop("expected a warning matching: ", pattern)
  invisible(value)
}
quiet <- function(expr) {
  value <- NULL
  invisible(utils::capture.output(value <- suppressWarnings(expr)))
  value
}
maxdiff <- function(x, y) terra::global(abs(x - y), "max", na.rm = TRUE)[1, 1]

extdata <- system.file("extdata", package = "water")
c2 <- file.path(extdata, "L8_C2")
MTL.L1 <- file.path(extdata, "LC82320832016040LGN00_MTL.txt")
MTL.L2 <- list.files(c2, pattern = "MTL.txt$", full.names = TRUE)
WS <- quiet(read.WSdata(WSdata = file.path(extdata, "INTA.csv"),
                        datetime.format = "%Y/%m/%d %H:%M",
                        columns = c("date" = 1, "time" = 1, "radiation" = 5,
                                    "wind" = 6, "RH" = 3, "temp" = 2, "rain" = 4),
                        lat = -33.00513, long = -68.86469, elev = 927, height = 2,
                        MTL = MTL.L2, tz = "America/Argentina/Mendoza"))
flatDEM <- function(image) {
  DEM <- terra::rast(image[[1]]); terra::values(DEM) <- 927; DEM
}

## METRICtopo(): slope in radians (terra::terrain() defaults to degrees);
## a constant DEM is flat everywhere, without an NA border
topo <- METRICtopo(terra::rast(DEM_Talca))
stopifnot(terra::global(topo$Slope, "max", na.rm = TRUE)[1, 1] < pi/2)
image.DN <- loadImage(path = extdata, sat = "L8")
flat <- METRICtopo(flatDEM(image.DN))
stopifnot(terra::global(is.na(flat$Slope), "sum")[1, 1] == 0,
          terra::global(flat$Slope, "max")[1, 1] == 0)

## solarAngles(): the hour angle is negative before solar noon
angles <- solarAngles(surface.model = flat, WeatherStation = WS, MTL = MTL.L1)
stopifnot(terra::global(angles$hour.angle, "max")[1, 1] < 0)

## calcTOAr(): L8 reflectance divided by cos(incidence angle), not the angle
toa <- calcTOAr(image.DN = image.DN, sat = "L8", MTL = MTL.L1,
                incidence.rel = angles$incidence.rel)
sun <- as.numeric(sub(".*= ", "", grep("SUN_ELEVATION", readLines(MTL.L1), value = TRUE)))
usgs <- (2e-5 * image.DN$NIR - 0.1) / sin(sun * pi / 180)
stopifnot(abs(terra::global(toa$NIR - usgs, "mean", na.rm = TRUE)[1, 1]) < 0.005)

## calcTOAr(): L7 bands by name (L7_Talca stores SWIR2 after Thermal1)
L7 <- terra::rast(L7_Talca)
MTL.L7 <- file.path(extdata, "L7.MTL.txt")
WS7 <- quiet(read.WSdata(WSdata = file.path(extdata, "apples.csv"),
                         date.format = "%d/%m/%Y", lat = -35.42222, long = -71.38639,
                         elev = 201, height = 2.2, tz = "America/Santiago",
                         columns = c(1, 2, 3, 4, 6, 7, 8), MTL = MTL.L7))
angles7 <- solarAngles(surface.model = METRICtopo(flatDEM(L7)), WeatherStation = WS7,
                       MTL = MTL.L7)
toa7 <- calcTOAr(image.DN = L7, sat = "L7", MTL = MTL.L7,
                 incidence.rel = angles7$incidence.rel)
stopifnot(stats::median(terra::values(toa7$SWIR2), na.rm = TRUE) < 0.2)

## LAI(method = "metric"): LAI = 6 for dense vegetation (it was NaN)
dense <- terra::rast(nrows = 1, ncols = 1, nlyrs = 6,
                     vals = c(0.03, 0.05, 0.05, 0.4786, 0.2, 0.1))
names(dense) <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2")
stopifnot(as.numeric(terra::values(LAI(method = "metric", image = dense))) == 6)

## calcAnchors(): masked reflectance pixels are never anchors, and plots
## work without a weather station
sr <- loadImageSR(c2)
Ts <- sr$Band_10
cloud <- terra::cellFromRowCol(sr, rep(60:70, each = 11), rep(100:110, 11))
sr[cloud] <- NA
Ts[cloud] <- 285  # cold cloud tops, not masked in the thermal band
alb <- albedo(sr, coeff = "Olmedo", sat = "L8")
lai <- LAI(method = "metric2010", image = sr)
zom <- momentumRoughnessLength(LAI = lai)
anchors <- quiet(calcAnchors(image = sr, Ts = Ts, LAI = lai, albedo = alb, Z.om = zom,
                             n = 1, plots = TRUE, anchors.method = "flexible"))
stopifnot(!any(anchors$pixel %in% cloud))

## "random" with a single candidate picks that candidate
stopifnot(all(replicate(20, water:::sampleOne(10585L)) == 10585L))

## METRIC.EB() with Collection 2 Level-2 data, cloud mask, print() and writeRaster()
image.SR <- cfmask(path = c2, image = loadImageSR(c2))
EB <- quiet(METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10,
                      WeatherStation = WS, MTL = MTL.L2, alb.coeff = "Olmedo",
                      n = 5, Z.om.ws = 0.03, ETp.coef = 1.2))
stopifnot(inherits(EB, "waterLSEB"), nrow(EB$anchors) == 10)
invisible(utils::capture.output(print(EB)))
f <- tempfile(fileext = ".tif")
writeRaster(EB, f)
stopifnot(terra::nlyr(terra::rast(f)) == 5)

## METRIC.EB(): weather station outside the image (#46)
WSfar <- WS
WSfar$location$lat <- WS$location$lat + 0.5  # about 55 km north
EBfar <- quiet(expect_warning(
  METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10, WeatherStation = WSfar,
            MTL = MTL.L2, alb.coeff = "Olmedo", n = 1, Z.om.ws = 0.03),
  "No pixel is within WSbuffer"))
stopifnot(inherits(EBfar, "waterLSEB"))

## METRIC.EB(): manual anchors with the default anchors.method, and
## Zom.method = "custom"
EBman <- quiet(METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10,
                         WeatherStation = WS, MTL = MTL.L2, alb.coeff = "Olmedo",
                         anchors = EB$anchors, Z.om.ws = 0.03, ETp.coef = 1.2))
stopifnot(maxdiff(EBman$EB$LatentHeat, EB$EB$LatentHeat) < 1e-6)
EBzom <- quiet(METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10,
                         WeatherStation = WS, MTL = MTL.L2, alb.coeff = "Olmedo",
                         Zom.method = "custom", extraParameters = c(a = 0.1, b = -5),
                         Z.om.ws = 0.03))
stopifnot(inherits(EBzom, "waterLSEB"))

## METRIC.Rn() and METRIC.G() give the same Rn and G as METRIC.EB(), also
## with split window Ts and another L
args <- list(image.DN = image.DN, image.SR = loadImageSR(c2), WeatherStation = WS,
             MTL = MTL.L1, sat = "L8", alb.coeff = "Olmedo", LST.method = "SW",
             L = 0.2)
EBsw <- quiet(do.call(METRIC.EB, c(args, n = 1, Z.om.ws = 0.03)))
Rn <- quiet(do.call(METRIC.Rn, args))
G <- quiet(do.call(METRIC.G, c(args, list(Rn = Rn))))
stopifnot(maxdiff(Rn, EBsw$EB$NetRadiation) < 1e-8,
          maxdiff(G, EBsw$EB$SoilHeat) < 1e-8)

cat("test-energy-balance.R: all tests passed\n")
