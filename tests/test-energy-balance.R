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
expect_error <- function(expr, pattern) {
  msg <- tryCatch({ expr; NULL }, error = function(e) conditionMessage(e))
  if (is.null(msg)) stop("expected an error matching: ", pattern)
  if (!grepl(pattern, msg)) stop("unexpected error: ", msg)
  invisible(msg)
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

## Kt (clearness of the atmosphere) is an argument (#84); it was fixed at 1
sw1 <- incSWradiation(flat, angles, WS)
stopifnot(maxdiff(incSWradiation(flat, angles, WS, Kt = 1), sw1) == 0,
          terra::global(incSWradiation(flat, angles, WS, Kt = 0.5) - sw1, "max")[1, 1] < 0)

## Weather at the overpass (#84): radiation against clear sky, period means
## and the clear-sky option. INTA: 587 W/m2 against a clear-sky 860 W/m2
ws.args <- list(WSdata = file.path(extdata, "INTA.csv"), datetime.format = "%Y/%m/%d %H:%M",
                columns = c("date" = 1, "time" = 1, "radiation" = 5, "wind" = 6, "RH" = 3,
                            "temp" = 2, "rain" = 4),
                lat = -33.00513, long = -68.86469, elev = 927, height = 2, MTL = MTL.L2,
                tz = "America/Argentina/Mendoza")
o <- water:::getDataWS(WS)
Rso <- water:::clearSkyRso(o$DOY, o$hours + 0.5, o$lat, o$long, o$long.z, o$elev) * 1e6 / 3600
stopifnot(abs(Rso - 860) < 1)
invisible(expect_warning(do.call(read.WSdata, ws.args), "68% of the clear-sky value"))
WSend <- quiet(do.call(read.WSdata, c(ws.args, timestamps = "end")))
stopifnot(abs(WSend$at.sat$radiation - 638) < 2)
WSclear <- quiet(do.call(read.WSdata, c(ws.args, overpass.radiation = "clear-sky")))
stopifnot(abs(water:::getDataWS(WSclear)$radiation - 860) < 1,
          WSclear$at.sat$radiation.measured == WS$at.sat$radiation,
          abs(hourlyET(WSclear) - 0.679) < 0.005)

## calcTOAr(): L8 reflectance divided by cos(incidence angle), not the angle
toa <- calcTOAr(image.DN = image.DN, sat = "L8", MTL = MTL.L1,
                incidence.rel = angles$incidence.rel)
sun <- as.numeric(sub(".*= ", "", grep("SUN_ELEVATION", readLines(MTL.L1), value = TRUE)))
usgs <- (2e-5 * image.DN$NIR - 0.1) / sin(sun * pi / 180)
stopifnot(abs(terra::global(toa$NIR - usgs, "mean", na.rm = TRUE)[1, 1]) < 0.005)

## calcTOAr(): L7 bands by name (L7_Talca stores SWIR2 after Thermal1)
L7 <- terra::rast(L7_Talca)
MTL.L7 <- file.path(extdata, "L7.MTL.txt")
w7 <- character()
WS7 <- withCallingHandlers(
  read.WSdata(WSdata = file.path(extdata, "apples.csv"),
              date.format = "%d/%m/%Y", lat = -35.42222, long = -71.38639,
              elev = 201, height = 2.2, tz = "America/Santiago",
              columns = c(1, 2, 3, 4, 6, 7, 8), cf = c(1, 0.2777778, 1),
              MTL = MTL.L7),
  warning = function(c) { w7 <<- c(w7, conditionMessage(c)); invokeRestart("muffleWarning") })
# a clear overpass: Rs/Rso is 0.94, no warning (#84)
stopifnot(!any(grepl("clear-sky", w7)))
angles7 <- solarAngles(surface.model = METRICtopo(flatDEM(L7)), WeatherStation = WS7,
                       MTL = MTL.L7)
toa7 <- calcTOAr(image.DN = L7, sat = "L7", MTL = MTL.L7,
                 incidence.rel = angles7$incidence.rel)
stopifnot(stats::median(terra::values(toa7$SWIR2), na.rm = TRUE) < 0.2)

## calcSR() detects the satellite from the MTL file (#78)
sr7 <- calcSR(image.TOAr = toa7, incidence.hor = angles7$incidence.hor,
              WeatherStation = WS7, surface.model = METRICtopo(flatDEM(L7)),
              MTL = MTL.L7)
stopifnot(terra::nlyr(sr7) == 6)

## LAI(method = "metric"): LAI = 6 for dense vegetation (it was NaN)
dense <- terra::rast(nrows = 1, ncols = 1, nlyrs = 6,
                     vals = c(0.03, 0.05, 0.05, 0.4786, 0.2, 0.1))
names(dense) <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2")
stopifnot(as.numeric(terra::values(LAI(method = "metric", image = dense))) == 6)

## LAI(method = "turner"): inverse of the NDVI-LAI relationship of Turner et al.
## (1999, Table 3); it was applied as LAI = f(NDVI) and never exceeded 0.75
turnerNDVI <- function(L) 0.5724 + 0.0989*L - 0.0114*L^2 + 0.0004*L^3
target <- c(0, 1, 3, 6)
red <- 0.05
nir <- red * (1 + turnerNDVI(target)) / (1 - turnerNDVI(target))
turner.img <- terra::rast(nrows = 1, ncols = 4, nlyrs = 6)
terra::values(turner.img) <- cbind(0.03, 0.04, red, nir, 0.2, 0.1)
names(turner.img) <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2")
stopifnot(isTRUE(all.equal(as.numeric(terra::values(LAI(method = "turner", image = turner.img))),
                           target, tolerance = 0.001)))

## Bands by name everywhere (#86): reordered layers give the same LAI and
## albedo (LAI dropped from 1.25 to 0.001); a missing band is an error, and
## layers without standard names are used by position, with a warning
sr0 <- suppressMessages(loadImageSR(c2))
rev6 <- sr0[[c("NIR", "R", "G", "B", "SWIR2", "SWIR1")]]
stopifnot(maxdiff(LAI(image = rev6), LAI(image = sr0)) == 0,
          maxdiff(albedo(rev6, coeff = "Olmedo", sat = "L8"),
                  albedo(sr0, coeff = "Olmedo", sat = "L8")) == 0,
          # a Landsat stack that starts with R and NIR is not MODIS
          maxdiff(albedo(sr0[[c("R", "NIR", "B", "G", "SWIR1", "SWIR2")]], coeff = "Olmedo"),
                  albedo(sr0, coeff = "Olmedo")) == 0)
noR <- sr0; names(noR)[3] <- "Red"
expect_error(LAI(image = noR), "no layer named \"R\"")
unnamed <- sr0[[1:6]]; names(unnamed) <- paste0("lyr", 1:6)
stopifnot(maxdiff(expect_warning(LAI(image = unnamed), "by position"), LAI(image = sr0)) == 0)
# MODIS stacks start with R and NIR: LAI used B and G
stopifnot(maxdiff(LAI(image = sr0[[c("R", "NIR", "B", "G", "SWIR1", "SWIR2")]]),
                  LAI(image = sr0)) == 0)

## Surface temperature in Kelvin (#86): Celsius stops (ET was +86%); a few
## cells below 200 K (cloud tops) warn, many stop
expect_error(quiet(METRIC.Rn(image.SR = sr0, surf.temp = sr0$Band_10 - 273.15,
                             WeatherStation = WS, MTL = MTL.L2, alb.coeff = "Olmedo")),
             "Celsius")
tops <- sr0$Band_10; tops[1:5] <- 190
invisible(expect_warning(water:::checkKelvin(tops), "5 cells .* below 200 K"))
tops[1:1000] <- 190
expect_error(water:::checkKelvin(tops), "below 200 K")

## calcAnchors(): masked reflectance pixels are never anchors, and plots
## work without a weather station
sr <- loadImageSR(c2)
Ts <- sr$Band_10
cloud <- terra::cellFromRowCol(sr, rep(60:70, each = 11), rep(100:110, 11))
sr[cloud] <- NA
Ts[cloud] <- 285  # cold cloud tops, not masked in the thermal band
alb <- albedo(sr, coeff = "Olmedo", sat = "L8")
# sat = "auto" no longer needs Landsat files in the working directory (#78)
stopifnot(maxdiff(albedo(sr, coeff = "Olmedo"), alb) == 0)
lai <- LAI(method = "metric2010", image = sr)
zom <- momentumRoughnessLength(LAI = lai)
anchors <- quiet(calcAnchors(image = sr, Ts = Ts, LAI = lai, albedo = alb, Z.om = zom,
                             n = 1, plots = TRUE, anchors.method = "flexible"))
stopifnot(!any(anchors$pixel %in% cloud))

## "random": candidates only, also with a single hot candidate (sample(x, 1)
## drew from 1:x)
syn <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 300, ymin = 0, ymax = 300,
                   crs = "EPSG:32719", nlyrs = 2, vals = 0)
names(syn) <- c("R", "NIR")
cold.cells <- 1:10; hot.cells <- 91:100
set.layer <- function(x, cold, hot, other) {
  v <- rep(other, terra::ncell(x)); v[cold.cells] <- cold; v[hot.cells] <- hot
  terra::rast(x[[1]], vals = v)
}
syn$R <- set.layer(syn, 0.03, 0.2, 0.15); syn$NIR <- set.layer(syn, 0.37, 0.3, 0.45)
Ts.syn <- set.layer(syn, 295, 320, 305); Ts.syn[95] <- 321
for (seed in 1:5) {
  set.seed(seed)
  a <- quiet(calcAnchors(image = syn, Ts = Ts.syn, LAI = set.layer(syn, 4, 0.1, 1),
                         albedo = set.layer(syn, 0.2, 0.14, 0.3),
                         Z.om = set.layer(syn, 0.05, 0.002, 0.01), n = 1,
                         plots = FALSE, anchors.method = "random"))
  stopifnot(a$pixel[a$type == "hot"] == 95, a$pixel[a$type == "cold"] %in% cold.cells)
}

## METRIC.EB() with Collection 2 Level-2 data, cloud mask, print() and writeRaster()
image.SR <- cfmask(path = c2, image = loadImageSR(c2))
EB <- quiet(METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10,
                      WeatherStation = WS, MTL = MTL.L2, alb.coeff = "Olmedo",
                      n = 5, Z.om.ws = 0.03))
stopifnot(inherits(EB, "waterLSEB"), nrow(EB$anchors) == 10)
invisible(utils::capture.output(print(EB)))
f <- tempfile(fileext = ".tif")
writeRaster(EB, f)
stopifnot(terra::nlyr(terra::rast(f)) == 5)

## Pinned results: means of the layers and anchor cells, so that changes in
## the results are deliberate (relative tolerance, for other platforms)
near <- function(x, y, tol = 1e-6) all(abs(x - y) <= tol * abs(y))
# (#85: surface reflectance outside 0-1 masked, mean LE +4.2%; #89: ASCE
# reference ET at the overpass, LE +0.01%)
pin.L8 <- c(576.7867019656, 89.5582897965, 232.5746736555, 259.7305664032, 311.2443506301)
stopifnot(near(terra::global(EB$EB, "mean", na.rm = TRUE)[, 1], pin.L8),
          setequal(terra::cellFromXY(EB$EB, sf::st_coordinates(EB$anchors)),
                   c(4517, 11142, 23763, 24554, 16283, 1166, 12315, 4794, 6342, 18343)))

## anchors.control and calcH.control (#102); the calcAnchors() parameters
## in extraParameters still work, with a message
EBctl <- quiet(METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10,
                         WeatherStation = WS, MTL = MTL.L2, alb.coeff = "Olmedo",
                         n = 5, Z.om.ws = 0.03, anchors.control = list(minDist = 500),
                         calcH.control = calcH.control(maxit = 20)))
EBold <- quiet(suppressMessages(METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10,
                                          WeatherStation = WS, MTL = MTL.L2,
                                          alb.coeff = "Olmedo", n = 5, Z.om.ws = 0.03,
                                          extraParameters = c(minDist = 500))))
stopifnot(maxdiff(EBctl$EB$LatentHeat, EB$EB$LatentHeat) == 0,
          maxdiff(EBold$EB$LatentHeat, EB$EB$LatentHeat) == 0)
expect_error(anchors.control(minDist = -1), "above 0")
expect_error(calcH.control(maxit = 0.5), "whole number")
expect_error(METRIC.EB(image.SR = "x", anchors.control = list(mindist = 1)),
             "Unknown elements in anchors.control")

## calcH(): the anchors are checked before iterating (#88)
xy <- sf::st_coordinates(EB$anchors)
df <- data.frame(X = xy[, 1], Y = xy[, 2], type = EB$anchors$type)
Ts.eb <- EB$EB$surfaceTemperature
H.args <- list(Ts = Ts.eb, Z.om = momentumRoughnessLength(LAI = LAI(image = image.SR)),
               WeatherStation = WS, DEM = flatDEM(image.SR), Rn = EB$EB$NetRadiation,
               G = EB$EB$SoilHeat, Z.om.ws = 0.03)
runH <- function(anchors, ...) {
  quiet(do.call(calcH, c(list(anchors = anchors), H.args, list(...))))
}
H0 <- runH(df)
stopifnot(maxdiff(H0$H, EB$EB$SensibleHeat) < 1e-6)
# hot-cold difference of the example, 1.6 K above the 5 K minimum
dT <- water:::anchorCells(df, Ts = Ts.eb, Ts.datum = Ts.eb, Rn = EB$EB$NetRadiation,
                          G = EB$EB$SoilHeat, Z.om = H.args$Z.om, DEM = H.args$DEM)$dT
stopifnot(near(dT, 6.671291436))
up <- df; up$type <- ifelse(df$type == "hot", "Hot", " COLD")
stopifnot(maxdiff(runH(up)$H, H0$H) == 0)
ll <- sf::st_transform(sf::st_as_sf(df, coords = c("X", "Y"), crs = terra::crs(Ts.eb)), 4326)
stopifnot(maxdiff(runH(ll)$H, H0$H) == 0)
sw <- df; sw$type <- ifelse(df$type == "hot", "cold", "hot")
expect_error(runH(sw), "swapped")
sw1 <- df; sw1$type[c(1, 6)] <- sw1$type[c(6, 1)]
expect_error(runH(sw1, method = "lm"), "pair")
far <- df; far$X[1] <- far$X[1] + 1e6
expect_error(runH(far), "anchor\\(s\\) 1 are outside the image")
na <- df; na[1, c("X", "Y")] <- terra::xyFromCell(Ts.eb, which(is.na(terra::values(Ts.eb)))[1])
expect_error(runH(na), "is NA at anchor\\(s\\) 1")
expect_error(runH(df[df$type == "hot", ]), "at least one hot and one cold")
expect_error(runH(transform(df, type = "warm")), "must be \"hot\" or \"cold\"")
# hot anchors 1-3 K warmer than the cold ones: the area lacks dry bare soil
cold.T <- mean(terra::extract(Ts.eb, xy[df$type == "cold", , drop = FALSE])[, 1])
warm <- which(terra::values(Ts.eb)[, 1] > cold.T + 1 & terra::values(Ts.eb)[, 1] < cold.T + 3)
flat.area <- df
flat.area[flat.area$type == "hot", c("X", "Y")] <- terra::xyFromCell(Ts.eb, warm[1:5])
expect_error(runH(flat.area), "minimum 5 K")

## Cold pixels and ET24h() use the same alfalfa reference ET (Allen et al.
## 2007): the reference ET fraction (ETrF) of the cold pixels is ETp.coef
cold <- terra::extract(EB$EB, terra::vect(EB$anchors[EB$anchors$type == "cold", ]))
lambda <- (2.501 - 0.00236 * (cold$surfaceTemperature - 273.15)) * 1e6
ETrF.cold <- mean(3600 * cold$LatentHeat / lambda) / hourlyET(WS, ET = "ETr")
stopifnot(abs(ETrF.cold - 1.05) < 0.03)

## ET24h() takes the METRIC.EB() result and returns ETinst, ETrF, ET24 and
## the H > Rn - G flag; the old call returns ET24 only, unchanged (#91)
et.old <- quiet(ET24h(EB$EB$NetRadiation, EB$EB$SoilHeat, EB$EB$SensibleHeat,
                      EB$EB$surfaceTemperature, WS, ETr.daily = suppressMessages(dailyET(WS))))
et.all <- quiet(suppressMessages(ET24h(EB)))
stopifnot(names(et.old) == "ET24",
          identical(names(et.all), c("ETinst", "ETrF", "ET24", "HgtRnG")),
          maxdiff(et.all$ET24, et.old) == 0,
          terra::global(et.all$HgtRnG, "sum", na.rm = TRUE)[1, 1] ==
            terra::global(EB$EB$SensibleHeat > EB$EB$NetRadiation - EB$EB$SoilHeat,
                          "sum", na.rm = TRUE)[1, 1],
          maxdiff(et.all$ETrF * hourlyET(WS, ET = "ETr"), et.all$ETinst) < 1e-12)

## "flexible": the grass-based Kc (Gonzalez et al. 2015) is converted to the
## alfalfa reference, ETp.coef = Kc * ETo / ETr
relaxed <- EB$anchors
relaxed$flex <- 0.05
relaxed$NDVI[relaxed$type == "cold"] <- 0.7
EBflex <- quiet(METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10,
                          WeatherStation = WS, MTL = MTL.L2, alb.coeff = "Olmedo",
                          anchors = relaxed, Z.om.ws = 0.03))
EBkc <- quiet(METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10,
                        WeatherStation = WS, MTL = MTL.L2, alb.coeff = "Olmedo",
                        anchors = relaxed, anchors.method = "best", Z.om.ws = 0.03,
                        ETp.coef = (2.11 * 0.7 - 0.4989) * hourlyET(WS, ET = "ETo") /
                          hourlyET(WS, ET = "ETr")))
stopifnot(maxdiff(EBflex$EB$LatentHeat, EBkc$EB$LatentHeat) < 1e-6)

## METRIC.EB(): weather station outside the image (#46)
WSfar <- WS
WSfar$location$lat <- WS$location$lat + 0.5  # about 55 km north
EBfar <- quiet(expect_warning(
  METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10, WeatherStation = WSfar,
            MTL = MTL.L2, alb.coeff = "Olmedo", n = 1, Z.om.ws = 0.03),
  "No pixel is within WSbuffer"))
stopifnot(inherits(EBfar, "waterLSEB"))

## METRIC.EB(): weather station outside the image but within WSbuffer: no
## fallback to the whole image
WSnear <- WS
WSnear$location$lat <- WS$location$lat + 0.06  # a few km north of the crop
w <- character()
invisible(utils::capture.output(withCallingHandlers(
  METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10, WeatherStation = WSnear,
            MTL = MTL.L2, alb.coeff = "Olmedo", n = 1, Z.om.ws = 0.03),
  warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })))
stopifnot(!any(grepl("WSbuffer", w)))

## METRIC.EB(): aoi in another CRS is reprojected to the image CRS
half <- terra::ext(image.SR); half[2] <- (half[1] + half[2]) / 2
aoi.ll <- sf::st_transform(sf::st_as_sf(terra::as.polygons(half, crs = terra::crs(image.SR))),
                           4326)
EBaoi <- quiet(METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10,
                         WeatherStation = WS, MTL = MTL.L2, alb.coeff = "Olmedo",
                         n = 1, Z.om.ws = 0.03, aoi = aoi.ll))
stopifnot(terra::ncell(EBaoi$EB) < terra::ncell(EB$EB))

## METRIC.EB(): manual anchors with the default anchors.method, and
## Zom.method = "custom"
EBman <- quiet(METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10,
                         WeatherStation = WS, MTL = MTL.L2, alb.coeff = "Olmedo",
                         anchors = EB$anchors, Z.om.ws = 0.03))
stopifnot(maxdiff(EBman$EB$LatentHeat, EB$EB$LatentHeat) < 1e-6)
EBturner <- quiet(METRIC.EB(image.SR = image.SR, surf.temp = image.SR$Band_10,
                            WeatherStation = WS, MTL = MTL.L2, alb.coeff = "Olmedo",
                            LAI.method = "turner", Z.om.ws = 0.03))
stopifnot(inherits(EBturner, "waterLSEB"))
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
          maxdiff(G, EBsw$EB$SoilHeat) < 1e-8,
          names(Rn) == "Rn", names(G) == "G")  # layers have names (#91)

## METRIC.Rn(): thermal layer named "thermal.low", detected as in METRIC.EB()
DN.low <- image.DN; names(DN.low)[names(DN.low) == "Thermal1"] <- "thermal.low"
Rn.low <- quiet(METRIC.Rn(image.DN = DN.low, image.SR = loadImageSR(c2), WeatherStation = WS,
                          MTL = MTL.L1, sat = "L8", alb.coeff = "Olmedo"))
stopifnot(inherits(Rn.low, "SpatRaster"))

## METRIC.EB() with Landsat 7 digital numbers, as in the METRIC_simple vignette
EB7 <- quiet(METRIC.EB(image.DN = L7, WeatherStation = WS7, MTL = MTL.L7, n = 5,
                       sat = "L7", thermalband = L7$Thermal1))
pin.L7 <- c(547.5577230608, 85.3656706699, 171.1774240015, 291.0163216654, 313.8947743373)
stopifnot(near(terra::global(EB7$EB, "mean", na.rm = TRUE)[, 1], pin.L7),
          setequal(terra::cellFromXY(EB7$EB, sf::st_coordinates(EB7$anchors)),
                   c(61345, 112247, 45977, 22960, 22009, 49292, 129603, 148901, 176846, 44752)))

cat("test-energy-balance.R: all tests passed\n")
