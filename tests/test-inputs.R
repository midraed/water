# Regression tests for input handling: weather station data, MTL files,
# file loaders, cloud masks and DEMs. Plain R, run by R CMD check.
library(water)

expect_error <- function(expr, pattern) {
  msg <- tryCatch({ expr; NULL }, error = function(e) conditionMessage(e))
  if (is.null(msg)) stop("expected an error matching: ", pattern)
  if (!grepl(pattern, msg)) stop("unexpected error: ", msg)
  invisible(msg)
}
expect_warning <- function(expr, pattern) {
  w <- NULL
  value <- withCallingHandlers(expr, warning = function(c) {
    if (grepl(pattern, conditionMessage(c))) w <<- conditionMessage(c)
    invokeRestart("muffleWarning")
  })
  if (is.null(w)) stop("expected a warning matching: ", pattern)
  invisible(value)
}

extdata <- system.file("extdata", package = "water")
MTL.L8 <- file.path(extdata, "LC82320832016040LGN00_MTL.txt")
tmp <- tempfile("water-tests-"); dir.create(tmp)

# Synthetic weather station: hourly records around a July overpass, in UTC
ws.july <- function() {
  t <- seq(as.POSIXct("2016-07-09 00:00", tz = "UTC"), by = "1 hour", length.out = 48)
  data.frame(datetime = format(t, "%Y-%m-%d %H:%M:%S"), radiation = 500, wind = 2,
             RH = 50, temp = 20, rain = 0)
}
mtl.july <- file.path(tmp, "july_MTL.txt")
writeLines(sub("DATE_ACQUIRED = 2016-02-09", "DATE_ACQUIRED = 2016-07-09",
               readLines(MTL.L8)), mtl.july)
read.ws <- function(data, ...) {
  read.WSdata(WSdata = data, columns = c(1, 1, 2, 3, 4, 5, 6), lat = -33,
              long = -68.9, elev = 927, height = 2, ...)
}

## read.WSdata(): the MTL overpass time is UTC (it was read as Europe/London)
ws <- read.ws(ws.july(), tz = "UTC", MTL = mtl.july)
overpass <- as.POSIXct("2016-07-09 14:27:29", tz = "UTC")
stopifnot(abs(as.numeric(difftime(ws$at.sat$datetime, overpass, units = "mins"))) < 1)

## read.WSdata(): clear errors and robust inputs (#41)
expect_error(read.ws(ws.july(), tz = "GTM-6"), "Unknown time zone")
expect_error(read.ws(ws.july(), tz = "UTC", datetime.format = "%d/%m/%Y %H:%M"),
             "No date could be read")
bad <- ws.july(); bad$datetime[3] <- "not a date"
ws.bad <- expect_warning(read.ws(bad, tz = "UTC"), "could not be read")
stopifnot(nrow(ws.bad$alldata) == 47)
tbl <- ws.july(); class(tbl) <- c("tbl_df", "tbl", "data.frame")
stopifnot(inherits(read.ws(tbl, tz = "UTC"), "waterWeatherStation"))
ws.int <- read.WSdata(WSdata = ws.july(), columns = c(1L, 1L, 2L, 3L, 4L, 5L, 6L),
                      lat = -33, long = -68.9, elev = 927, tz = "UTC")
stopifnot(inherits(ws.int, "waterWeatherStation"))
chr <- ws.july(); chr$temp <- sub("\\.", ",", format(chr$temp + 0.5))
expect_error(read.ws(chr, tz = "UTC"), "not numeric")

## read.WSdata(): a station without rain data (all-NA column, read as logical)
norain <- ws.july(); norain$rain <- NA
stopifnot(all(is.na(read.ws(norain, tz = "UTC")$alldata$rain)))

## read.WSdata(): any time zone that R knows; an empty MTL means no MTL
stopifnot(inherits(read.ws(ws.july(), tz = "Etc/GMT+3"), "waterWeatherStation"))
stopifnot(is.null(read.ws(ws.july(), tz = "UTC", MTL = character(0))$at.sat))

## read.WSdata(): overpass outside the records warns, no at.sat (#28)
ws.feb <- expect_warning(read.ws(ws.july(), tz = "UTC", MTL = MTL.L8),
                         "outside the weather station records")
stopifnot(is.null(ws.feb$at.sat))
expect_error(water:::getDataWS(ws.feb), "no conditions at the satellite overpass")

## MTL must be a file name, not its contents (#45)
expect_error(read.ws(ws.july(), tz = "UTC", MTL = utils::read.delim(mtl.july)),
             "MTL must be the path")

## read.WSdata2() reads semicolon / decimal comma files
csv2 <- file.path(tmp, "ws.csv"); d2 <- ws.july(); d2$temp <- 20.5
utils::write.csv2(d2, csv2, row.names = FALSE)
ws2 <- read.WSdata2(WSdata = csv2, columns = c(1, 1, 2, 3, 4, 5, 6), lat = -33,
                    long = -68.9, elev = 927, height = 2, tz = "UTC",
                    datetime.format = "%Y-%m-%d %H:%M:%S")
stopifnot(all(ws2$alldata$temp == 20.5))
# ... and read.WSdata() points to it (it failed with "duplicate 'row.names'")
expect_error(read.ws(csv2, tz = "UTC"), "read.WSdata2")

## read.WSdata(): named columns are matched by name, in any order (#83)
july <- ws.july()
july$temp <- 20 + seq_len(48) %% 5
july$RH <- 40 + seq_len(48) %% 7
ref <- read.ws(july, tz = "UTC", MTL = mtl.july)
by.name <- function(data, columns, ...) {
  read.WSdata(WSdata = data, columns = columns, lat = -33, long = -68.9, elev = 927,
              height = 2, tz = "UTC", ...)
}
swapped <- july[, c("datetime", "temp", "RH", "wind", "radiation", "rain")]
ws.named <- expect_warning(by.name(swapped, c(datetime = 1, temp = 2, RH = 3, wind = 4,
                                              radiation = 5, rain = 6), MTL = mtl.july),
                           "matched by name")
stopifnot(isTRUE(all.equal(ws.named$alldata, ref$alldata)),
          isTRUE(all.equal(ws.named$at.sat, ref$at.sat)))
expect_error(by.name(july, c(datetime = 1, humidity = 4)), "Unknown names in columns")
expect_error(by.name(july, c(datetime = 1, radiation = 2, wind = 3, RH = 4)),
             "columns has no \"temp\"")

## read.WSdata(): units and plausible ranges (#83)
f <- july; f$temp <- f$temp * 1.8 + 32
expect_error(read.ws(f, tz = "UTC"), "Fahrenheit")
k <- july; k$temp <- k$temp + 273.15
expect_error(read.ws(k, tz = "UTC"), "Kelvin")
fr <- july; fr$RH <- fr$RH / 100
expect_error(read.ws(fr, tz = "UTC"), "fraction")
kmh <- july; kmh$wind <- 30
invisible(expect_warning(read.ws(kmh, tz = "UTC"), "km/h"))
mj <- july; mj$radiation <- mj$radiation / 277.8
invisible(expect_warning(read.ws(mj, tz = "UTC"), "MJ"))

## read.WSdata(): records are sorted; duplicates removed or averaged (#83).
## Reversed records failed, shuffled ones gave a wrong at.sat
set.seed(1)
stopifnot(isTRUE(all.equal(read.ws(july[48:1, ], tz = "UTC", MTL = mtl.july)$at.sat, ref$at.sat)),
          isTRUE(all.equal(read.ws(july[sample(48), ], tz = "UTC", MTL = mtl.july)$at.sat,
                           ref$at.sat)))
ws.dup <- expect_warning(read.ws(rbind(july, july[10, ]), tz = "UTC", MTL = mtl.july),
                         "duplicated records removed")
stopifnot(isTRUE(all.equal(ws.dup$at.sat, ref$at.sat)))
conflict <- july[c(1:48, 10), ]; conflict$temp[49] <- 30
ws.avg <- expect_warning(read.ws(conflict, tz = "UTC"), "averaged")
stopifnot(nrow(ws.avg$alldata) == 48, ws.avg$alldata$temp[10] == (july$temp[10] + 30) / 2)
twice <- rbind(july, transform(july, temp = temp + 1))
expect_error(read.ws(twice, tz = "UTC"), "share their time")

## hourlyET(): solar time correction for the centre of the time zone of the
## station clock, from its UTC offset (#77); it used the station longitude
apples <- read.WSdata(WSdata = file.path(extdata, "apples.csv"), date.format = "%d/%m/%Y",
                      lat = -35.42222, long = -71.38639, elev = 201, height = 2.2,
                      columns = c(1, 2, 3, 4, 6, 7, 8), cf = c(1, 0.2777778, 1),
                      tz = "America/Santiago", MTL = file.path(extdata, "L7.MTL.txt"))
overpass <- water:::getDataWS(apples)
stopifnot(overpass$long.z == -45)  # UTC-3 in February (daylight saving time)
et.default <- hourlyET(apples)
stopifnot(all.equal(et.default, hourlyET(overpass, hours = overpass$hours, DOY = overpass$DOY,
                                         ET.instantaneous = TRUE, long.z = -45)),
          abs(et.default - hourlyET(overpass, hours = overpass$hours, DOY = overpass$DOY,
                                    ET.instantaneous = TRUE, long.z = overpass$long)) > 0.005)
stopifnot(all.equal(dailyET(apples, MTL = file.path(extdata, "L7.MTL.txt")),
                    dailyET(apples, MTL = file.path(extdata, "L7.MTL.txt"), long.z = -45)))

## dailyET(): without MTL, the day is date or the day of the overpass
## conditions (it summed the first 24 records, another day) (#87)
et.mtl <- dailyET(apples, MTL = file.path(extdata, "L7.MTL.txt"))
stopifnot(all.equal(dailyET(apples), et.mtl),
          all.equal(dailyET(apples, date = "2013-02-15"), et.mtl))
expect_error(dailyET(read.ws(ws.july(), tz = "UTC")), "needs the day")
expect_error(dailyET(apples, date = "2001-01-01"), "No weather station records")

## print() returns the weather station invisibly
printed <- withVisible(print(apples))
stopifnot(!printed$visible, identical(printed$value, apples))

## mtlValues(): keys repeated in several groups need a group (#78)
mtl.l2 <- list.files(file.path(extdata, "L8_C2"), pattern = "MTL.txt$", full.names = TRUE)
expect_error(water:::mtlValues(mtl.l2, "REFLECTANCE_MULT_BAND_4"), "more than once")
stopifnot(water:::mtlValues(mtl.l2, "REFLECTANCE_MULT_BAND_4",
                            group = "LEVEL2_SURFACE_REFLECTANCE_PARAMETERS") == 2.75e-05,
          water:::mtlValues(mtl.l2, "REFLECTANCE_MULT_BAND_4",
                            group = "LEVEL1_RADIOMETRIC_RESCALING") == 2e-05)

## Method arguments: exact names (no partial matching), with a hint for
## typos, checked before any raster work (#82)
expect_error(hourlyET(ws, ET = "ETor"), "Unknown ET")
expect_error(hourlyET(ws, ET = "ETR"), "did you mean \"ETr\"")
expect_error(dailyET(ws, ET = "E"), "Unknown ET")
one <- terra::rast(nrows = 2, ncols = 2, vals = 1)
expect_error(momentumRoughnessLength(method = "short.crop", LAI = one),
             "did you mean \"short.crops\"")
expect_error(momentumRoughnessLength(method = "Perrier", LAI = one, fLAI = 2, h = 1),
             "fLAI must be")
expect_error(momentumRoughnessLength(method = "custom", LAI = one, NDVI = one, albedo = one),
             "needs a and b")
expect_error(soilHeatFlux(method = "T"), "Unknown method = \"T\"")
expect_error(LAI(method = "METRIC2010", image = one), "did you mean \"metric2010\"")
expect_error(albedo(one, coeff = "olmedo"), "did you mean \"Olmedo\"")
expect_error(calcH(method = "LM"), "did you mean \"lm\"")
expect_error(surfaceTemperature(sat = "L7", method = "SW"), "split window")
# image.SR is never used: the errors come first
for (a in list(list(alb.coeff = "olmedo"), list(LAI.method = "METRIC2010"),
               list(Zom.method = "perrier"), list(anchors.method = "Best"),
               list(LST.method = "sw"), list(G.method = "T"), list(sat = "L 8"))) {
  expect_error(do.call(METRIC.EB, c(list(image.SR = "not used", MTL = MTL.L8), a)),
               "Unknown")
}
expect_error(METRIC.EB(image.SR = "x", MTL = MTL.L8, Zom.method = "custom"),
             "needs extraParameters a and b")
expect_error(METRIC.EB(image.SR = "x", MTL = MTL.L8, extraParameters = c(minDist = "500")),
             "named numeric")
expect_error(METRIC.EB(image.SR = "x", MTL = MTL.L8, extraParameters = c(WSBuffer = 1000)),
             "Unknown extraParameters: WSBuffer")
expect_error(METRIC.Rn(image.SR = "x", MTL = MTL.L8, sat = "L7", LST.method = "SW"),
             "split window")
expect_error(METRIC.G(image.SR = "x", MTL = MTL.L8, G.method = "bastiaanssen"),
             "did you mean \"Bastiaanssen\"")

## getSat(): satellite from the MTL file
stopifnot(water:::getSat(MTL = MTL.L8) == "L8")
mtl9 <- file.path(tmp, "L9_MTL.txt")
writeLines(sub('"LANDSAT_8"', '"LANDSAT_9"', readLines(MTL.L8)), mtl9)
stopifnot(water:::getSat(MTL = mtl9) == "L9")

## loadImage(): real-time tier and renamed files, clear error (#50)
make.scene <- function(dir, prefix, bands = c(2:7, 10, 11)) {
  dir.create(dir)
  for (b in bands) {
    terra::writeRaster(terra::rast(nrows = 3, ncols = 3, vals = b),
                       file.path(dir, paste0(prefix, b, ".TIF")))
  }
  dir
}
for (prefix in c("LC08_L1TP_232083_20160209_20200907_02_RT_B",
                 "LC08_L1TP_232083_20160209_20200907_02_T1_2020-03-29_B",
                 "LC09_L1TP_232083_20230209_20230310_02_T1_B")) {
  img <- loadImage(path = make.scene(tempfile(tmpdir = tmp), prefix), sat = "L8")
  stopifnot(terra::nlyr(img) == 8, terra::global(img$Thermal2, "max")[1, 1] == 11)
}
expect_error(loadImage(path = tmp, sat = "L8"), "No Landsat Level-1 band files")
# DN 0 is fill (#85)
z.dir <- make.scene(tempfile(tmpdir = tmp), "LC08_L1TP_232083_20160209_20200907_02_T1_B")
terra::writeRaster(terra::rast(nrows = 3, ncols = 3, vals = c(0, rep(2, 8))), overwrite = TRUE,
                   file.path(z.dir, "LC08_L1TP_232083_20160209_20200907_02_T1_B2.TIF"))
z.img <- terra::values(loadImage(path = z.dir, sat = "L8"))
stopifnot(is.na(z.img[1, "B"]), !anyNA(z.img[-1, ]))

## loadImageSR(): all six bands are required
sr.dir <- tempfile(tmpdir = tmp); dir.create(sr.dir)
for (b in c(2:4, 6:7)) {
  terra::writeRaster(terra::rast(nrows = 3, ncols = 3, vals = 1),
                     file.path(sr.dir, sprintf("LC08_X_SR_B%d.TIF", b)))
}
expect_error(loadImageSR(sr.dir), "Missing surface reflectance bands")
# only Collection 2 Level-2 files (the ESPA format is no longer read)
espa.dir <- tempfile(tmpdir = tmp); dir.create(espa.dir)
for (b in 2:7) {
  terra::writeRaster(terra::rast(nrows = 3, ncols = 3, vals = 1),
                     file.path(espa.dir, sprintf("LC82320832016040LGN00_sr_band%d.tif", b)))
}
expect_error(loadImageSR(espa.dir), "No Collection 2 Level-2")
# fill and values outside the valid range (reflectance 0 to 1) are masked in
# every layer, surface temperature included (#85)
sr.ok <- tempfile(tmpdir = tmp); dir.create(sr.ok)
for (b in 2:7) {
  v <- rep(10000, 9); if (b == 2) v[1:3] <- c(5000, 0, 50000)
  terra::writeRaster(terra::rast(nrows = 3, ncols = 3, vals = v), datatype = "INT2U",
                     file.path(sr.ok, sprintf("LC08_X_SR_B%d.TIF", b)))
}
terra::writeRaster(terra::rast(nrows = 3, ncols = 3, vals = 40000), datatype = "INT2U",
                   file.path(sr.ok, "LC08_X_ST_B10.TIF"))
msg <- character()
sr3 <- withCallingHandlers(loadImageSR(sr.ok), message = function(m) {
  msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage")
})
v3 <- terra::values(sr3)
stopifnot(all(is.na(v3[1:3, ])), !anyNA(v3[4:9, ]),
          all(v3[4:9, "B"] == 10000 * 0.0000275 - 0.2), any(grepl("^3 pixels", msg)))

## cfmask(): Collection 2 QA_PIXEL bits (#55)
img <- terra::rast(nrows = 1, ncols = 5, vals = 1:5)
qa <- terra::rast(img, vals = c(21824, 21952, 22280, 23888, 1)); names(qa) <- "QA_PIXEL"
masked <- cfmask(image = img, cfmask = qa, buffer = 0)
stopifnot(isTRUE(all.equal(as.numeric(terra::values(masked)), c(1, 2, NA, NA, NA))))
img5 <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 150, ymin = 0, ymax = 150,
                    vals = 1, crs = "EPSG:32719")
qa5 <- terra::rast(img5, vals = 21824); qa5[13] <- 22280; names(qa5) <- "QA_PIXEL"
stopifnot(terra::global(is.na(cfmask(image = img5, cfmask = qa5, buffer = 0)), "sum")[1, 1] == 1,
          terra::global(is.na(cfmask(image = img5, cfmask = qa5, buffer = 30)), "sum")[1, 1] == 5,
          terra::global(is.na(cfmask(image = img5, cfmask = qa5, buffer = 45)), "sum")[1, 1] == 9)
names(qa) <- "qa_pixel"  # e.g. renamed by other tools
stopifnot(isTRUE(all.equal(as.numeric(terra::values(cfmask(image = img, cfmask = qa, buffer = 0))),
                           c(1, 2, NA, NA, NA))))
qa.old <- terra::rast(img, vals = c(0, 0, 2, 0, 4))
stopifnot(isTRUE(all.equal(as.numeric(terra::values(cfmask(image = img, cfmask = qa.old, buffer = 0))),
                           c(1, 2, NA, 4, NA))))

## prepareSRTMdata(): DEM on the image grid (#31)
image <- terra::rast(file.path(extdata, "LC82320832016040LGN00_band2.tif"))
srtm <- terra::rast(xmin = -69, xmax = -68, ymin = -34, ymax = -33, res = 1/360,
                    crs = "EPSG:4326", vals = 900)
srtm.dir <- tempfile(tmpdir = tmp); dir.create(srtm.dir)
terra::writeRaster(srtm, file.path(srtm.dir, "s34_w069_1arc_v3.tif"))
DEM <- prepareSRTMdata(path = srtm.dir, extent = image)
stopifnot(terra::compareGeom(DEM, image, stopOnError = FALSE))
expect_error(prepareSRTMdata(path = srtm.dir, extent = createAoi(c(0, 1), c(1, 0))),
             "no CRS")
# a polygon gives a DEM on the Landsat 30 m grid (cell edges at 15 m + k * 30 m)
aoi.utm <- createAoi(topleft = c(511000, -3651000), bottomright = c(513000, -3653000),
                     EPSG = 32619)
DEM.aoi <- prepareSRTMdata(path = srtm.dir, extent = aoi.utm)
stopifnot(all(terra::res(DEM.aoi) == 30), terra::xmin(DEM.aoi) %% 30 == 15)
expect_error(prepareSRTMdata(path = srtm.dir, extent = sf::st_transform(aoi.utm, 4326)),
             "longitude/latitude")

## calcRadiance(): Landsat 7 MTL values have no "E-" notation (they gave NaN)
rad7 <- calcRadiance(terra::rast(L7_Talca), sat = "L7",
                     MTL = file.path(extdata, "L7.MTL.txt"))
stopifnot(!anyNA(terra::global(rad7, "mean", na.rm = TRUE)[, 1]))

## checkSRTMgrids(): file names of the tiles needed (it built outdated links)
stopifnot("s34_w069_1arc_v3.tif" %in% suppressMessages(checkSRTMgrids(image)))

## waterOptions(): output folder (trim() was a raster function)
waterOptions(outputFolder = tmp)
stopifnot(getOption("waterOutputFolder") == paste0(tmp, "/"))
waterOptions(default = TRUE)

## No hidden inputs: an object called aoi in the workspace and MTL files in
## the working directory are never used (#87)
assign("aoi", createAoi(c(0, 1), c(1, 0), EPSG = 32719), envir = .GlobalEnv)
stopifnot(terra::ncell(loadImageSR(file.path(extdata, "L8_C2"))) == 134 * 184)
rm("aoi", envir = .GlobalEnv)
expect_warning(waterOptions(autoAoi = TRUE), "no longer used")
old <- setwd(tmp)
invisible(file.copy(c(MTL.L8, mtl.july), c("a_MTL.txt", "b_MTL.txt")))
msg <- character()
ws0 <- withCallingHandlers(read.ws(ws.july(), tz = "UTC"), message = function(m) {
  msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage")
})
stopifnot(is.null(ws0$at.sat), any(grepl("at.sat", msg)))
expect_error(water:::getSat(), "needs the MTL")
expect_error(solarAngles(surface.model = NULL, WeatherStation = ws), "MTL is required")
expect_error(METRIC.Rn(image.SR = "x", WeatherStation = ws), "needs the MTL")
setwd(old)

unlink(tmp, recursive = TRUE)
cat("test-inputs.R: all tests passed\n")
