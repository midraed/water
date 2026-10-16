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

## hourlyET(): only "ETr" and "ETo"
expect_error(hourlyET(ws, ET = "ETor"), "should be one of")

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

## aoiCrop(): honours waterOptions(autoAoi = FALSE)
aoi <- createAoi(topleft = c(510495, -3651000), bottomright = c(512000, -3652000), EPSG = 32619)
assign("aoi", aoi, envir = .GlobalEnv)
waterOptions(autoAoi = FALSE)
stopifnot(terra::ncell(water:::aoiCrop(image)) == terra::ncell(image))
waterOptions(autoAoi = TRUE)
stopifnot(terra::ncell(water:::aoiCrop(image)) < terra::ncell(image))
rm("aoi", envir = .GlobalEnv)

unlink(tmp, recursive = TRUE)
cat("test-inputs.R: all tests passed\n")
