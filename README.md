# water

Actual Evapotranspiration with Energy Balance Models. Repository for the
development version of the water package.

water is not on CRAN at the moment (the CRAN version was archived on
2021-10-05). Install it from GitHub:

    install.packages("remotes")
    remotes::install_github("midraed/water", build_vignettes = TRUE)

Building the vignettes needs [pandoc](https://pandoc.org/) (RStudio includes
it). Without pandoc use `remotes::install_github("midraed/water")`; the
vignettes are then not installed, but you can read their sources on GitHub:
[Landsat 8 and 9](https://github.com/midraed/water/blob/master/vignettes/Landsat8.Rmd),
[METRIC simple](https://github.com/midraed/water/blob/master/vignettes/METRIC_simple.Rmd) and
[METRIC advanced](https://github.com/midraed/water/blob/master/vignettes/METRIC_advanced.Rmd).

On Windows, if the installation fails with "cannot open file" errors from a
dependency, the path of your R library is probably too long: install into a
short folder, e.g. `remotes::install_github("midraed/water", lib = "C:/Rlib")`
(create it first and add it with `.libPaths("C:/Rlib")`).

## Quick start (Landsat 8 and 9, Collection 2 Level-2)

    library(water)
    scene <- system.file("extdata", "L8_C2", package = "water")
    MTL <- list.files(scene, pattern = "MTL.txt$", full.names = TRUE)
    image.SR <- cfmask(path = scene, image = loadImageSR(scene))   # reflectance, Ts, clouds masked
    WS <- read.WSdata(system.file("extdata", "INTA.csv", package = "water"),
                      datetime.format = "%Y/%m/%d %H:%M",
                      columns = c(datetime = 1, radiation = 5, wind = 6, RH = 3, temp = 2, rain = 4),
                      lat = -33.00513, long = -68.86469, elev = 927, height = 2,
                      tz = "America/Argentina/Mendoza", MTL = MTL)
    EB <- METRIC.EB(image.SR = image.SR, WeatherStation = WS, MTL = MTL)
    summary(EB)      # calibration checks
    ET <- ET24h(EB)  # ETinst, ETrF, ET24 (mm/day)

You will see a warning if the radiation of the station at the overpass differs much from the clear-sky value (see the vignette). Start with the vignette "Landsat 8 and 9": `vignette("Landsat8", package = "water")`.
Two guides show how to search and download Level-2 scenes: with the
[USGS M2M API](https://github.com/midraed/water/blob/master/vignettes/Landsat_download_USGS.qmd)
and with [Microsoft Planetary Computer](https://github.com/midraed/water/blob/master/vignettes/Landsat_download_Planetarycomputer.qmd).

## Overview

water implements the METRIC surface energy balance model (Allen et al., 2007)
for Landsat 7, 8 and 9 images, including Landsat 8/9 Collection 2 Level-2
surface reflectance and surface temperature products, together with weather
station data processing and reference evapotranspiration (ASCE
Penman-Monteith). It uses the terra and sf packages.

The vignettes show the workflow (`vignette(package = "water")`):

- "Landsat 8 and 9": Collection 2 Level-2 example data. Start here.
- "METRIC simple" and "METRIC advanced": step by step, with Landsat 7 example
  data (Level-1 digital numbers).

Authors: Guillermo Federico Olmedo, Samuel Ortega-Farías, David Fonseca-Luengo,
Daniel de la Fuente-Sáiz and Fernando Fuentes-Peñailillo. Contributors: María
Victoria Munafó and Francisco Martín Corvalán.

Maintainer: Guillermo Federico Olmedo <guillermo.olmedo@gmail.com>

License: GPL (>= 2)

## Contributing

Contributions are accepted by [fork](https://docs.github.com/articles/fork-a-repo/) and [pull-request](https://docs.github.com/articles/using-pull-requests/)
