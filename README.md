# water

Actual Evapotranspiration with Energy Balance Models. Repository for the
development version of the water package.

[![](https://cranlogs.r-pkg.org/badges/grand-total/water)](https://cran.r-project.org/package=water)

water package is available at [CRAN](https://cran.r-project.org/package=water).

You can install the development version with:

    install.packages("remotes")
    remotes::install_github("midraed/water")

## Overview

water implements the METRIC surface energy balance model (Allen et al., 2007)
for Landsat 7, 8 and 9 images, including Landsat 8/9 Collection 2 Level-2
surface reflectance and surface temperature products, together with weather
station data processing and reference evapotranspiration (ASCE
Penman-Monteith). It uses the terra and sf packages.

The vignettes show the workflow (`vignette(package = "water")`):

- "METRIC simple" and "METRIC advanced": Landsat 7 example data.
- "Landsat 8 and 9": Landsat 8 Collection 2 Level-2 example data.

Authors: Guillermo Federico Olmedo, Samuel Ortega-Farías, David Fonseca-Luengo,
Daniel de la Fuente-Sáiz and Fernando Fuentes-Peñailillo. Contributors: María
Victoria Munafó and Francisco Martín Corvalán.

Maintainer: Guillermo Federico Olmedo <guillermo.olmedo@gmail.com>

License: GPL (>= 2)

## Contributing

Contributions are accepted by [fork](https://help.github.com/articles/fork-a-repo/) and [pull-request](https://help.github.com/articles/using-pull-requests/)
