## Submission water 1.0

water was archived on CRAN on 2021-10-05 (check problems not corrected in
time), so this is a new submission of a package that was on CRAN. Version
1.0 is a large update: it was migrated from raster/sp to terra/sf, supports
Landsat 8 and 9 Collection 2 Level-2 products, and fixes many errors in the
energy balance (see NEWS).

## Test environments
* local Windows 11, R 4.5.3 (R CMD check --as-cran, with --run-donttest and
  the vignettes built with pandoc)
* win-builder, R 4.6.1 (release): 1 NOTE (see below)
* R-hub (GitHub Actions), R-devel: Linux and Windows, no errors (commits
  faa5078 and b13b4e8); macOS Intel and Apple Silicon passed on a587ab9
  (the later commits changed only tests, documentation and package
  metadata). The Intel macOS job of R-hub fails afterwards in its own setup
  step, before the package is built (it cannot download R-devel).

## R CMD check results
There were no ERRORs or WARNINGs.

NOTEs:

* CRAN incoming feasibility: "New submission" and "Package was archived on
  CRAN" (see above). The words flagged as possibly misspelled in DESCRIPTION
  (ASCE, Monteith, evapotranspiration, reflectance, "et al") are correct.
