## Submission water 1.0

water was archived on CRAN on 2021-10-05 (check problems not corrected in
time), so this is a new submission of a package that was on CRAN. Version
1.0 is a large update: it was migrated from raster/sp to terra/sf, supports
Landsat 8 and 9 Collection 2 Level-2 products, and fixes many errors in the
energy balance (see NEWS).

## Test environments
* local Windows 11, R 4.5.3 (R CMD check --as-cran, with --run-donttest)

## R CMD check results
There were no ERRORs or WARNINGs.

NOTEs:

* CRAN incoming feasibility: "New submission", "Package was archived on
  CRAN" (see above).
* On the local machine pandoc is not installed, so README.md and the
  vignettes could not be checked there.

## Not checked locally
win-builder, r-hub and macOS checks, and the vignette build, are to be run
before submission.
