## Submission water 1.0

water was archived on CRAN on 2021-10-05 (check problems not corrected in
time), so this is a new submission of a package that was on CRAN. Version
1.0 is a large update: it was migrated from raster/sp to terra/sf, supports
Landsat 8 and 9 Collection 2 Level-2 products, and fixes many errors in the
energy balance (see NEWS).

## Test environments
* local Windows 11, R 4.5.3 (R CMD check --as-cran, with --run-donttest)
* R-hub (GitHub Actions), R-devel: Linux and Windows, no errors
  (commit faa5078); macOS (Intel and Apple Silicon) passed on the previous
  commit (a587ab9, tests only changed afterwards).
* macOS arm64 R-devel on the final commit: TO FILL IN before submission.

## R CMD check results
There were no ERRORs or WARNINGs locally.

NOTEs:

* CRAN incoming feasibility: "New submission", "Package was archived on
  CRAN" (see above).
* On the local machine pandoc is not installed, so README.md and the
  vignettes could not be checked there.

## Not checked
win-builder and the vignette build (pandoc is not installed locally) are to
be run before submission. The Intel macOS R-devel job of R-hub currently fails
in its own setup step (it cannot download R-devel), before the package is
built.
