## Resubmission 0.3.2

Previous submission had an error:
Found the following (possibly) invalid URLs:
  URL: https://cran.r-project.org/web/packages/water/
    From: README.md
    Status: 200
    Message: OK
    CRAN URL not in canonical form
  The canonical URL of the CRAN page for a package is
  https://cran.r-project.org/package=pkgname

* URL was corrected.


## Test environments
* local linux install, 3.2.5
* win-builder (devel and release)

## R CMD check
  results
There were no ERRORs or WARNINGs.

There
  was 1 NOTE:

* Possibly mis-spelled words in
  DESCRIPTION:
  Evapotranspiration (2:15)

  evapotranspiration (9:54)

Is not an error. The
  first use is upper case because it's on the title.




################################################################################


## Resubmission 0.3
This is a resubmission. In this version I have:

* Found (possibly) invalid URLs

SOLVED. Was a github.com repository marked as private, now marked as public.

* Title not in title case: SOLVED

* Author field differs from that derived from Authors@R: SOLVED

* Problems in NAMESPACE: ALL SOLVED

- ET24h: no visible global function definition for 'colorRampPalette': OK
- calcAnchors: no visible global function definition for 'points': OK
- calcH: no visible global function definition for 'points': OK
- plot.waterWeatherStation: no visible global function definition for
  'par': OK
- plot.waterWeatherStation: no visible global function definition for
  'abline': OK
- plot.waterWeatherStation: no visible global function definition for
  'points': OK
- plot.waterWeatherStation: no visible global function definition for
  'axis': OK
- plot.waterWeatherStation: no visible global function definition for
  'mtext': OK
- plot.waterWeatherStation: no visible global function definition for
  'axis.POSIXct': OK
- read.WSdata: no visible global function definition for 'read.csv': OK