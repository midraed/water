# water 0.8.5 (master 70027b3): code review by Gabe (R package developer)

Scope: whether the code is correct, how it performs on full scenes, how consistent the API is, and package etiquette. Every finding below was reproduced with R 4.5.3, terra 1.9.11 and sf 1.1.0 on Windows 11 (16 cores, 64 GB RAM; other reviewers' jobs were running at the same time). Before testing I checked that the installed `lib_new` matches the repository: all 53 functions are identical. The reproduction scripts are in `scratchpad/review2/Gabe/` (`rNN_*.R`; each sources `pre.R`). File:line references are to the repository.

## R CMD check and tests

- `R CMD build --no-build-vignettes` followed by `R CMD check --as-cran --no-manual --ignore-vignettes`: **0 ERRORs, 0 WARNINGs, 3 NOTEs**.
  - CRAN incoming: "New submission / Package was archived on CRAN (2021-10-05)". The same NOTE says "VignetteBuilder field but no prebuilt vignette index", which only appears because the build skipped vignettes.
  - "unable to verify current time" (environment).
  - "README.md/NEWS.md cannot be checked without pandoc" (environment).
- The install is clean (no "replacing previous import"). "R code for possible problems": OK.
- Timings:
  - `tests/test-energy-balance.R`: **193 s**. `test-inputs.R`: about 10 s.
  - The `METRIC.EB` example takes 40 s under `--run-donttest`. All other examples take under 2.5 s.
  - `getDEM` is `\dontrun` because it needs the internet. That is acceptable.
- **Test coverage (covr, tests only): 68.2 %**. See finding 18 for the gaps.
- Because the package was archived, CRAN will review it as a new submission (see finding 19).

---

## Findings

### 1. HIGH. `waterOptions(writeResults = TRUE)` silently corrupts results
**Evidence**
- `saveLoadClean()` (R/water_internals.R:93-114) always writes to `<outputFolder>/<fixed name>.tif` and overwrites the file (callers pass `overwrite = TRUE`; `waterOverwrite = TRUE` only drops the time stamp). It then returns `rast(file)`.
- Any later call that writes the same name replaces the file under every SpatRaster object that already points to it.

**Reproduction** (`r05b_write.R`, `r05c_write.R`)
- `lai1 <- LAI(image = sr)` has mean 1.2545. After `LAI(method = "turner", image = sr)`, the same `lai1` object has mean **0.5679**.
- `METRIC.EB(..., Zom.method = "Perrier", extraParameters = c(fLAI = 0.5, h = 1.8))` gives mean LE 194.8 W/m2 without writing and **226.9 W/m2** with writing (max difference 81 W/m2).
  - Cause: R/water_METRIC.R:431-440 writes `Z.om.tif` twice (the Perrier `Z.om`, then the short-crops `Z.om.sc`).
  - With writing on, the "Perrier" run is identical to a short.crops run (max difference 1e-4).
- Running two scenes in one session has the same problem: the first result's `EB.tif` is overwritten by the second.

**Why it matters:** this option exists for large scenes. Using it gives wrong numbers with no error or warning.

**Proposed change**
- Write to unique file names, for example `tempfile(pattern = file, tmpdir = outputFolder, fileext = ".tif")`, or let terra manage temporary files (`terraOptions(todisk = TRUE)`).
- Write user-facing files only for final outputs, and never overwrite a file that a live object references.
- Add a test (the writeResults branch currently has 17 % coverage).

### 2. HIGH. A misspelled method name silently returns stale objects from the user's workspace, or fails with a cryptic error
**Evidence**
- No exported function validates its method, coefficient or satellite strings (there is no `match.arg`). The only exception is `hourlyET(ET=)`.
- Every branch is an `if (method == "...")` without an `else`. When nothing matches, the result variable is never created, so R looks it up in the enclosing environments, which end at `.GlobalEnv`.

**Reproduction** (`r15c_global.R`, `r15_args.R`, `r25_extra.R`). The objects below have the names the METRIC_advanced vignette creates and were computed for a 24,656-cell image; the new image has 12,328 cells.
- `momentumRoughnessLength(method = "short.crop", ...)` returns the **old global `Z.om`** (identical, 24,656 cells).
- `soilHeatFlux(..., method = "tasumi")` returns the **old global `G`**.
- `surfaceTemperature(sat = "L 8", ...)` returns the old global `Ts`.
- Other typos inside `METRIC.EB()`:

  | Argument | Error message | After |
  |---|---|---|
  | `alb.coeff = "olmedo"` | `object 'wb' not found` | |
  | `LAI.method = "METRIC2010"` | `comparison (<) is possible only for atomic and list types` | |
  | `Zom.method = "perrier"` | `object 'Z.om' not found` | 4 s on the crop (minutes on a scene) |
  | `anchors.method = "Best"` | `object 'hot' not found` | |
  | `LST.method = "sw"` | `object 'Ts' not found` | |
  | `G.method = "tasumi"` | an unrelated error inside calcH, or a stale global `G` used | |

- Missing required extra parameters also fail cryptically:
  - `Zom.method = "custom"` without `a`/`b`: `need finite 'ylim' values`.
  - `"Perrier"` without `fLAI`/`h`: `missing value where TRUE/FALSE needed`.
  - `extraParameters = c(minDist = "500")`: `[<=] only '==' is supported with categorical comparisons`.

**Proposed change**
- Call `match.arg()` for `method`, `coeff`, `sat`, `anchors.method`, `Zom.method`, `G.method`, `LST.method`, `LAI.method` and `alb.coeff` at the top of each exported function. Validate them all in `METRIC.EB()` before any raster work.
- Use `switch(method, ..., stop("unknown method"))`.
- Validate `extraParameters` (numeric; required names present).

### 3. HIGH. `calcH()`: the r_ah iteration oscillates at low wind, and the result depends on where the loop stops
**Evidence**
- R/water_sensibleHeatFlux.R:651-738 is an undamped fixed-point iteration.
- Convergence requires a change of less than 1 % in r_ah at both the hot and the cold anchors. On failure the loop just `break`s with a warning and returns the last iterate.

**Reproduction** (`r13b/c/d_*.R`, example crop, station wind set to 0.53-0.70 m/s, u200 about 1.0-1.4 m/s)
- The cold-pixel r_ah alternates in a 2-cycle: 134.6, 0.5, 127.3, 0.8, 118.7 ... s/m.
- dT_cold alternates in step: 27.4, 0.1, 26.0 K.
- After 19 iterations the function warns "No convergence after 20 iterations".
- Mean LE by `maxit`: **199.1 (18), 208.2 (19), 201.2 (20), 207.3 (21)** W/m2. The output is whichever half of the oscillation the loop stopped on.

**Why it matters:** wind of 0.5-1 m/s at a 10:30 overpass is common in irrigated valleys.

**Proposed change**
- Damp the update, e.g. `r.ah <- (r.ah.new + r.ah.old)/2`, or iterate in log space.
- Test convergence on dT_hot/dT_cold (a and b), not only on r_ah.
- Return a `converged` flag and the iteration history in the result.

**Also in the loop logic (`r12_calcH.R`):**
- `i == maxit` is off by one: `maxit = 3` runs 2 iterations, and the warning says "after 3".
- `maxit = 1` is ignored: the loop ran 13 iterations. With a non-converging case it would loop forever.
- `method = "lm"` uses `>=`, which is inconsistent with `"mean"`.

### 4. HIGH. `dailyET()` integrates the wrong day east of UTC, and the `date` argument uses UTC days
**Evidence:** R/water_ET.R:236-249.
- `date = "auto"` compares the local station date with the MTL `DATE_ACQUIRED`, which is a UTC date.
- An explicit `date` uses `as.Date(POSIXct)`, which works in UTC.
- The day is integrated with `for(i in 1:24)` and `sum(na.rm = TRUE)`.

**Reproduction** (`r06_dailyET_tz.R`, `r06b_dailyET.R`). Synthetic New Zealand station (Pacific/Auckland); overpass 2016-02-09 10:30 NZDT, which is 2016-02-08 21:30 UTC:
- `read.WSdata()` finds the overpass correctly.
- `dailyET(auto)` returns **4.85 mm**, the 8 February total. The true 9 February (overpass day) total is **14.66 mm**.
- `date = "2016-02-09"` selects 11 rows (13:00-23:00 local).
- INTA example (UTC-3): `date = "2016-02-09"` uses 21 of the 24 hours (5.311 vs 5.340 mm).
- A date with no records returns **0 with no warning**.
- `date = as.Date(...)` errors with "character string is not in a standard unambiguous format".

**Why it matters:** this affects every user from about UTC+11 eastwards (New Zealand; eastern Australia in summer) and every explicit `date`.

**Proposed change**
- Select hours by the local date of `at.sat$datetime`: `format(datetime, "%Y-%m-%d", tz = attr(datetime, "tzone"))`.
- Accept Date or character for `date`.
- Stop when no rows are found, and warn when fewer than 24 hourly values exist.

### 5. HIGH. Full-scene performance: one scene takes about 1 hour and about 33 GB of RAM
**Reproduction** (`r08_full.R`, `r08_full.log`). LC08 232/083 Collection 2 Level-2, 7801 x 7741 = 60.4 M cells. One double layer is 461 MiB; a logical vector is 230 MiB. Steps run one by one, as `METRIC.EB()` does:

| Step | Time | Process peak working set |
|---|---|---|
| `loadImageSR` | 51 s | 14.2 GB |
| `cfmask` | 31 s | 14.6 GB |
| `metricSurface` | 297 s | 30.4 GB |
| `metricRn` | 297 s | 33.3 GB |
| `soilHeatFlux` | 50 s | |
| `calcAnchors` (n = 1) | 251 s | |
| `calcAnchors` (n = 5) | 461 s | |
| **`calcH`** | **2315 s** | |

- Total for n = 5: **about 58 min**.
- `calcH` also wrote 5.9 GB of terra temporary files.

**Causes (each measured or counted)**
- **calcH**: every iteration recomputes about 20 full rasters (dT, rho, H, L, three masked `phi` replacements, u*, r_ah). On top of that there are about 20 `r[cells]` extractions per iteration, which explains the 7.6 s on the 24k-cell crop: profiling shows the time is terra S4 object overhead.
- **calcAnchors** makes **97** whole-raster `values()` passes for n = 1 and **205 passes plus 8 whole-raster `terra::buffer()` calls** for n = 5 (`r26_count.R`). One pass costs about 2.6 s on the full scene (`r27_plot.R`), so 97 × 2.6 ≈ 250 s, matching the measurement. Each relaxation step of the "flexible" method recomputes about 8 comparisons and `values()` calls.
- On an all-masked or unsuitable scene, "flexible" runs all 200 relaxation steps before failing (`r29_allna.R`: 18 s and 200 printed lines on the crop). At 8 passes × 2.6 s per step that is roughly 70 min on a full scene.
- Constant fields are stored as full in-memory rasters via `values<-`:
  - DEM for `plain = TRUE` (R/water_METRIC.R:171-172).
  - Declination and hour angle (R/water_netRadiation.G.R:264-274).
- `albedo()` uses `app()` with an R closure for each band (:370-373) instead of `sum(image.SR[[1:6]] * wb)`.
- `SWtrasmisivity()` is computed twice (once for incoming short-wave, once for long-wave) and uses `[`/`[<-` extractions.
- The `METRIC.EB` example takes 40 s and the tests 193 s on a 134 x 184 crop.

**Proposed change**
- **calcH:** iterate the calibration on the anchor values only (scalars), store the sequence of (a, b), and then run the per-pixel stability iterations in a single `terra::lapp()` pass with a vectorised function. This gives the same result with one pass instead of about 20 × k.
- **calcAnchors:** read LAI, albedo, Z.om, NDVI, Ts and the buffer **once** into vectors. Compute for each pixel the relaxation it needs (the maximum over the criteria), then pick the smallest flex that leaves more than n candidates. Enforce `minDist` with the candidates' coordinates instead of `terra::buffer()`. Stop early when there are no non-NA pixels or no pixel with NDVI > 0.5.
- Use scalars for constants, plain arithmetic in `albedo()`, and compute transmissivity once.

### 6. MEDIUM. Hidden dependence on a global object called `aoi`
**Evidence:** `aoiCrop()` (R/water_internals.R:117-129) is used by `loadImage`, `loadImageSR` and `albedo`. When `waterAutoAoi = TRUE` (the default), it uses any `aoi` object it finds in `.GlobalEnv`. The METRIC vignettes create such an object.

**Reproduction** (`r03_aoi.R`)
- With the Talca vignette `aoi` still in the workspace, `loadImageSR(c2)` and `albedo(sr)` fail with `[crop] extents do not overlap`.
- With an overlapping `aoi`, `albedo(sr)` has 3,600 cells while `LAI(sr)` has 24,656, and `netRadiation()` fails with `[*] extents do not match`.

**Proposed change:** remove `autoAoi`, or make its default FALSE, and never read `.GlobalEnv`. `metricSurface()` already has to pass `aoi = ext(image.SR)` to work around this.

### 7. MEDIUM. Hidden dependence on `getwd()` for the MTL file
**Evidence:** when `MTL` is missing, these functions scan the working directory:
- `read.WSdata` (R/water_WeatherStation.R:199)
- `dailyET` (:230)
- `METRIC.*`, `solarAngles`, `calcTOAr`, `calcRadiance`, `surfaceTemperature`, `getSat`

**Reproduction** (`r10_getwd.R`): `read.WSdata(apples.csv, ...)` without `MTL`:
- With one unrelated MTL file in the working directory, it warns "overpass (2016-02-09) is outside the weather station records" and drops `at.sat`.
- With two MTL files, it errors: "More than one MTL file in <wd> ... Pass the one to use with MTL".
- `dailyET(WS)` gives the same error.

**Proposed change:** use `MTL = NULL` and do not search unless the user passes an explicit `path`. Do the same for `loadImage(path = getwd())` and `cfmask(path = getwd())`, which are documented but surprising defaults.

### 8. MEDIUM. `read.WSdata()` assumes records are sorted and timestamps are unique
**Evidence:** R/water_WeatherStation.R:174-181 uses `datetime[1]` and `tail(datetime, 1)` to build the hourly sequence. Lines 222-229 use the "previous" and "next" records in data order. The daily summary combines `unique(date)` with `tapply` results, which come back sorted.

**Reproduction** (`r16_ws.R`, `r16b_ws.R`, INTA.csv)
- Reversed order (some loggers export newest first): "wrong sign in 'by' argument".
- Middle rows shuffled: **silently wrong `at.sat`**. Radiation 517.6 instead of 587.3 W/m2; temperature 24.97 instead of 25.31.
- A duplicated 11:00 record: "'-' only defined for equally-sized data frames".

**Proposed change:** after parsing, `order()` the records, drop exact duplicates with a warning, and interpolate `at.sat` with `approx()` on the sorted series.

### 9. MEDIUM. `calcH()` does not validate the anchors
**Reproduction** (`r12_calcH.R`, `r11_side.R`)
- A cold anchor on a masked (NA) cell: "missing value where TRUE/FALSE needed".
- Labels "Hot"/"Cold": "need finite 'ylim' values".
- Anchors outside the image: the same error.
- **Hot and cold swapped runs without any complaint**: a = −0.553 (dT decreases with Ts) and mean H 378 W/m2 instead of 264.
- Anchors given as sf in EPSG:4326 work (checked).

**Proposed change:** before iterating, check:
- that at least one hot and one cold anchor are present (with case-insensitive `type`);
- that the anchors are inside the raster;
- that Ts, Rn, G and Z.om are non-NA at the anchors;
- that the mean hot Ts is greater than the mean cold Ts.

Give clear errors.

### 10. MEDIUM. `par()` is changed and not restored (CRAN policy)
**Reproduction** (`r11_side.R`)
- `par(mfrow = c(2,2)); METRIC.EB(...)` leaves mfrow at `1 1`. R/water_METRIC.R:441 and :476 set it, with no `on.exit`.
- After an error inside METRIC.EB, mfrow stays at `1 2`.
- `plot(WS)` leaves `mar = 5.1 7.1 4.1 9.6` (R/water_WeatherStation.R:315).

**Proposed change:** `op <- par(...); on.exit(par(op), add = TRUE)`.

### 11. MEDIUM. Console output cannot be silenced, and plots are mandatory
**Reproduction** (`r11_side.R`, `r27_plot.R`)
- `METRIC.EB(verbose = FALSE)` inside `suppressMessages()` still prints 22 lines.
  - `relaxCriteria` prints with `print()` (R/water_sensibleHeatFlux.R:348), up to 400 lines.
  - `print(anchors)` and `print(paste0(...))` in METRIC.EB (:453, :460, :469).
  - `dailyET` prints "not implemented yet" (ET.R:260).
- METRIC.EB, METRIC.Rn, METRIC.G, calcH and ET24h always plot, and no function has a `plots` argument. METRIC.EB hard-codes `plots = T` for calcAnchors.
- `ET24h()` plots with `maxcell = 30 %` of the cells: **53.5 s on a full scene, versus 2.3 s** with terra's default.
- The calcH plot uses a fixed `xlim = c(0, 15)` while `maxit` defaults to 20.

**Proposed change**
- Use `message()` gated by `verbose` (or a `quiet` argument).
- Add `plots = interactive()`, or `plots = FALSE` by default, to the wrappers.
- Use terra's default `maxcell`.

### 12. MEDIUM. `loadImage()` does not mask Level-1 fill (DN = 0)
**Reproduction** (`r24_fill.R`): a 20 x 20 corner of the L8 Level-1 example set to 0.
- No cells become NA.
- A fill pixel gives TOA NIR = −0.125, **Ts = −53.6 K** and **Rn = 705.5 W/m2**.
- On full Level-1 scenes about a third of the raster is fill (19.6 M of 60.4 M cells in this scene). Those cells end up with plausible-looking flux values unless the user runs `cfmask()` on the DN image.

**Proposed change:** add `raw.image[raw.image == 0] <- NA` in `loadImage()` (R/water_RSfunctions.R:89-101), as `loadImageSR()` already does.

### 13. MEDIUM. u200 floor: values below 1 m/s are replaced by 4 m/s, which makes results jump
**Evidence:** R/water_sensibleHeatFlux.R:593-597.

**Reproduction** (`r13_u200.R`, same anchors). Mean LE:

| Station wind | u200 used | Mean LE |
|---|---|---|
| 0.50 m/s | 4.00 m/s | 207.6 W/m2 |
| 0.53 m/s | 1.03 m/s | **164.7 W/m2** |
| 0.55 m/s | | 187.0 W/m2 |
| 0.70 m/s | | 208.8 W/m2 |

A tiny increase in wind changes LE by 21 %. Combined with finding 3, the low-wind regime is ill-behaved.

**Proposed change:** clamp to the threshold, `u200 <- max(u200, u_min)`. Sam can confirm the value.

### 14. MEDIUM. Bands are selected by position in some functions and by name in others
**Evidence**
- `albedo()`, `LAI()`, `soilHeatFlux()` and `calcSR()` use layers `[[1:6]]` / `[[3]]` / `[[4]]` by position.
- `calcAnchors()`, `calcTOAr()` and `calcRadiance()` use names (`landsatBands()`).

**Reproduction** (`r28_bands.R`): the same named bands reordered as NIR, R, G, B, SWIR1, SWIR2.
- `LAI` mean drops from 1.254 to **0.001**.
- Albedo changes from 0.160 to 0.147.
- The NDVI that calcAnchors uses does not change.
- No warning in either case.

**Proposed change:** use `landsatBands()` everywhere, and stop if R or NIR is missing.

### 15. MEDIUM. `surfaceTemperature(method = "SC")` does not apply the method it documents
**Evidence:** R/water_netRadiation.G.R:591-602.
- `AF.j` (the Jiménez-Muñoz et al. 2014 water-vapour coefficients) is computed and then never used; codetools confirms "AF.j assigned but may not be used". The code comment also notes that ea is not water vapour content.
- The atmospheric functions actually used are METRIC's τNB, Rp = 0.91 and R_sky (Landsat 5/7 band-6 constants), applied to L8 band 10.
- `delta_e` is also unused.

**Proposed change:** either implement Jiménez-Muñoz (needs water vapour content in g/cm2) or document what is actually computed, and delete the dead code. **Flag for Sam.**

### 16. LOW. `getDEM()`: a tile that fails to download produces a partial DEM
**Evidence:** R/water_netRadiation.G.R:101-105 uses `tryCatch(..., error = function(e) NULL)`.

**Reproduction** (`r22_getDEM.R`, one tile URL forced to fail): the returned DEM has 5,800 of 13,400 cells NA (43 %). The only signal is a GDAL "HTTP 404" warning, which is the same signal an ocean tile produces. Network errors are treated like ocean tiles.

**Proposed change:** treat only a 404 as "no tile". Stop or warn on other errors, list the missing tiles, and report the NA fraction.

### 17. LOW. Return types and values
- `calcAnchors()` returns **NDVI as a character column**: `str()` shows `NDVI : chr "0.24" ...` (`r04_anchors.R`). R/water_sensibleHeatFlux.R:466-472 converts only columns 1-5; METRIC.EB works around it with `as.numeric`. All rows go through `c(..., "hot")` and so become character before being converted back. Build the data.frame column-wise instead.
- `flex` reports only the cold relaxation; `flex.hot` is computed but unused.
- The cold and hot follow-up loops use different thresholds (`<= 2` vs `< 2`).
- METRIC.Rn, METRIC.G, albedo and Z.om return layers with **empty names** (`r32_names.R`), because `saveLoadClean(stack.names = NULL)` clears them.
- `hourlyET()` (`r14_hourlyET.R`):
  - The default `height = 2` is never applied (ET.R:70 only assigns it when given), so a data.frame without `height` returns `numeric(0)`.
  - A waterWeatherStation without `at.sat` also returns `numeric(0)` silently.
- `dailyET(<data.frame>)` prints "not implemented yet" and then errors with "object 'ET.daily' not found". Use `stop()`.

### 18. MEDIUM. Tests miss the main paths where the bugs above live
**Evidence** (`r21_covr.R`, 68.2 % overall). Uncovered:
- `ET24h` (0 %, the final product)
- the plot methods
- the `saveLoadClean` write branch (finding 1)
- `dailyET(date=)` (finding 4)
- the whole **Landsat 7 path through METRIC.EB** (calcSR inside metricSurface, L7 Ts; only the vignettes run it, and `--ignore-vignettes` skips them)
- `plain = FALSE` / mountainous corrections
- Perrier Z.om
- `calcH(method = "lm")`
- calcAnchors "best" and "random" with n > 1, and the deprecated names

Separately, `test-energy-balance.R` runs METRIC.EB 12 times (193 s).

**Proposed change**
- Add small tests for the paths above, including regression tests for findings 1-4 and 8.
- Make the tests cheaper: run calcH on a 20 x 20 window, and skip the heaviest tests unless `NOT_CRAN` is set.

### 19. MEDIUM. CRAN resubmission readiness
- About 24 exported functions have **no `\examples`**: LAI, albedo, calcAnchors, calcH, cfmask, ET24h, loadImage, surfaceTemperature, soilHeatFlux, METRICtopo, solarAngles, calcTOAr/SR/Radiance, waterOptions and others. CRAN usually asks for them on a new submission.
- `cran-comments.md` still describes version 0.3.2.
- The tests write `aoi` to `.GlobalEnv` (tests/test-inputs.R:205). That is tolerated, but avoidable once finding 6 is fixed.
- `waterOptions()` contains a dead `save` block (R/zzz.R:158-176). If it were ever enabled, it would write a `rasterOptions_` file into the startup directory, which CRAN forbids.

### 20. LOW. API consistency
**Argument order**
- `METRIC.EB(image.DN, image.SR, WeatherStation, ...)` versus `METRIC.Rn` and `METRIC.G(image.DN, WeatherStation, MTL, ...)`. An L7 positional call `METRIC.EB(DN, WS, MTL)` passes WS as `image.SR`.
- `LAI(method, image, ...)`: `LAI(sr)` errors with "[==] Can only match character values".
- `calcH(anchors, method, Ts, ...)`.
- `albedo(image.SR, aoi, coeff, ...)`: `albedo(sr, "Olmedo")` errors with "[crop] cannot get a SpatExtent".

**Unused arguments and options**
- An unused `aoi` argument sits early in the argument list of LAI, calcTOAr, calcSR, soilHeatFlux, surfaceTemperature and calcAnchors.
- `verbose` exists only in METRIC.EB.
- calcH's `method` and `maxit` cannot be reached from METRIC.EB.

**Defaults**
- `read.WSdata` defaults to `date.format = "%Y-%m-%d"`, but `read.WSdata2` to `"%d/%m/%Y"`, and this is undocumented.
- The default `alb.coeff = "Tasumi"` uses Landsat 7 coefficients, while L8/L9 is the primary path and the package's own docs recommend "Olmedo" (an `"auto"` option would help).
- `sat = "auto"` means two different things: layer names in `albedo()`, MTL or files in getwd() elsewhere.

**Other**
- `aoi` only crops to the bounding box; it does not mask to the polygon.
- No input checks `compareGeom()`: a DEM one cell larger or shifted half a cell fails with "[*] extents do not match" after the first steps (`r30_grid.R`).

### 21. LOW. Dead, duplicated and commented-out code
**Unused variables and arguments** (codetools)
- `sun.elevation`, `sun.azimuth` and a duplicated SCENE_CENTER_TIME parse in `solarAngles`
- `d2` and `bands` in `calcTOAr`; `bands` in `calcSR`
- `AF.j` and `delta_e` (finding 15)
- `latent.heat.vaporization` in calcH (LE.cold uses 0.002361 while ET24h uses 0.00236)
- `flex.hot`
- `r` in waterOptions
- `DOY` in dailyET

**Dead assignments**
- `tauD <- rast(tauB)` (R/water_netRadiation.G.R:502)
- `G <- rast(NDVI)` (:758)
- `calcSR`'s fallback `calcTOAr()` with no arguments (R/water_RSfunctions.R:328) can never work.

**Duplicated code**
- 6 copies of the "next anchor" loop (about 140 lines; R/water_sensibleHeatFlux.R:179-436).
- The calcH "mean" and "lm" loops (about 90 lines each).
- 4 copies of the MTL date/time regex parsing.
- METRIC.EB:378-396 repeats metricSurface:146-169 (the image.SR check and thermal-layer detection).
- The NDVI and ea formulas are repeated across functions.

**Commented-out and redundant code**
- 28 commented-out code blocks, e.g. zzz.R:16-19, 52-60, 109-124; METRIC.EB `#setTxtProgressBar`; the old `cfmask` at RSfunctions.R:434-445.
- `try(..., silent = TRUE)` around the anchor sampling hides errors.

### 22. LOW. Style and messages
- 13 assignments use `=`; `T`/`F` is used 28 times (tools::checkTnF); 81 lines are over 100 characters; 21 `if()` conditions use `&`/`|`; 5 use `1:length()`.
- Naming mixes dot.case, camelCase and snake_case, sometimes inside one function (calcH: `cell_vals`, `r.ah.hot.previous`).
- terra calls are mixed: 60 qualified with `terra::` and about 200 unqualified.
- `import(sf)` is unnecessary because every sf call is `sf::`-qualified.
- The exported name is misspelled: `SWtrasmisivity` (add `SWtransmissivity` as an alias).
- Some messages contain multi-line string literals with embedded indentation: R/water_sensibleHeatFlux.R:127-128 and R/water_WeatherStation.R:73-75 and :84-85.
- `paste()` produces double spaces: "I can only find  1  anchors".
- `prepareSRTMdata`'s error tells users to save tiles in the "working directory" although `path` is an argument.
- Use `exp()` instead of `e <- 2.71828`, and use 273.15 consistently (hourlyET uses 273.16).

---

## Checked and found not to be bugs
- terra `[`/`[<-` with NA masks in soilHeatFlux and SWtrasmisivity keeps cells aligned (`r02_index.R`).
- The weather-station marker in the calcAnchors plot works (`points.sf` exists).
- `cfmask()` with no clouds and `buffer > 0` works.
- Anchors given as sf in another CRS work (`r33_anccrs.R`).

---

## Top 5 priorities
1. **Finding 1:** fix `saveLoadClean()` file naming. writeResults = TRUE currently corrupts results and silently turns Perrier runs into short-crops runs.
2. **Finding 2:** validate all method/coefficient arguments with `match.arg()` up front. Typos currently return stale objects from the user's workspace or fail after minutes of work.
3. **Findings 3, 9, 13:** make `calcH()` robust. Damp the r_ah iteration, validate anchors (including hot/cold swaps), fix `maxit`, clamp u200 at the threshold, and report convergence.
4. **Finding 4 (with 8):** time handling in weather data. Have `dailyET()` select the local overpass day (correct east of UTC and for explicit `date`), and sort and de-duplicate records in `read.WSdata()`.
5. **Finding 5:** full-scene performance (about 58 min and 33 GB for one scene now). calcH should iterate on the anchors and then do one `lapp()` pass; calcAnchors should read its inputs once into vectors and use coordinate distances; replace constant rasters with scalars. Then add tests for the uncovered main paths (finding 18) and examples for CRAN (finding 19).
