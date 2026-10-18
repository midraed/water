# Review of **water** 0.8.5 (master 70027b3) — usability for irrigation practice

Reviewer: **R** (agricultural extensionist; agronomy, crop water use, basic R).

## What I did

I worked as a new user would, with R 4.5.3 and the installed master (`lib_new`):

- I read README, DESCRIPTION, NEWS 0.8.5, the three vignettes and the two Spanish download guides, and about 30 help pages as text.
- I ran the vignette code with `knitr::purl()` + `source()`: Landsat8.Rmd, then METRIC_simple.Rmd, then METRIC_advanced.Rmd.
- I tried about 40 common mistakes with the example data and the full scene LC08_L2SP_232083_20160209 cropped with `createAoi()`. For each one I checked how far the error travels into the final daily ET map.
- I timed runs from the example crop up to the full scene.

All scripts and outputs are in `<scratch>` (files `s01_*.R` to `s27_*.R`, `timing_*.txt`).

Reference values for the Landsat 8 example (INTA station, 9 Feb 2016):

| Quantity | Value |
|---|---|
| Daily ETr | 5.34 mm/day |
| Daily ETo | 4.48 mm/day |
| Hourly ETr at the overpass | 0.50 mm/h |
| Mean ET24 (vignette run: cfmask, n = 5) | 3.98 mm/day |

Most "mistake" results below are compared against these values.

Overall: the package runs end to end on the examples. Many of the error messages added in 0.8.5 are clear, for example for a wrong date format, an unknown time zone, an MTL file that does not exist, MTL contents passed instead of the file name, Level-1 versus Level-2 folders, and a missing `surf.temp`. But many realistic mistakes still pass **silently** and change daily ET by 30–100% or more. The package also gives the practitioner almost no help to judge whether a map can be trusted, or to turn it into numbers per field.

---

## Findings

### A. Getting started

**1. The README says water is on CRAN, but it was archived in 2021. The install line does not build the vignettes. — High**

*Evidence.*
- README.md:8 says "water package is available at [CRAN]" and shows a CRAN downloads badge. `available.packages(repos = "https://cloud.r-project.org")` does not list water. The CRAN page says: "Archived on 2021-10-05 as check problems were not corrected in time" (checked with `s18b_cran.R`).
- README.md:12-13 installs with `remotes::install_github("midraed/water")`. `remotes::install_github()` defaults to `build_vignettes = FALSE` and `--no-build-vignettes` (remotes reference manual on CRAN), but README.md:23 then sends the user to `vignette(package = "water")`.

*Why it matters.* `install.packages("water")` is the first thing a new user types, and it fails with "package 'water' is not available". After the GitHub install, the vignettes the README points to are not there.

*Proposed change.*
- Write: "water is not on CRAN at the moment (archived in 2021). Install the current version from GitHub: `remotes::install_github("midraed/water", build_vignettes = TRUE)` (needs pandoc, included in RStudio), or read the vignettes online at https://midraed.github.io/water/."
- Remove the CRAN badge, or label it as historical.

**2. There is no single entry point: no reading order, no quick start, no `?water`, and a search for "evapotranspiration" finds nothing. — Medium**

*Evidence.*
- README.md:23-26 lists "METRIC simple" and "METRIC advanced" first. Both are Landsat 7, and Landsat 7 is no longer the main use case.
- Landsat8.Rmd is the one vignette for current data, and it never explains METRIC. The equations (LE = Rn − G − H, the 24 h scaling) are only in the Landsat 7 vignettes. The anchor pixels get one sentence, only in METRIC_advanced.Rmd:181. No vignette uses the terms "ETrF", "crop coefficient" or "Kc".
- There is no `water-package.Rd`, so `?water` returns nothing.
- `help.search("evapotranspiration", package = "water")`, `"crop coefficient"`, `"irrigation"` and `"Kc"` all return no topics (`s02_search.R`).
- The README has no code example.

*Why it matters.* A Landsat 8/9 user has to read three documents, written for two different satellites, to put one workflow together.

*Proposed change.*
- In the README, add "Start here: vignette 'Landsat 8 and 9'", followed by a 10-line quick start: `loadImageSR` → `cfmask` → `read.WSdata` → `METRIC.EB` → `dailyET` → `ET24h` → `writeRaster`.
- Move the METRIC explanation (energy balance, anchors, ETrF, ET24) into the Landsat 8/9 vignette, and keep the Landsat 7 vignettes as "advanced / step by step".
- Add a package help page (roxygen `"_PACKAGE"`) with a map of the workflow, and `@concept evapotranspiration` / `@concept crop coefficient` on the main topics.

**3. The Spanish download guides are not installed with the package, and the rest of the documentation is only in English. — Medium**

*Evidence.*
- `.Rbuildignore:7` (`^vignettes/.*\.qmd$`) excludes both guides from the package. They can only be reached through GitHub links in Landsat8.Rmd:16.
- The guides are in Spanish. The README, vignettes and help pages are in English only, and there is no Spanish quick start.
- The guides need Quarto (`quarto render`, Landsat_download_USGS.qmd:51-53). Quarto is not part of an R installation.

*Why it matters.* The intended users (Mendoza and Talca examples; INTA and CITRA stations) are mostly Spanish-speaking extension staff. Today they get Spanish instructions for downloading, then English for everything else. Someone working offline after installing the package gets no download instructions at all.

*Proposed change.*
- Publish both guides as pkgdown articles and link them from README and `?loadImageSR`.
- Add a Spanish version of the Landsat 8/9 vignette, or at least a Spanish "Guía rápida". The vignette is short, so this is cheap.
- Say that the Planetary Computer guide can also be run in plain R with `knitr::purl()` + `source()`, without Quarto.

**4. The download guides have developer defaults, download whole scenes, and stop before the hand-off to water. — Medium**

*Evidence (from reading the code; Quarto is not installed here).* The Planetary Computer search itself works: I reproduced the guide's id logic with a direct STAC query (`s07_net.R`).
- **USGS guide, `ids`:** the default is `ids: "LC08_L2SP_232083_..."` (Landsat_download_USGS.qmd:15), while step 1 says "manteniendo `ids: ""`" (line 33). With the user's own polygon and dates, line 246 stops ("Elegir IDs de la tabla de esta búsqueda."), and with `execute: error: false` (line 9) the render aborts, so the user never sees the table.
- **Destination folders:** the USGS guide uses `destino: "../inst/extdata/L8"`, a path inside the package sources. The Planetary Computer guide uses `destino <- "../../L8"` (line 50), two folders above the file.
- **Download size:** the Planetary Computer guide downloads every `.TIF` of the scene (line 143): 19 files (Landsat8.Rmd:23: "A full scene is almost 1 GB"). water reads only 8 of them (SR_B2–B7, ST_B10, QA_PIXEL) plus the MTL, and those already take 639 MB here. The USGS guide keeps the TAR *and* the extracted files.
- **Cloud cover:** both guides only show `eo:cloud_cover` for the whole scene (they say so).
- **Hand-off:** neither guide ends with the call to water. The Planetary Computer guide ends with manual scaling (`b4 * 0.0000275 - 0.2`, line 222), which `loadImageSR()` already does.

*Why it matters.* Rural extension offices have slow connections and little disk space. A scene can be "7% cloudy" and still have clouds over the valley you care about.

*Proposed change.*
- Defaults: `ids: ""` and `destino: "landsat"` (next to the .qmd).
- Planetary Computer: read only the 9 needed files, cropped to the AOI with `/vsicurl/` + `crop()`, as `getDEM()` already does. For a 30 × 30 km AOI this is tens of MB instead of a whole scene.
- Report the clear fraction *inside the AOI* from QA_PIXEL.
- End both guides with:

  ```r
  library(water)
  sr <- loadImageSR(carpeta, aoi = aoi)
  sr <- cfmask(carpeta, image = sr)
  ```

  plus a note: "no aplicar los factores de escala a mano: loadImageSR() ya lo hace".

### B. Weather station input (read.WSdata, dailyET)

**5. The names in `columns` are ignored: the same columns written in file order give 587 °C and an ETr of 408 mm/day, with no warning. — High**

*Evidence.* `s24_columns.R`, INTA.csv. The columns in the file are datetime, temp, RH, pp, radiation, wind.
- Writing `columns = c(date=1, time=1, temp=2, RH=3, rain=4, radiation=5, wind=6)`, the same assignment in the file's order, gives:
  - at the overpass: radiation 25.3 W/m², wind 58.3 m/s, RH 0, temp 587.3 °C;
  - `dailyET()` = 407.7 mm/day (correct value: 5.34).
- `?read.WSdata` does say "Names in this vector are ignored", but the Usage line shows a *named* vector, which invites users to name the columns.
- Giving column names as strings (`radiation = "radiation"`) instead produces a confusing deprecation warning about "the old syntax", then "undefined columns selected".

*Why it matters.* This is the most natural way to describe your own logger file. The result is nonsense and nothing stops it.

*Proposed change.*
- When `columns` has names, match by name (`columns[c("date","time","radiation","wind","RH","temp","rain")]`).
- Also accept column *names* of the CSV (`columns = c(datetime = "fecha", temp = "temp", ...)`).
- Stop with a clear message if a required role is missing.

**6. Weather data are not checked for plausibility or units, so common unit mistakes silently change daily ETr by −73% to +111%. — High**

*Evidence (`s09`, `s10`, `s23`).* None of these cases gives an error or a warning:

| Mistake | Daily ETr (correct: 5.34 mm/day) | Effect on mean ET24 |
|---|---|---|
| Wind logged in km/h, `cf` not set | 7.45 | 3.98 → 5.24 mm/day (+32%) |
| Radiation in MJ m⁻² h⁻¹ | 1.43 | |
| Temperature in °F | 11.25 | |
| RH as a fraction (0–1) | 6.31 | |
| Radiation 1500 W/m² at the overpass (above clear sky) | 6.69 | |
| RH and temp columns swapped (temp = 58 °C) | 8.38 | METRIC.EB later crashes (finding 10) |

On top of that:
- `print(WeatherStation)` shows no units (water_WeatherStation.R:371-384).
- `$daily$radiation_sum` is a sum of W/m² readings, which has no physical meaning: 5663 for hourly records, four times more for 15-minute records.
- `cf` is positional and covers only radiation, wind and temperature. `?read.WSdata` Details documents only the km/h case.

*Why it matters.* The maintainer's own example had exactly this km/h bug (#72, ET of the vignettes 45% too high). Every station network logs in different units.

*Proposed change.*
- Add range checks in `read.WSdata()`, with messages that say what to do:
  - wind > 20 m/s → "is wind in km/h? use cf = c(1, 0.2777778, 1)";
  - daytime maximum radiation < 50 → "is radiation in MJ/m2? convert to W/m2";
  - temp outside −30…55 °C, or RH outside 0–100 (or max RH ≤ 1) → error;
  - Rs/Rso > 1.1 at midday, using the clear-sky Rso that `hourlyET()` already computes (ASCE-EWRI 2005, App. D).
- Print units in the header: `radiation (W/m2) wind (m/s) RH (%) temp (C)`.
- Replace `cf` with a named `units = c(radiation = "W/m2", wind = "km/h", temp = "C")`.

**7. A wrong time zone is accepted silently, and the overpass time is never checked against a Landsat morning. — Medium**

*Evidence (`s09`).*
- `tz = "UTC"` for a logger on local time: `at.sat` is taken at 14:27 local, giving radiation 788.9 W/m² instead of 587.3, temp 27.5 °C instead of 25.3, and hourly ETr 0.72 mm/h instead of 0.50 (+44%). No warning. Mean ET24 changes only −4% (`s10`) because of the internal calibration, but every condition printed "at satellite overpass" is wrong.
- `tz = "Etc/GMT-3"` (the POSIX sign trap) puts the overpass at 17:27 local.
- On the positive side, `tz = ""` warns and names the computer's time zone, and a misspelled zone is a clear error.

*Why it matters.* Station clocks in UTC, in local standard time, or with daylight saving time are a constant source of confusion in the field.

*Proposed change.*
- In `read.WSdata()`, warn when the overpass converted to the station clock falls outside about 09:00–12:30 local solar time. For example: "The Landsat overpass is 14:27 in the station clock (tz = 'UTC'); Landsat passes around 10:00–11:30 local time. Check tz."
- Always print the overpass time in both UTC and the station time.

**8. Incomplete or gappy days are summed without warning: records up to 13:00 give daily ETr 2.04 instead of 5.34 mm/day. — High**

*Evidence (`s09`, `s10`).*
- `dailyET()` sums whatever hourly values exist (`sum(ET.daily, na.rm = TRUE)`, water_ET.R:261).
- A file that ends at 13:00 (logger downloaded at noon) gives ETr 2.04 mm/day and mean ET24 1.52 mm/day instead of 3.98 (−62%), silently.
- Four missing midday hours are interpolated linearly: ETr 4.76 (−11%). Four `NA` midday radiation values: 4.83 (−10%).
- On the positive side, records that end *before* the overpass do give a clear warning.

*Why it matters.* ET24 is directly proportional to daily ETr. Partial days are common: a field download, a power cut, a rejected sensor reading.

*Proposed change.*
- Count the hours that are really measured on the image day. If fewer than 24, or a gap is longer than 2 h, warn or stop. For example: "Only 14 of 24 hours on 2016-02-09 (last record 13:00): daily ETr would be underestimated. Provide the whole day or ETr.daily from another source."
- Stop dropping missing hours silently: today `na.rm = TRUE` is hard-coded at water_ET.R:261.

**9. With a multi-day station file, `dailyET()` without `MTL` returns the ETr of the first day, not the image day. — High**

*Evidence (`s21_multiday.R`).*
- Three days of records (8–10 Feb), read with `MTL`; `WS$at.sat` correctly shows 9 Feb.
- `dailyET(WS)` returns 3.14 mm/day, the ETr of 8 Feb. The correct value is 5.34, which is what `dailyET(WS, MTL = MTL)` gives. The working directory had no MTL file.
- `dailyET(WS, date = "2016-02-09")` without MTL also returns 3.14: the help says `date` "is ignored" without MTL.
- A `date` that is not in the records returns 0 (`"2016-03-01"`) or 0.03 (`"2016-02-10"`), silently (`s10`).
- `date = "2016-02-09"` *with* MTL gives 5.31 instead of 5.34, because that branch builds days with `as.Date()` in UTC (water_ET.R:242).

*Why it matters.* Practitioners download weeks or months of station data. Forgetting `MTL` in the second call, while the station object already knows the overpass date, gives an ET map 41% too low with no message.

*Proposed change.*
- `dailyET()` should default to the date of `WeatherStation$at.sat` when present.
- Honour `date` without MTL, using the station time zone.
- Stop when the requested day has no records.
- Warn when the records cover more than one day and neither MTL, at.sat nor date tells which day to use.

**10. Wrong inputs end in a cryptic crash: "missing value where TRUE/FALSE needed". — Medium**

*Evidence (`s11_long.R`).*
- Station longitude typed without the minus sign (`long = 68.86469`): only the warning "No pixel is within WSbuffer (30000 m) of the weather station", then `Rn = NaN` at the anchors. `calcH` fails at water_sensibleHeatFlux.R:733 with `Error in if (abs(delta.r.ah.hot) < 1 & ...): missing value where TRUE/FALSE needed`.
- RH and temp swapped (temp 58 °C): the H iteration diverges (a = 2.4e126) and crashes with the same message.

*Why it matters.* The user cannot tell that the cause is a weather input.

*Proposed change.*
- In `read.WSdata()` / `METRIC.EB()`, compare the station position with the scene corners in the MTL. For example: "The weather station (lat −33.01, long 68.86) is about 11,400 km from the image. Check the signs: negative for south and west."
- In `calcH()`, stop as soon as `a`, `b`, `r.ah`, Rn or G at the anchors are not finite, and print those values: "The H iteration diverged: check the weather data at the overpass (temp, RH, wind) and the anchors".

**11. A semicolon / decimal-comma CSV gives "duplicate 'row.names' are not allowed". — Medium**

*Evidence (`s09` case 6).*
- A file written with `write.csv2()` and read with `read.WSdata()` fails inside `read.csv()` with "duplicate 'row.names' are not allowed".
- The package's own hint ("If the file uses decimal commas, use read.WSdata2()", water_WeatherStation.R:89-91) is never reached.
- `read.WSdata2()` reads the same file correctly.

*Why it matters.* This is the default CSV format of Excel in Spanish locales, so most practitioners in Chile and Argentina will hit it first.

*Proposed change.*
- Wrap `read.csv()` in `tryCatch()`, or sniff the first line for `;`, and say: "The file seems to use ';' and decimal commas: use read.WSdata2() (or sep = ';', dec = ',')."

**12. The u200 fallback warning is cryptic, and one vignette hides it. Both Talca vignettes now trigger it. — Medium**

*Evidence.*
- With the corrected km/h wind, the Talca station has 0.31 m/s at the overpass. Both METRIC_simple and METRIC_advanced give: "u200 less than threshold value = 0.5903m/s. using u200 = 4m/s" (water_sensibleHeatFlux.R:593-597).
- METRIC_simple hides it with `warning=FALSE` (METRIC_simple.Rmd:107, also line 130).
- Allen et al. (2007) do not mention a 4 m/s replacement. It is a package rule, and it multiplies the wind about 7 times.

*Why it matters.* The flagship example runs on a wind speed that is not the measured one, and the reader is not told.

*Proposed change.*
- Reword the warning: "Wind at the station at the overpass is 0.31 m/s, so the wind extrapolated to 200 m (u200 = 0.59 m/s) is below the 1 m/s minimum for the H iteration; u200 = 4 m/s is used instead. Check the wind units (cf) and the sensor; results are sensitive to this value."
- Document the rule in `?calcH`, and make the threshold and the replacement value arguments.
- Remove `warning=FALSE` from the vignette and explain the warning there.

**13. Surface temperature in °C is accepted and almost doubles ET. — Medium**

*Evidence (`s14` u2).*
- `surf.temp = image.SR$Band_10 - 273.15` runs without complaint. The anchor table shows Ts 43.8 and 35.0, and mean ET24 goes from 3.62 to 6.73 mm/day (+86%).

*Why it matters.* Surface temperature exported from other tools (GEE, QGIS) is often in °C.

*Proposed change.*
- In `METRIC.EB()` / `METRIC.Rn()` / `METRIC.G()`: `if (global(surf.temp, "max", na.rm = TRUE) < 200) stop("surf.temp must be in Kelvin (values look like Celsius: add 273.15)")`.

### C. Image, area of interest and calibration

**14. A hidden `aoi` from the workspace is used: the Landsat 8 vignette fails at its first line after METRIC_simple. — High**

*Evidence (`s15_aoitrap.R`).*
- `waterAutoAoi = TRUE` by default (zzz.R:9). `loadImage()`, `loadImageSR()` and `albedo()` crop to any object called `aoi` in `.GlobalEnv` (water_internals.R:117-129).
- The METRIC_simple and METRIC_advanced vignettes create `aoi <- createAoi(...)` (Talca, EPSG 32719) at lines 59 and 61, and never use it.
- Running the Landsat 8 vignette afterwards in the same session: `loadImageSR(path = scene_folder)` fails with "[crop] extents do not overlap". Nothing mentions the `aoi` object.
- With an `aoi` that partly overlaps, the image is silently cropped (100 × 67 instead of 134 × 184 pixels).
- With an `aoi` of the user's farm in the workspace, a call meant to load the whole scene loads only the farm, which leads straight to finding 16.

*Why it matters.* `aoi` is the name every example teaches. The failure happens at the first step and does not point to the cause.

*Proposed change.*
- Set `autoAoi = FALSE` by default, or at least `message("Cropping to the object 'aoi' found in your workspace (waterOptions(autoAoi = FALSE) to disable)")`.
- In the Landsat 7 vignettes, either pass `aoi` explicitly or remove the unused `createAoi()` chunk.

**15. Nothing warns when the clouds were not masked: a cold anchor landed on a cloud pixel. — Medium**

*Evidence (`s13_clouds.R`).*
- AOI of 40 × 24 km around the station; 22.6% of the area is flagged as cloud or shadow in QA_PIXEL; the MTL says CLOUD_COVER = 7.55.
- Without `cfmask()` and with `n = 5`, the coldest "cold" anchor (Ts 301.2 K) is a cloud pixel. Mean ET24 on clear pixels is 0.19 mm/day lower than with the mask.
- With the default `n = 1`, the same cold pixel was chosen in both runs. The safety of skipping the mask depends on luck.

*Why it matters.* `cfmask()` is a separate, easy-to-forget step, and the user is not told that it is mandatory.

*Proposed change.*
- In `METRIC.EB()`: read `CLOUD_COVER` from the MTL, and warn when it is > 0 and the image has no masked pixels: "The scene has 7.6% clouds but image.SR has no masked pixels: run cfmask() first".
- Alternatively, mask by default in `loadImageSR(mask = TRUE)`.

**16. An AOI the size of a farm runs without complaint and underestimates ET by about 40%. — High**

*Evidence (`s14` q, `s25_smallaoi.R`).*
- A 20 × 20 pixel AOI (600 m × 600 m, the size of a farm around the station) runs normally. The anchors are 1.0 K apart (cold 309.1 K, hot 310.2 K). The "flexible" method relaxes the criteria and replaces ETp.coef with 0.892.
- For **the same 400 pixels**, mean ET24 is:

  | AOI | Mean ET24 of the 400 pixels |
  |---|---|
  | Farm-sized (600 m) | 2.87 mm/day |
  | 3 km | 5.01 mm/day (hot–cold ΔT 7.7 K) |
  | 10 km | 4.60 mm/day (ΔT 12.8 K) |
  | 30 km | 4.77 mm/day (ΔT 23.3 K) |

- No warning at all.

*Why it matters.* Cropping to "my farm" is the first thing a practitioner will do, and it gives a believable but wrong map. The docs say nothing about the minimum area. Landsat8.Rmd:23 only says to crop "an area of interest".

*Proposed change.*
- Stop, or warn strongly, when the hot–cold Ts difference is below a threshold (for example 5–10 K), and when the AOI is smaller than a few km² or has fewer than N candidates.
- Document in `?METRIC.EB` (aoi) and in the vignette: "Run METRIC on a large area (tens of km) that contains both well-watered full-cover fields and dry bare soil, then crop the *results* to your fields."

**17. `n` has a default of 1, while every example uses 5, and the choice changes ET by up to 10%. Nothing explains it. — Medium**

*Evidence (`s14` n).*
- `METRIC.EB(n = 1)` is the default (water_METRIC.R:374). The vignettes and examples all use `n = 5`. `?METRIC.EB` only says "number of pairs of anchor pixels to select".
- On the example crop, mean ET24 is 3.62 / 3.86 / 3.98 / 3.71 mm/day for n = 1 / 3 / 5 / 10.
- With `n = 10` the function warns "I can only find  6  anchors with cold pixel conditions".
- With `n = 5`, the five cold anchors span up to 5 K (302.7–307.7 K in `s13`), because they are the 1st to 5th coldest candidates at least 500 m apart.

*Why it matters.* Users copy the examples without knowing what they trade off.

*Proposed change.* Explain `n` in `?METRIC.EB` and in the vignette ("more pairs average out noisy pixels but include warmer or less ideal fields"), and use the same default as the examples.

**18. Many pixels have ETrF above 1.05, and the package gives no quality summary. — Medium**

*Evidence (`s05_etrf.R`, `s08`).*
- In the Landsat 8 vignette run, ETrF (ET24 / ETr) has median 0.80, 95th percentile 1.17 and maximum 1.57. 16.6% of the pixels are above 1.05, 3.4% above 1.2 and 0.9% above 1.3.
- The Talca vignettes: 6.9% and 6.1% of pixels above 1.1, maxima 1.38 and 1.56.
- Allen et al. (2007, after Eq. 54): "One should generally expect ETrF values to range from 0 to about 1.0 … perhaps up to 1.1 if it has been recently wetted".
- Also in the Landsat 8 run: SensibleHeat reaches 1031 W/m², which is more than Rn − G, and 4.8% of pixels have ET24 = 0.

*Why it matters.* An advisor would read ETrF = 1.3–1.5 as "this vineyard uses more than alfalfa" and over-irrigate.

*Proposed change.*
- Have `METRIC.EB()` / `ET24h()` print a short QC summary: share of pixels with ETrF > 1.1 and < 0, the anchor Ts, and the hot–cold ΔT.
- Optionally return a QC flag layer.
- Add a histogram of ETrF to the vignette, with an explanation of what values outside 0–1.1 mean (water, wet soil, anchor problems).

**19. There is no guidance or tool for checking the anchor pixels. — Medium**

*Evidence.*
- `calcAnchors(plots = TRUE)` plots the anchors over LAI only (water_sensibleHeatFlux.R:457-462). There is no RGB background and no table with albedo or distance to the station.
- In the 40 × 24 km and 30 × 30 km AOIs, the hottest candidate (330.3 K, NDVI 0.12) was about 7 K hotter than the next hot candidates (322–323 K) (`s13`). It is a mixed pixel at the edge of a non-vegetated surface (NDVI −0.10 to 0.01, Ts 331–332 K; 5 × 5 window in `s26_hotpixel.R`), probably not a dry agricultural field. With the default `n = 1` it is the only hot pixel used (`s25`).
- Mean ET24 of the same fields varies between 4.60 and 5.01 mm/day depending on the AOI size, because different anchors are chosen (`s25`).
- Allen et al. (2007, Summary): METRIC "is highly reliant on the operator's ability to select appropriate cold and hot pixels. It is currently not a hands-off model."
- The default is fully automatic, with criteria relaxed by "flexible".

*Why it matters.* The step the paper says needs expert judgement is the one the practitioner gets no help with.

*Proposed change.*
- Add a vignette section "Checking the anchors":
  - plot the anchors on an RGB composite (`plotRGB(image.SR, 4, 3, 2)`) and export them with `sf::st_write(Energy.Balance$anchors, "anchors.gpkg")` to look at them in Google Earth or QGIS;
  - check that the cold pixels are in large, uniform, well-watered fields (not tree rows, ponds or cloud edges) and the hot pixels in dry bare fields (not roads or roofs);
  - re-run with `anchors =` if needed.
- Add a helper `checkAnchors(Energy.Balance)` that prints per-anchor albedo, NDVI, Ts, distance to the station and the standard deviation of Ts in the 3 × 3 neighbourhood.

**20. The console output is noisy and partly mislabelled, and `verbose = FALSE` does not silence it. — Low**

*Evidence (`s18`, `s11`).*
- With `verbose = FALSE` + `suppressMessages()` + `suppressWarnings()`, METRIC.EB still prints 12 lines: "relaxing criteria … 1%" … "11%" and the anchor table. These come from `print()` calls at water_sensibleHeatFlux.R:348 and water_METRIC.R:453/460/469.
- The relaxation warning shows impossible or empty limits: "LAI 2.67 to 6.66" (LAI is capped at 6), "Z.om NA to 0.005", "Ts 322.0862 to NA".
- The `flex` column repeats the cold-pixel value in the hot rows (water_sensibleHeatFlux.R:473-476).
- With `verbose = TRUE`, the cold and hot tables have wrong column names: "LAI NIR V1 V1.1 Band_10" instead of LAI, NDVI, albedo, Z.om, Ts (lines 447-454).
- METRIC.EB also draws two plots as a side effect.

*Why it matters.* When you process 20 dates in a loop, real warnings get lost in the noise.

*Proposed change.*
- Use `message()` once ("Cold-pixel criteria relaxed by 11%: LAI 2.67–6, albedo 0.160–0.278, Z.om 0.027–0.089"), only when `verbose`.
- Name the columns correctly, put `flex` only on the cold rows (or add `flex.hot`), and add `plots = TRUE/FALSE` to METRIC.EB.

### D. Names, arguments and help pages

**21. Several argument names are specialist jargon, the defaults do not match the examples, and the vignettes contradict each other on `Z.om.ws`. — Medium**

*Evidence.*
- **`ETp.coef`:** "p" stands for the outdated "potential ET". It is really the ETrF (alfalfa-basis Kc) of the cold pixel. When "flexible" changes it, the user reads "ETp.coef updated to: 0.892 (Kc = 1.02 …)".
- **`Z.om.ws`:** the help text is good (0.123 × vegetation height). But METRIC_simple uses the default 0.015 for the Talca orchard station, while METRIC_advanced uses 0.03 for the same station (METRIC_advanced.Rmd:196). Landsat8.Rmd:64 also uses 0.03 for INTA. None of them says why.
- **`plain`:** reads as "plain output"; the meaning is flat terrain.
- **`extraParameters`:** hides the practical knobs `WSbuffer`, `minDist` and `deltaTemp` inside a named vector.
- **`cf`:** positional (see finding 6).
- **`read.WSdata(height = 2.2)`:** this default is the height of the example station, while the standard is 2 m and `hourlyET()` uses 2.
- **`METRIC.EB` title:** "Estimates Energy Balance using METRIC2010 Model", while everything else cites METRIC 2007.

*Why it matters.* Each of these names makes a non-specialist stop and guess.

*Proposed change.*
- Keep the old names as aliases, and add or document clearer forms:
  - `ETrF.cold` (instead of `ETp.coef`);
  - `ws.vegetation.height` in m, from which `Z.om.ws = 0.123 * h` is computed;
  - `flat` (instead of `plain`);
  - `anchors.control = list(WSbuffer, minDist, deltaTemp)`.
- Make `height` required, or default it to 2.
- Use one `Z.om.ws` value in all vignettes and say why.
- Change the title to "… with the METRIC model (Allen et al. 2007)".

**22. `ET24h()` needs four separate layers, drops the ETrF it computes, and has no example. — Medium**

*Evidence.*
- Every vignette calls `ET24h(Rn = Energy.Balance$EB$NetRadiation, G = …$SoilHeat, H = …$SensibleHeat, Ts = …$surfaceTemperature, WeatherStation, ETr.daily)`.
- ETrF is computed at water_ET.R:165, but only ET24 is returned. The layer has no name and no unit (`names(ET24) == ""`, `s17`).
- `?ET24h`, `?calcAnchors`, `?calcH` and `?cfmask` have no Examples section.
- There is no single call from the inputs to daily ET.

*Why it matters.* ETrF, the crop coefficient, is *the* number an irrigation advisor uses. ET24h is the last step of the workflow and, like `cfmask()`, it has no example.

*Proposed change.*
- Add `ET24h(Energy.Balance, ETr.daily = …)` for `waterLSEB` objects, returning `c(ET24 = …, ETrF = …)` with `units()` "mm/day" and "-".
- Copy the vignette code into the Examples of `?ET24h`.
- Optionally add a wrapper `METRIC.ET()` that returns the energy balance, ETrF and ET24 together.

**23. The help pages say what the outputs are, not what they mean or which values are reasonable. — Medium**

*Evidence.*
- `?METRIC.EB` Value lists "NetRadiation, SoilHeat, SensibleHeat and LatentHeat (W/m2)". It does not say that these are instantaneous values at the overpass, how LE relates to mm/h, or what `anchors` contains (pixel, Ts, LAI, NDVI, type, flex). `flex` is only explained in `?calcAnchors`.
- `?ET24h` does not give typical ranges: 0 for dry soil, about 1.0–1.05 × ETr for full-cover well-watered crops.
- `?ET24h` says "C.rad … See Allen (2007)", but no function computes it.
- In `?calcH`, Rn and G have no units, and Z.om has no units in `?calcH` / `?calcAnchors`.

*Why it matters.* Without reference ranges, a non-specialist cannot tell a good map from a bad one.

*Proposed change.*
- One sentence per output layer, with a typical range at a Landsat overpass. For example: "LatentHeat: energy used by ET at the overpass, typically 0–600 W/m2; ET (mm/h) ≈ LE / 680".
- Describe the columns of `anchors` in `?METRIC.EB`.
- In `?ET24h`, add a "How to read the result" paragraph that links to the QC of finding 18.

### E. From the map to irrigation decisions

**24. ETrF is an alfalfa-basis crop coefficient, and nothing explains how to compare it with FAO-56 Kc. — Medium**

*Evidence (`s17`).*
- The example has daily ETr 5.34 and ETo 4.48 mm/day (ratio 1.19). A field with ETrF 0.88 therefore has a grass-basis Kc of 1.05.
- `dailyET(ET = "ETo")` exists, but no vignette or help page shows the conversion. The only mention is in the "flexible" Details.
- Allen et al. (2007, Eq. 54) state that ETrF "is the same as the well-known crop coefficient, Kc, when used with an alfalfa reference basis".

*Why it matters.* Advisors in Chile and Argentina use FAO-56 Kc tables (grass basis) and ETo from the national networks (INIA/Agromet, INTA). Comparing ETrF with those Kc values directly is off by about 20%.

*Proposed change.*
- In `?ET24h` and the vignette: "ETrF = ET24 / ETr24 is a crop coefficient on an alfalfa basis. To compare with FAO-56 Kc (grass basis), use ET24 / ETo24, i.e. ETrF × ETr24 / ETo24 (1.19 in this example)."
- Optionally return `Kc.ETo` as well.

**25. There is no example of summarising ET per field or exporting to GIS, although both work. — Medium**

*Evidence (`s17_outputs.R`).*
- `writeRaster(Energy.Balance, "EB.tif")` keeps the layer names, `terra::writeRaster(ET.24, "ET24.tif")` works, and `sf::st_write(Energy.Balance$anchors, "anchors.shp")` works.
- `terra::extract(c(ET24, ETrF, Kc_grass), vect(fields), fun = mean)` gives the table an advisor needs, for example:

  | Field | ET24 (mm/day) | ETrF | Kc (grass) |
  |---|---|---|---|
  | vines A | 2.76 | 0.52 | 0.62 |
  | olives B | 4.70 | 0.88 | 1.05 |

- None of this appears in the vignettes or help pages.
- `plot(Energy.Balance)` panels have no units (water_LSEBclass.R:11).

*Why it matters.* This last step is what turns a map into an irrigation recommendation.

*Proposed change.* Add a vignette section "From the map to your fields":
- `fields <- sf::st_read("fields.shp")`;
- `terra::extract(..., fun = mean)` → a table of mm/day and Kc;
- `writeRaster()` for QGIS and `st_write()` for the anchors;
- a short note on pixels mixed with the field border (`exact = TRUE`, or a negative buffer).

**26. There is no way to get weekly or seasonal ET between images. — Medium**

*Evidence.*
- Each run gives one day. Allen et al. (2007, Eq. 58-59) compute ET for a period by interpolating ETrF between image dates and multiplying by the daily ETr. This is the paper's method for monthly and seasonal ET maps.
- No function or example does this, and `dailyET()` handles one day at a time.

*Why it matters.* Irrigation is scheduled weekly, from cumulative ET between Landsat dates (8–16 days apart), not from one clear day.

*Proposed change.* Add `ETperiod(ETrF = list(rast1, rast2, …), dates, WeatherStation)` that interpolates ETrF linearly per pixel and sums ETrF × ETr24 per day. Or at least add a vignette example with two dates and a loop over the station's days.

**27. There is no guidance on choosing and checking the weather station. — Medium**

*Evidence.*
- Neither the help pages nor the vignettes say what the station should be like:
  - over well-watered short vegetation (reference conditions);
  - at a reasonable distance (Allen et al. 2007: ETr is "often representative for areas of 50–100 km");
  - with hourly or shorter data;
  - checked for "faulty and biased weather measurements" (Allen et al. 2007, citing ASCE-EWRI 2005).
- The Talca example station is "in an apple orchard" (`?water-extdata`) and is used without comment.
- With an AOI 80 km from the station, the only message is about where the anchors are searched (`s14` p).

*Why it matters.* Daily ETr scales the whole map. Choosing the station is the practitioner's main decision, and it is where their expertise is useful.

*Proposed change.*
- Add a vignette box "Choosing the weather station": reference-condition checklist, maximum distance, plotting the station day with `plot(WeatherStation, hourly = TRUE)` against clear-sky Rso.
- Have `read.WSdata()` print the distance from the station to the image centre.

**28. There is no section on when the results should not be trusted. — Medium**

*Evidence.*
- The hot pixel is always assumed to have ET = 0 (water_sensibleHeatFlux.R:660-661). Allen et al. (2007, Eq. 47) run a soil water balance to add residual evaporation after rain. There is no argument for LE at the hot pixel and no warning, even though the station file has a rain column.
- Early-season and non-growing-season images are handled only by the "flexible" Kc rule, without any explanation for the user.
- Terrain: `plain = TRUE` is the default, with no rule of thumb for when to use a DEM. On the gentle example slopes (median 1°, 95th percentile 2.9°), `getDEM()` + `plain = FALSE` gave mean ET24 4.20 instead of 3.98 mm/day (+5.6%) and a 1-pixel NA border (640 cells instead of 13; `s19`).
- Small AOIs (finding 16), unmasked clouds (finding 15), partial station days (finding 8) and distant stations (finding 27) are not mentioned either.

*Why it matters.* An ET map always looks convincing. A practitioner needs a short list of "red flags" to decide whether to use it.

*Proposed change.*
- Add a "When not to trust the result" checklist to the vignette and to `?METRIC.EB`: rain in the previous ~3–5 days; low vegetation cover; scene or AOI clouds; AOI too small, or ΔT hot–cold < 5–10 K; station far away or with gaps; slopes above about 5° or a large elevation range without a DEM; many pixels with ETrF > 1.1.
- Add an `ET.hot` (or `LE.hot`) argument to `calcH()` / `METRIC.EB()`.

**29. The default roughness method is for short crops, but both examples are vineyards and orchards, and nothing says when to change it. — Medium**

*Evidence.*
- `Zom.method = "short.crops"` (Z.om = 0.018 LAI) is the default. Allen et al. (2007, after Eq. 33): it "tends to underestimate zom for corn … and a different function would be required for vines and trees".
- The examples are Mendoza (INTA, vineyards) and Talca (an apple orchard station), and all vignettes keep the default.
- `?momentumRoughnessLength` Details only gives a definition of Zom. The `Perrier` parameters `fLAI` and `h` come with no typical values, apart from the bare `extraParameters = c(fLAI = 0.5, h = 1.8)` in the `?METRIC.EB` example.
- On the Landsat 8 example, Perrier with those values gives mean ET24 3.49 instead of 3.98 mm/day (−12%; per-pixel 5–95%: −1.05 to +0.02 mm/day; `s27_zom.R`).

*Why it matters.* Vineyards, orchards and olive groves are the main irrigated crops in the regions the package targets. A 12% systematic difference is as large as many irrigation decisions.

*Proposed change.*
- In `?METRIC.EB` / `?momentumRoughnessLength` and the vignette, add a short table of which method to use for which cover (short crops; maize; vines; orchards and olives). Give typical `h` and `fLAI` values with their sources (Santos et al. 2012 and Pôças et al. 2014, already cited).
- Say explicitly that "short.crops" is not meant for trees and vines.

**30. The vignettes contain contradictions and leave the user without a recommendation. — Medium**

*Evidence.*
- **Level-1 claim:** METRIC_simple.Rmd:53 and METRIC_advanced.Rmd:54 list "A raw (Level-1) Landsat 7, 8 or 9 image" as the input. For Landsat 8/9, `METRIC.EB(image.DN = …)` without `image.SR` stops with "image.SR is missing. It can only be computed from image.DN for L7 or MODIS" (`s18`).
- **Two Landsat 8 workflows, no advice:** Landsat8.Rmd runs both, Level-2 ST and Level-1 split window. Mean Ts differs by 6.4 K (311.2 vs 304.9 K), and per-pixel ET24 differs by −0.86 to +1.00 mm/day (5th–95th percentile; r = 0.947; `s22`). The vignette does not say which one to use.
- **Unexplained warning:** the first METRIC.EB run of Landsat8.Rmd produces the warning "Criteria used for cold pixels: …".
- **Unexplained parameters:** the vignette passes `n = 5, Z.om.ws = 0.03` without explanation.
- **Unused aoi:** see finding 14.

*Why it matters.* A new user cannot tell which advice applies to their satellite.

*Proposed change.*
- Fix the input list ("Landsat 7: Level-1; Landsat 8/9: Collection 2 Level-2 SR and ST").
- In Landsat8.Rmd, recommend the Level-2 path for routine use and show how much the two paths differ.
- Explain each non-default argument where it is first used.

**31. Run times are not documented, there is no progress information, and a full scene did not finish in 50 h. — Medium**

*Evidence.* Laptop measurements (`s04`, `s13`, `s16`, `s20`):
- The Landsat 8 vignette, 134 × 184 pixels and two METRIC.EB runs: 42 s.
- 1.07 M pixels (40 × 24 km): METRIC.EB 1.2–1.6 min.
- 4 M pixels (60 × 60 km): load 3 s, cfmask 2 s, **METRIC.EB 280 s**, ET24h 9 s. That is roughly linear, about 70 µs per pixel.
- Full scene (60.4 M pixels): `loadImageSR()` 141 s (it already spills to disk), `cfmask()` 54 s. **My METRIC.EB run had not finished after about 50 hours** (16.8 GB working set, while other reviewers' jobs ran on the same machine), and I stopped it on 18 Oct. G, on the same machine (64 GB RAM) with his steps run one by one, measured about 58 min and a 33 GB peak (G.md finding 5). Either way, a full scene does not fit on a typical 8–16 GB extension-office laptop.
- During METRIC.EB the only signs of progress are the iteration plot and the relaxation lines.

*Why it matters.* Practitioners need to plan a season of images: 20 dates × full scene.

*Proposed change.*
- Add a sentence on expected times and memory to the vignette, and recommend an AOI of 30–60 km around the irrigated area (which also addresses finding 16).
- Print short step messages ("Net radiation… done (12 s)") when `verbose`.

---

## Top 5 priorities

1. **Make the weather input fail loudly.** Match `columns` by name (finding 5). Add unit and range checks with "what to do" messages, plus a check that the overpass falls in a Landsat morning, for km/h, MJ, °F, RH as a fraction, swapped columns, the wrong tz sign, and Kelvin for `surf.temp` (findings 6, 7, 13). Catch the decimal-comma CSV (finding 11), and check the station position and non-finite values in `calcH` (finding 10).
2. **Make daily ETr and ET24 safe and complete.** `dailyET()` should use the overpass date stored in the station object, honour `date`, stop on a missing day, and warn on incomplete or gappy days (findings 8, 9). `ET24h()` should accept the METRIC.EB result and return **ET24 and ETrF** with names and units, and its help should explain ETrF as an alfalfa-basis Kc and how to convert it to FAO-56 Kc (findings 22, 24).
3. **Guard the calibration.** Stop or warn when the hot–cold ΔT is small or the AOI is too small (finding 16). Warn when the MTL reports clouds but nothing was masked (finding 15). Print a QC summary with the share of ETrF > 1.1 (finding 18). Document `n` (finding 17), give the anchor checks a vignette section and a helper (finding 19), and say which roughness method to use for vines and orchards (finding 29).
4. **Remove the hidden global `aoi`.** Set `autoAoi = FALSE` by default, or at least print a message when it is used, and remove the unused `createAoi()` from the Landsat 7 vignettes (finding 14).
5. **Write a practitioner path, in Spanish too.** Fix the README install instructions (the CRAN package is archived; build the vignettes) (finding 1). Make the Landsat 8/9 vignette the entry point (finding 2). It should go from download (guides shipped with the package, with AOI-cropped downloads and a hand-off to `loadImageSR()`; findings 3, 4) to a **per-field table of ET (mm/day) and Kc plus GeoTIFF export** (finding 25), with "Choosing the station" and "When not to trust the result" boxes (findings 27, 28). Seasonal ET (finding 26) can follow next.
