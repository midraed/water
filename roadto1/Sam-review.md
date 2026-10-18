# Scientific review of water 0.8.5 (master 70027b3) against Allen, Tasumi & Trezza (2007)

Reviewer: Sam (land-surface energy balance, METRIC/SEBAL, ASCE-EWRI reference ET)

**Scope and method.** I checked every METRIC step in `R/` against Allen et al. (2007, J. Irrig. Drain. Eng. 133(4):380–394; equation and page numbers below refer to it) and the hourly reference ET against ASCE-EWRI (2005). I used the full report text, which I downloaded. I verified the numbers by running the installed master (`lib_new`) on the shipped examples:

- the Landsat 8 C2 L2 crop with INTA.csv (Mendoza);
- the L7_Talca and DEM_Talca data with apples.csv (Talca);
- the full L8 QA_PIXEL scene.

I also wrote an independent strict-ASCE hourly implementation, and used patched copies of single functions inside my own R session. No repository file was changed. The scripts (`s01`–`s20`) are in `<scratch>`.

The positive part first: the core physics is a faithful transcription of the paper. The checked items are listed at the end. Most of the problems are in what surrounds that core:
- anchor selection,
- the reference-ET bookkeeping,
- a few defaults that differ from the paper,
- missing diagnostics.

These matter because METRIC moves all the errors into the calibration of H.

---

## Findings

### 1. Anchor-pixel criteria are site-specific, hard-coded, and decide the ET map — **high**

**Evidence.**
- `calcAnchors()` (`R/water_sensibleHeatFlux.R:321-336`, also `:156-169`, `:233-246`) fixes these windows:
  - cold: albedo 0.18–0.25, LAI 3–6, z0m 0.03–0.08;
  - hot: albedo 0.13–0.15, NDVI 0.10–0.28, z0m ≤ 0.005.
- The windows come from "CITRA and MCB (personal communication)". They are not in Allen et al. (2007). The paper defines the cold pixel as full-cover, well-watered agricultural vegetation with LAI > 4 (p. 387–388) and the hot pixel as dry bare agricultural soil (p. 387).
- Because z0m = 0.018·LAI (Eq. 33), the z0m upper bound of 0.08 caps the cold LAI at **4.44**. The effective cold window is therefore LAI 3–4.44, which excludes the fullest canopies the paper asks for. In the L8 example this excludes 272 of the 1849 pixels with LAI ≥ 3.

L8 example (vignette call, `s01`/`s02`/`s14`):
- The windows do not match this scene:
  - median albedo of full-cover pixels is 0.153–0.156 (Tasumi/Olmedo), and only 20–22% of them fall inside 0.18–0.25;
  - median albedo of bare pixels is 0.194, and only 6.5–6.8% fall inside the hot window.
- The search had to relax the criteria (flex = 0.11). The cold anchors (Ts 308.2–310.1 K) sit in the middle of the vegetated Ts distribution: **23% of all pixels are colder than their mean**.
- Result: **10.9% of pixels have ETrF > 1.1** (3.4% > 1.2, maximum 1.57). The paper expects ETrF from 0 to about 1.0, "perhaps up to 1.1" (p. 389).
- The coldest 1% of full-cover pixels (albedo 0.11–0.17) get ETrF 1.33–1.52.
- At the dry end, 1044 pixels are hotter than the mean hot anchor, and **4.8% of pixels have H > Rn − G**, so LE is clipped to 0.
- Sensitivity: using the 5 coldest full-cover pixels (any albedo) as cold anchors, with the hot anchors unchanged, lowers mean ETrF from 0.745 to 0.515. Mean ET24 drops from 3.98 to **2.75 mm/d (−31%)**, and ETrF > 1.1 disappears.

Other problems in the selection logic:
- The "flexible" relaxation multiplies the hot Ts threshold, in Kelvin, by (1 − flex) (`:323-324`, `:349`). Each 1% step lowers it by about 3 K.
- "best" and "flexible" take the single coldest and hottest 30 m pixel (`:382`, `:413`). There is no homogeneity test, although TIRS is 100 m, and no check of ST_QA or distance to clouds.

**Why it matters.** METRIC moves all its biases into the dT(Ts) line through the two anchors. The paper warns that METRIC "is highly reliant on the operator's ability to select appropriate cold and hot pixels. It is currently not a hands-off model" (p. 393). A fixed albedo window tuned for one valley can move scene ET by tens of percent without any warning.

**Proposed change.**
- Make all windows arguments, and default to scene-relative (percentile) rules in the spirit of Allen et al. (2013, JAWRA 49:563–576, doi:10.1111/jawr.12056).
- Set cold LAI ≥ 4 (or a high NDVI percentile) and drop the z0m upper bound for cold pixels.
- Relax temperature thresholds in Kelvin, not percent.
- Average a few homogeneous candidates instead of taking the extreme pixel. Optionally screen with ST_QA and ST_CDIST.
- Report diagnostics in `print(waterLSEB)`: fraction of ETrF > 1.1, fraction of LE < 0, and the percentile of the anchors in the Ts distribution of vegetated and bare pixels.

### 2. Daily (and night-time hourly) reference ET is not ASCE-EWRI 2005: +9% on a clear day — **high**

**Evidence.** In `hourlyET()` (`R/water_ET.R:109-123`):
- **(a) Night-time cloudiness.** At night Rso ≈ 0, so Rs/Rso is set to 0 (`:113`) and fcd = 0.05 (`:114`). This means "overcast" every night, and night-time net longwave is about 0. ASCE Eq. 46 instead carries forward fcd from the last period with sun angle β > 0.3 rad. It applies that value to dusk, night and dawn (pp. 35–36).
- **(b) Rs/Rso clamp.** Rs/Rso ≤ 0.3 is set to 0 instead of 0.3 (ASCE Eq. 45 limits it to 0.3–1.0).
- **(c) Daytime definition.** Daytime for Cd is decided by Rs > 0 (`:123`). ASCE defines daytime as hourly Rn > 0 (p. 5 and p. 44; Table 1).
- **(d) Ra integration.** Ra uses clamped ω1, ω2 in the first term but unclamped values in sin ω2 − sin ω1 (`:110-111`), against ASCE Eq. 56.

Reproduction (`s04`, `s15`), strict ASCE implementation on the package's own hourly table:

| Day | Reference | dailyET() | Strict ASCE | Bias |
|---|---|---|---|---|
| Talca, 15 Feb 2013 (clear) | ETr | 7.09 mm | 6.50 mm | **+9.1%** |
| Talca, 15 Feb 2013 (clear) | ETo | — | — | +7.0% |
| Mendoza, 9 Feb 2016 (cloudy afternoon) | ETr | — | — | −1.2% |

- On the Talca day the hour 19–20 h alone is 0.26 vs 0.12 mm.
- For Talca the bias ranges from +4.1% (negative hours set to 0) to +11.6% (night fcd = 1, as in the refet package), depending on convention.
- Minor: dailyET() sums samples taken at the hour mark rather than hourly means of the 15-min data. This changes the Talca ETr24 by −1.2% (`s16`).

**Why it matters.** ET24 = ETrF × ETr24 (Eq. 55), so the bias passes 1:1 into every pixel of the daily map. It is largest on clear days with dry air, which is exactly when Landsat images are usable.

**Proposed change.**
- Implement Eq. 46: one fcd>0.3 per day, applied from β < 0.3 rad until β > 0.3 rad again.
- Limit Rs/Rso to [0.3, 1].
- Switch Cd and G on Rn > 0.
- Use clamped ω in both Ra terms.
- Add a regression test against ASCE Appendix C or the refet package.

### 3. The weather data at the overpass are not checked; the shipped L8 example calibrates the cold pixel with a cloud-affected Rs — **medium**

**Evidence** (`s04`, `s12`, `s13`).
- At the 11:27 overpass, the INTA station's interpolated Rs = 587 W/m², against ASCE Rso = 860 W/m² (**Rs/Rso = 0.68**). The package's own clear-sky Rs↓ is 832 W/m².
- QA_PIXEL shows no cloud or shadow flag within 5 km of the station, but 10.8% flagged within 20 km. Hourly means were likely depressed by passing clouds or a low-reading pyranometer: the morning ratio stays around 0.74 every hour.
- Matching 587 W/m² with Eq. 3–4 would require Kt ≈ 0.17. The paper's "extremely turbid" bound is Kt = 0.5 (p. 382).
- Consequence for the cold pixel (Eq. 48):

| Rs used for ETr_inst | ETr_inst | H_cold | Bowen ratio |
|---|---|---|---|
| Station, 587 W/m² | 0.499 mm/h | 181 W/m² | 0.52 |
| Clear-sky Rso, 860 W/m² | 0.679 mm/h | 54 W/m² | 0.11 |

  A Bowen ratio of 0.52 is implausible for a well-watered full-cover pixel; 0.11 is typical.
- In addition, `read.WSdata()` interpolates `at.sat` treating hourly records as instantaneous at their timestamps (`R/water_WeatherStation.R:222-231`). `hourlyET()`, on the other hand, treats the same records as means of the previous hour. As period means centred at :30, Rs at the overpass would be 638 W/m² (ETr_inst 0.532 mm/h).
- Kt is fixed at 1 (`R/water_netRadiation.G.R:321`, `:688`).

**Why it matters.**
- Rn uses computed clear-sky Rs↓ (allowed by the paper, p. 382), but ETr uses measured Rs, so the two are inconsistent on such days.
- Allen et al. (2007) require "high-quality hourly (or shorter) weather data" (p. 381, p. 393). ASCE-EWRI (2005, App. D, with Allen 1996) prescribes Rs-versus-Rso integrity checks.

**Proposed change.**
- Warn when Rs/Rso at the overpass of a clear image is below about 0.85.
- Offer ETr_inst computed with the modelled clear-sky Rs↓ for the calibration.
- Interpolate at period midpoints when records are period means.
- Expose Kt.
- Add a station-QC plot (Rs vs Rso envelope).

### 4. Collection 2 SR outside the valid range is kept, and the cold NDVI threshold uses the image maximum — **medium**

**Evidence.**
- `loadImageSR()` masks only the fill value 0 (`R/water_RSfunctions.R:148`). USGS gives the valid SR range as 7273–43636 DN, i.e. reflectance 0–1 (USGS FAQ "How do I use a scale factor with Landsat Level-2 science products").
- `calcAnchors()` requires NDVI ≥ max(NDVI) − 0.15 (`:325`, `:241`). The maximum is taken after discarding only NDVI > 1.
- In the L8 example, 109 pixels (0.44%) have negative SR, mostly in blue, plus 7 in red. They raise max NDVI from 0.923 (valid pixels) to 0.987. The cold threshold rises from NDVI 0.77 to 0.84.
- Masking them (`s09`):
  - the relaxation disappears (flex 0.11 → 0);
  - the cold anchors change (NDVI 0.84 → 0.78–0.81);
  - mean ETrF goes from 0.745 to 0.777, and mean ET24 from **3.98 to 4.15 mm/d (+4.3%)**;
  - median pixel ETrF changes by +0.036.

**Why it matters.** Under half a percent of artefact pixels moves the whole calibration, and the switch into "flexible" relaxation also decides whether the cold-pixel ETrF is replaced (finding 14).

**Proposed change.**
- Set SR outside 0–1 to NA in `loadImageSR()`, keeping the >1 tolerance for snow if wanted.
- Use a high percentile (e.g. the 99th) of NDVI instead of the maximum.

### 5. Thermal path: τNB is not the paper's value, and the "SC" method is not Jiménez-Muñoz et al. (2014) — **medium**

**Evidence.**

(a) τNB.
- `surfaceTemperature()` uses τNB = 0.75 + 2×10⁻⁵ z (`R/water_netRadiation.G.R:572`). That is the ASCE/FAO broadband clear-sky shortwave transmissivity (ASCE Eq. 47).
- For Eq. 21 the paper gives τNB = 0.866, together with Rp = 0.91 and Rsky = 1.32 W m⁻² sr⁻¹ µm⁻¹ (p. 384). Rp and Rsky are implemented correctly: Rsky ≈ 1.33 at 30 °C.
- In the L7 Talca example (`s06`), τNB = 0.754 makes Ts warmer by **+9.6 to +11.0 K** (median +10.1 K):
  - the cold anchors are 12.4 K above air temperature, against +2.3 K with 0.866;
  - Rn changes by −13 to −18 W/m²;
  - bare-soil G (Eq. 27b) changes by up to +18 W/m²;
  - the calibration absorbs most of this: pixel ETrF changes by −0.04 to +0.02 (5–95%).

(b) Landsat 8/9 "SC" method.
- It computes the Jiménez-Muñoz atmospheric functions from ea in kPa, where the method needs water vapour in g/cm² (`:595-598`, the code comment admits it). It then never uses them (`AF.j`).
- Instead it uses the Landsat 7, Idaho defaults with the wrong τNB.
- Against the USGS ST_B10 for the same scene (`s07`):

| Method | Ts − ST_B10 (median) |
|---|---|
| Package SC | +3.2 K |
| Package SW (band 11, fixed w = 1.2 g/cm², `:608`) | −6.3 K |
| Proper JM2014 SC (w = 2.6 g/cm²) | −5.9 K |

- USGS advises against relying on TIRS band 11 for split-window retrievals (stray light).

**Why it matters.** Ts is an output layer of METRIC.EB(). It drives RL↑, RL↓ (each pixel's Ts is used as Ta), G and the cold-pixel plausibility. A 10 K bias is physically wrong even where ET survives.

**Proposed change.**
- Use τNB = 0.866 or let the user set τNB, Rp and Rsky.
- For L8/L9, take the per-pixel ST_ATRAN, ST_URAD and ST_DRAD bands that ship with C2 L2 (these are τNB, Rp and Rsky), or simply default to ST_B10.
- Fix the SC method (water vapour from Garrison & Adler, W/10) or remove it.
- Document the band-11 caveat for SW.

### 6. The hot pixel always has ET = 0; Eq. 47 cannot be used — **medium**

**Evidence.**
- `calcH()` sets dT_hot = (Rn − G)·r_ah/(ρ cp) (`R/water_sensibleHeatFlux.R:660-661`, `:763`). There is no LE_hot term and no argument for it.
- The paper runs an FAO-56 bare-soil water balance to confirm ET_hot = 0 "or to supply a nonzero value for ET if there is residual evaporation", giving H_hot = (Rn − G)_hot − LE_hot (Eq. 47, p. 387).
- The station "rain" column is read but never used.

**Why it matters.** After rain or irrigation, assuming ET_hot = 0 steepens dT(Ts) and underestimates ET for every low- and mid-ET pixel, by roughly the ignored ETrF_hot.

**Proposed change.**
- Add an `ETrF.hot` argument (default 0), symmetric to `ETp.coef`.
- Warn when the station reports rain in the days before the image. Optionally provide an FAO-56 Ke helper.

### 7. Low-wind handling: u200 is replaced by 4 m/s, and the stability iteration oscillates — **medium**

**Evidence.**
- `calcH()` replaces any u200 below 1 m/s by 4 m/s (`:593-597`). This is discontinuous: u2 = 0.50 m/s gives u200 = 0.97, replaced by 4; u2 = 0.52 m/s gives 1.01, kept. The paper prescribes nothing.
- **Both L7 vignettes trigger it.** The Talca wind at the overpass is 0.31 m/s, so u200 = 0.59 m/s is used as 4.0 (×6.8).
- With a floor of 1 m/s instead (`s18`, `s19`):
  - mean ETrF changes from 0.819 to 0.842 (pixels −0.06 to +0.10);
  - the iteration does not converge in 20 iterations. r_ah,cold oscillates 141 → 2.6 → 98 → 6.6 → 65 → … → 25.3 → 26.7.
- The L8 example converges, but with a large damped oscillation (r_ah,hot 77.5 → 5.2 → 21.3 → 12.4 → … → 14.5).
- Convergence is tested only at the anchors (`:724-737`).

**Why it matters.** Calm mornings are common in irrigated valleys. Under free convection the Monin–Obukhov iteration (Eqs. 38–45) needs damping. An arbitrary 4 m/s hides a solver problem and changes r_ah everywhere.

**Proposed change.**
- Under-relax the iteration (e.g. average successive r_ah or H).
- Test convergence on a high percentile of |Δr_ah| over all pixels.
- Use a continuous, documented wind floor.
- Report the u200 actually used in `$methods`.

### 8. C_rad (Eq. 56) is never computed for sloping pixels — **medium (mountain mode)**

**Evidence.**
- `ET24h()` takes a scalar `C.rad = 1` (`R/water_ET.R:157`, `:166`). No function computes the per-pixel ratio [Rso(inst)_hor/Rso(inst)_pixel]·[Rso(24)_pixel/Rso(24)_hor] (Eqs. 56–57, p. 390), even though `METRICtopo()`, `solarAngles()` and `incSWradiation()` already contain every ingredient.
- Integrating Eqs. 3–4 over the day on DEM_Talca (`s08`): for the 11.7% of pixels with slope > 5°, C_rad = **0.88–1.29** (5–95%), and 0.81–1.61 (1–99%).

**Why it matters.** With `plain = FALSE`, daily ET on slopes is off by that factor. The instantaneous ETrF on a sun-facing morning slope is not representative of its 24-h energy.

**Proposed change.** Add a `Crad()` function (Rso_i from Eq. 3 with cos θrel clamped ≥ 0) and use it by default in `ET24h()` when a surface model is given.

### 9. dailyET() picks the day by UTC date and does not check that 24 hours are present — **medium**

**Evidence.**
- The day is chosen by `trunc(datetime, "days") == as.Date(sat.date)`, where `sat.date` is the UTC DATE_ACQUIRED (`R/water_ET.R:238-240`).
- For stations at UTC+11 to +13 (New Zealand, south-east Australia in summer), the local morning overpass falls on the previous UTC date.
- Synthetic test in Pacific/Auckland (`s10`): with the overpass at 10:30 NZDT, dailyET() returned the previous local day (4.85 mm) instead of the image day (14.66 mm).
- The `for(i in 1:24)` loop (`:245-249`) with `sum(na.rm = TRUE)` silently underestimates ETr24 when hours are missing.

**Why it matters.** ETr24 multiplies every pixel (Eq. 55).

**Proposed change.**
- Select the local date of `at.sat$datetime`.
- Require 24 records (or warn) before summing.

### 10. MODIS albedo uses Liang's direct near-infrared formula instead of the shortwave one — **medium (MODIS path)**

**Evidence.**
- `albedo()` uses c(0.037, 0.479, −0.068, 0.0976, 0.266, 0.0757, 0.107) for MODIS (`R/water_netRadiation.G.R:369`; the code comment says "direct-NIR"). In Liang (2001, Eq. 15) this is exactly α_direct-NIR.
- The shortwave formula is α_short = 0.160α1 + 0.291α2 + 0.243α3 + 0.116α4 + 0.112α5 + 0.081α7 − 0.0015. Allen et al. (2007, Table 2) also give Tasumi weights for MODIS.
- `coeff` is ignored for MODIS.
- Illustration with typical green-vegetation reflectances: α ≈ 0.32 (package) vs 0.18 (Liang shortwave) or 0.16 (Tasumi). That is about −110 to −130 W/m² of Rn at Rs↓ = 800 W/m².

**Proposed change.** Use Tasumi Table 2 (default) or Liang shortwave for MODIS, and honour `coeff`.

### 11. The z0m minimum of 0.005 m (text after Eq. 33) is missing; LAI = 0 gives u* = 0.1 regardless of wind — **low**

**Evidence.**
- `momentumRoughnessLength()` returns 0.018·LAI with no floor (`:44`). Perrier also gives 0 at LAI = 0 (`:52`).
- The paper states "a minimum value for zom = 0.005 m represents roughness typical of bare agricultural soils" (p. 386).
- With z0m = 0, u* = 0 is replaced by 0.1 (`:606`, `:707`). That u* no longer depends on wind or stability.
- In the L8 example:
  - the hot anchors have z0m = 0.0005–0.002 m (neutral r_ah 74–83 s/m, against 68 with the floor);
  - adding the floor (`s03`) changes mean ETrF by −0.012, and by −0.019 for LAI 0.28–1;
  - for the few pixels with LAI = 0 (0.02–1.1% depending on the LAI method), mean H is 23–37% higher than for pixels with 0 < LAI < 0.3 at the same Ts (`s17`).

**Proposed change.**
- Apply `pmax(z0m, 0.005)` for all methods.
- Note in the docs that Eq. 33 is for short crops: the paper says "a different function would be required for vines and trees". A single image-wide crop height h in Perrier is equally limited.

### 12. Other mountain-mode inconsistencies — **low**

**Evidence** (`s08`).
- **(a) Unnormalised C2 SR on slopes.** LaSRC SR is normalised to a horizontal surface, but `albedo()` uses it directly while Rs↓ uses θrel (Eq. 3). Eq. 10 normalises by cos θrel. Albedo is therefore scaled by cos θrel/cos θhor:
  - 0.93–1.04 (5–95%), with extremes of 0.25–1.29;
  - for α = 0.20, the absorbed-shortwave error is −6.7 to +10 W/m² (5–95%) and −59 to +39 W/m² at the extremes.
- **(b) Anchors on Ts.** Anchors are selected on Ts, not on Ts_datum. With the lapse correction, 2 of the 5 hot anchors change.
- **(c) Eq. 35 on flat pixels.** As coded (`:55`), Eq. 35 reduces z0m by up to 25% on slopes < 5° (88% of the Talca pixels).
- **(d) Shaded slopes.** cos θrel < 0 would give negative Rs↓ (`R/water_netRadiation.G.R:324`). There is none in Talca, but it is possible in steep terrain.
- **(e) Fixed lapse rate.** The lapse rate is fixed at 6.49 K/km (`:585`). The paper uses a "customized lapse rate".

**Proposed change.**
- Apply a slope correction to SR when `plain = FALSE`.
- Select anchors on Ts_datum.
- Use max(1, Eq. 35 factor). This is my reading of the intent; please check against the METRIC manual.
- Clamp cos θrel ≥ 0.
- Expose the lapse rate.

### 13. Energy-balance closure and outputs — **low**

**Evidence.**
- METRIC.EB() sets LE < 0 to 0 but leaves H unchanged (`R/water_METRIC.R:481`). In the L8 example 379 pixels (1.5%) have H > Rn, up to 1031 W/m², and 4.8% have H > Rn − G. The returned layers therefore do not close.
- ETrF and ET_inst, the products the paper uses for interpolation between images (Eqs. 54, 58–59), are not returned. There is no helper for seasonal ETrF × ETr24.

**Proposed change.**
- Return ETrF and ET_inst.
- Either cap H at Rn − G where LE is clipped, or flag those pixels.
- Consider an Eq. 58 helper.

### 14. "Flexible" cold-pixel ETrF update: the trigger and continuity — **low**

Context: the maintainer's ETo/ETr conversion is accepted.

**Evidence.**
- The Reyes-González Kc replaces ETp.coef only if `flex.cold > 0 & NDVIcold < 0.8` (`R/water_METRIC.R:456-472`). So it depends on whether the *search* was relaxed, not on the state of the vegetation. With flex = 0, a cold pixel with NDVI 0.75 keeps 1.05.
- At the 0.8 boundary the coefficient jumps from 1.021 to 1.05 (with ETo/ETr = 0.874 here). The jump is larger where ETo/ETr is lower.
- NDVI is TOA with image.DN and SR without it. The paper's alternative for partial cover is ETr,cold/ETr ≈ 1.25–1.30·NDVI on TOA (pp. 387–388).

**Proposed change.**
- Decide on the cold pixels' NDVI or LAI alone.
- Use min(1.05, f(NDVI)) so the coefficient is continuous.
- State which reflectance the NDVI must come from.

### 15. Water, snow and desert special cases are not implemented — **low**

**Evidence.**
- For NDVI ≤ 0 the paper sets εNB = ε0 = 0.985 (p. 384) and G/Rn = 0.5 for water and snow (p. 385).
- The package gives water ε0 = 0.95 and Tasumi G. The 56 water pixels of the L8 example get G/Rn = 0.285.
- The paper's desert G limitation and the reduced dT slope above Ts_hot (pp. 385, 389) are absent.

**Proposed change.** Add the NDVI ≤ 0 branch (water or snow by albedo and Ts < 277 K, as in the paper), or mask those pixels.

### 16. Minor constant and formula deviations — **low**

**Evidence.**
- **(a) Eq. 4.** `SWtrasmisivity()` uses −0.00149 instead of the paper's −0.00146 (`R/water_netRadiation.G.R:500`): τsw −0.18%, Rs↓ about −1.5 W/m².
- **(b) Eq. 37 air density.**
  - `rho` divides by Ts instead of Ts − dT, and builds pressure from Ts − dT instead of 293 K (`R/water_sensibleHeatFlux.R:672`).
  - `air.density` (`:587`, used for dT_cold, dT_hot and L) also uses Ts.
  - Effect: ρ is −1.2% at cold and −2.2% at hot pixels, and H_hot = 1.0046·(Rn − G) at 927 m. Most of this is self-cancelled by the calibration.
- **(c) Constants.** 273.16 in Eq. 27b (`:761`) and 0.002361 in λ (`:611`) are negligible.
- **(d) Sunrise and sunset hours.** `hourlyET` sets Ra ≈ 0 whenever β at the period midpoint is ≤ 0, even if the sun is up for part of the hour (`R/water_ET.R:109`).

**Proposed change.** Align with the equations as written. Each item is cosmetic on its own, but together they make the code harder to audit.

### 17. "Olmedo" albedo weights are undocumented — **low**

**Evidence.**
- The weights c(0.246, 0.146, 0.191, 0.304, 0.105, 0.008) (`:366`) are described only by a code comment ("SMARTS … Direct Normal Irradiance"). Broadband albedo weights should come from global (direct plus diffuse) at-surface irradiance (Eq. 15, Tasumi et al. 2008).
- The SWIR2 weight is 0.008, against Tasumi's 0.036.
- The effect in the example is small: median −0.002 compared with Tasumi.

**Proposed change.** Publish or derive the weights reproducibly with global horizontal irradiance. Until then, document Tasumi on OLI as the reference choice.

### 18. Landsat 7 TOA reflectance uses hard-coded low-gain calibration — **low**

**Evidence.**
- `calcTOAr()` hard-codes the ETM+ low-gain G and B values (`R/water_RSfunctions.R:256-257`) and ignores the MTL RADIANCE_MULT_BAND_x and RADIANCE_ADD_BAND_x keys. `calcRadiance()` already reads those keys.
- The example scene is low gain in every band, so the problem is invisible here. High-gain scenes would get wrong red and NIR reflectance, and therefore wrong SAVI, LAI and anchors.

**Proposed change.** Read the gains from the MTL (or use REFLECTANCE_MULT and REFLECTANCE_ADD from Collection 2 Level-1).

### 19. LAI from SR (comment on a known decision) — **low**

**Evidence** (`s11`).
- The METRIC SAVI–LAI relations (Eq. 18–19 and the 2010 cubic) were fitted on TOA reflectance.
- On the same scene, SR raises LAI (median +0.21, 95th percentile +1.23). It triples the pixels with LAI ≥ 3 (608 → 1849) and gives 272 pixels with LAI > 4.44 instead of 32. This shifts anchor eligibility (finding 1) and also changes G/Rn, ε0 and z0m.
- The about 1% ET effect the maintainer measured is plausible after calibration. The anchor criteria are the more sensitive part.

**Proposed change.** Tie the anchor thresholds to the reflectance type, or offer an SR-based LAI relation for crops.

### 20. calcH(method = "lm") uses only the last anchor pair — **low**

**Evidence.**
- It iterates each pair and returns the a and b of the last one (`R/water_sensibleHeatFlux.R:740-842`). r_ah carries over between pairs. This is documented, but it is not a regression, and it discards n − 1 pairs.

**Proposed change.** Fit dT against Ts_datum across all converged pairs, or remove the option.

---

## Checked and consistent with the paper

- Eq. 2 (Rn).
- Eq. 3 and Eq. 9 (d²).
- Eq. 5–6 (P, W).
- Eq. 7 incidence angle, including the aspect convention after `aspect − π`. Sanity check (`s20`): with the sun at azimuth 64.6°, E/NE slopes get 1.12–1.16×, W slopes 0.75× and S-facing slopes 0.83× the horizontal illumination.
- Hour angle and equation of time.
- Eq. 10 (L8 constants, cos θrel); Eq. 11–14 (all Table 1 constants); Eq. 15 (Tasumi weights).
- Eq. 17 and Eq. 22 (emissivities); Eq. 18–19 (LAI with the 0.687 and 0.1 limits); METRIC 2010 cubic.
- Eq. 21 Rp and Rsky.
- Eq. 24–25. Per-pixel Ts as Ta is allowed (p. 384); Ts_cold as Ta would change Rn by −22 to +54 W/m² here.
- Eq. 26–27.
- Eq. 28–31; Eq. 32 (u200).
- Eq. 34b; Eq. 36 (ω).
- Eq. 38–45 (Paulson and Webb functions; stable case with z = 2 m).
- Eq. 46, 48–50. The code fits the intercept correctly and does not reproduce the paper's misprinted Eq. 51.
- Eq. 52–55; λ (Eq. 53); ETr24 from summed hourly ETr.
- ASCE hourly Cn/Cd/G constants, u2 conversion and Δ.
- Solar-time correction with long.z including daylight saving time (correct for clock-time records).
- Z.om.ws = 0.015 m (consistent with ASCE's 4.87/ln(67.8z − 5.42), whose 1/67.8 = 0.0148 m).
- Cold-pixel ET = 1.05·ETr (p. 387).
- Using ST_B10 instead of METRIC Ts is justified, because dT is indexed to Ts and needs no absolutely accurate surface temperature (p. 381). Its uncertainty layers (ST_QA, ST_CDIST) should be used to screen anchors.

## Top 5 priorities

1. **Anchor selection (finding 1, with 4 and 19).** Use scene-relative, configurable criteria (Allen et al. 2013): cold LAI ≥ 4 with no z0m cap, a percentile NDVI threshold, relaxation in Kelvin, and homogeneous candidates. Print calibration diagnostics: ETrF > 1.1, LE < 0, anchor percentiles. In the example, this one choice moves mean ET24 by up to 31%.
2. **ASCE-compliant hourly and daily reference ET (findings 2 and 9).** Night-time fcd per Eq. 46, Rn > 0 daytime, Rs/Rso limits, clamped Ra, local-date selection and a completeness check. These remove a +9% clear-day bias that passes 1:1 into ET24.
3. **Mask Collection 2 SR outside 7273–43636 DN (finding 4).** A two-line change that moved mean ET24 by 4.3% in the shipped example.
4. **Weather-data QC at the overpass (finding 3).** Check Rs/Rso, offer clear-sky ETr_inst for calibration, interpolate period means at midpoints, expose Kt. The shipped L8 example has a cold pixel with H = 181 W/m² (Bowen 0.52) from a cloud-affected Rs.
5. **A robust stability iteration at low wind (finding 7).** Damping, convergence over all pixels, and a continuous documented u200 floor instead of 0.59 → 4 m/s, which both L7 vignettes trigger.

Next in line: τNB = 0.866 or the per-pixel ST_ATRAN/URAD/DRAD (finding 5); an ETrF_hot option per Eq. 47 (finding 6); C_rad (finding 8); the z0m floor (finding 11); MODIS albedo (finding 10).
