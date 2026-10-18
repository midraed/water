# Sam: discussion round

## 1. Science-relevant findings of Gabe (G) and Rosa (R)

- **G3, calcH oscillation: agree, high.** Same mechanism as my finding 7. With u200 = 1 m/s I see r_ah,cold go 141 → 2.6 → 98 → … → 26.7 without converging. Damping is the real fix: average successive r_ah (or iterate on log r_ah), test convergence on a and b and on a high percentile of |Δr_ah| over all pixels, and return a `converged` flag.
- **G13, u200 floor: agree, but do not rely on the floor alone.**
  - Recommended: a continuous clamp, `u200 <- max(u200, 2)` m/s. With Z.om.ws = 0.015 that equals u2 ≈ 1 m/s at the station.
  - 1 m/s is too low: the undamped loop oscillates there (finding 7).
  - 4 m/s is arbitrary and makes results jump.
  - The paper gives no value, so this is a judgement call. Document it and report the u200 actually used in `$methods`.
  - With damping in place, the floor only guards against nonsense inputs.
- **G15, SC method: agree.** This is my finding 5b. Preferred: drop "SC" for L8/L9 and default to the Level-2 ST_B10. If a DN path is kept, use METRIC Eq. 20–21 with the per-pixel ST_ATRAN, ST_URAD and ST_DRAD (or τNB = 0.866). Do not keep a hybrid labelled Jiménez-Muñoz.
- **G12, Level-1 fill: agree, raise to high for DN users.** Ts = −53.6 K pixels enter the cold-anchor search (min Ts) and the Ts quantiles. It is the same class of bug as my finding 4 (out-of-range SR), and the fix is equally cheap. Mask DN = 0 and SR outside 7273–43636 DN.
- **G4 and R8/R9, dailyET day and gaps: agree.** This is my finding 9. Rank it with my ASCE fix (finding 2): both scale ET24 1:1.
- **G9, anchor validation: agree.** Add "mean Ts_hot − Ts_cold > ~5 K" to the checks. That also covers R16.
- **G14, bands by position: agree, high impact.** LAI 1.25 → 0.001 changes ε, G, z0m and the anchors silently.
- **R13, °C surf.temp: agree.** It is also an absolute-Ts error: Rn and G use σTs⁴.
- **R16, small AOI: agree, high.** METRIC needs both extremes inside the calibration area (paper p. 393).
  - Recommend: warn when ΔTs between anchors is < 5 K, or when there are fewer than about 20 cold or hot candidates.
  - Document "calibrate on tens of km, then crop the results".
- **R17, `n`: agree.** Prefer a default n of about 3–5 homogeneous anchors, but averaged over candidates within a temperature band rather than "the 1st…5th coldest at 500 m". This ties into my finding 1.
- **R18, ETrF QC: agree.** It is the same evidence as my finding 1 (16.6% > 1.05, 10.9% > 1.1).
  - The QC is necessary but not sufficient: the cause is the anchor windows.
  - Wording for the histogram: values of 1.1–1.2 can be wet surfaces; values above 1.2 over vegetation mean "calibration problem", not "vineyard uses more than alfalfa".
- **R24, ETrF to Kc: agree.** Kc(grass) = ETrF × ETr24/ETo24 for the daily value. Note: use daily ratios, not the instantaneous ETo/ETr at the overpass.
- **R26, seasonal ET: agree.** This is Eq. 58–59 (my finding 13). Interpolate ETrF, never ET.
- **R29, roughness for vines and orchards: agree, medium.**
  - Caution on Perrier: one scene-wide h applies tree height to bare soil and alfalfa too, and gives z0m = 0 at LAI = 0 (my finding 11).
  - Recommend a land-cover mask (or an h raster), the 0.005 m floor, and the documentation table Rosa asks for.
  - The −12% she measured is real sensitivity, not a correction to apply blindly.
- **R12, cryptic u200 warning: agree.** It is resolved by the G13 and G3 fixes plus a plain-language message.
- **R27, station choice: agree.** Add my overpass Rs/Rso check (finding 3) to that guidance.

## 2. Conflicts between perspectives and proposed resolutions

- **Scene-percentile anchor criteria (Sam) vs simple defaults (Rosa).**
  - No real conflict: percentiles are a simpler default for the user than the current hidden albedo windows.
  - Proposal:
    - one default (cold = NDVI ≥ 95th percentile and LAI ≥ 4, coldest 20% of those; hot = NDVI ≤ 10th percentile, hottest 20%; average of homogeneous candidates);
    - a single `anchors.control = list(...)` for experts (this also matches R21);
    - an always-printed QC line.
  - The windows stay available as `anchors.method = "citra"`.
- **Failing loudly (Rosa, Gabe) vs a "flexible" method that never fails.**
  - Flexible relaxation hides exactly the cases where METRIC is invalid (R16, my findings 1 and 4).
  - Keep it, but make it warn with the reason and record it in `$methods` (flex, criteria used).
  - Stop when ΔTs < 5 K.
  - Make the ETrF,cold change depend on vegetation state, not on flex (my finding 14).
- **Performance (Gabe 5: iterate on anchors, then one pass) vs damping and an all-pixel convergence test (me).**
  - Compatible: iterate the anchors to convergence, then iterate each pixel's stability locally with damping. This is cheap per pixel (vectorised, about 10 passes).
- **Exposing physical knobs (my Kt, τNB, lapse rate, ETrF.hot) vs fewer arguments (Rosa R21).**
  - Group them in `physics.control` with paper defaults: Kt = 1, τNB = 0.866, lapse 6.49 K/km, ETrF.hot = 0.
  - Only ETrF.hot gets a vignette mention ("after rain").
- **Returning ETrF (Rosa) vs closure (me).**
  - Return ETrF and ET_inst unclipped for statistics, plus clipped maps for display.
  - Flag pixels where H > Rn − G.

## 3. Updated team top-5 (science-weighted)

1. **Anchor selection and calibration guards.** Percentile defaults (R17, R19), ΔTs/AOI checks (R16, G9), the QC summary (R18), masking of invalid input (G12, my finding 4). This is the largest effect: ±31% in mean ET24, about 40% for farm-sized AOIs.
2. **Reference ET correctness.** ASCE Eq. 46 night fcd, Rn > 0 daytime, Rs/Rso limits (my finding 2: +9%), local-date and complete-day dailyET (G4, R8, R9), and the overpass Rs/Rso QC (my finding 3).
3. **A robust calcH.** Damped iteration, convergence flag, continuous u200 ≥ 2 m/s clamp, anchor validation (G3, G13, G9).
4. **Silent-input safety.** Bands by name (G14), Kelvin check (R13), `columns` by name and unit checks (R5, R6), writeResults corruption (G1), typo'd methods (G2).
5. **Outputs a practitioner can use correctly.** ETrF/ET24 from `ET24h(waterLSEB)`, the ETrF → grass-Kc note with daily ratios (R22, R24), and roughness guidance for vines and orchards with the z0m floor (R29, my finding 11).
   - Next: C_rad for slopes, ETrF.hot (Eq. 47), MODIS albedo, ST_B10 as the only L8/L9 thermal default.
