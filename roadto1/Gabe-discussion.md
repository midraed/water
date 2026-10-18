# Gabe: discussion round

Effort to implement: S = under 1 day, M = a few days, L = a redesign or a week or more.

## 1. Code-related findings from Sam and Rosa

**Sam**
- **S1, anchor criteria hard-coded and site-specific:** agree, high, effort L.
  - Do this together with my finding 5 (calcAnchors rewrite). The new design would compute, for each pixel, the relaxation it needs. That makes configurable, scene-relative (percentile) criteria and relaxation in Kelvin almost free.
  - Put the criteria in one exported `anchorCriteria()` list instead of six hard-coded copies (my finding 21).
- **S2, ASCE hourly and daily reference ET:** agree, high, effort M.
  - The changes are pure scalar functions and easy to unit-test against ASCE Appendix C.
  - Overlaps my finding 22 (273.16).
- **S3, weather QC at the overpass:** agree, medium, effort M.
  - Overlaps Rosa R6 and R7. Use one `checkWS()` validator called inside `read.WSdata()`.
- **S4, Collection 2 SR valid range:** agree, raise to high, effort S. It is two lines in `loadImageSR()`, the same kind of fix as my finding 12 (DN = 0 fill in `loadImage()`).
- **S5, τNB and the "SC" method:** agree; this is the same issue as my finding 15. Effort S to document, M to implement Jiménez-Muñoz.
- **S6, hot-pixel ET always 0; S8, C_rad never computed; S12, other mountain-mode issues:** agree, medium, effort M.
  - Mountain mode has 0 % test coverage (my finding 18). Add a DEM_Talca test before touching it.
- **S7, low-wind u200 floor and oscillation:** this is my findings 3 and 13, and I agree. Effort M: damping, a `converged` flag, and a continuous floor. My parity demonstration (LE 199 vs 208 W/m2 depending on `maxit`) supports it.
- **S9, `dailyET()` uses the UTC date:** this is my finding 4 and Rosa R8 and R9. Agree, high, effort S.
  - The fix has three parts: select by the local date of `at.sat`, stop on a missing day, and warn when fewer than 24 hours are present.
- **S10, MODIS albedo; S16-S18, minor constants:** agree, low, effort S.
- **S11, z0m minimum 0.005 m:** agree, effort S. It also removes the `u* == 0 -> 0.1` hack in calcH.
- **S13, energy-balance closure; S14, flexible ETrF update trigger:** agree, low.
- **S20, `calcH(method = "lm")` uses only the last pair:** agree. Either deprecate it or fit the regression properly. Effort S to deprecate.

**Rosa**
- **R5, names in `columns` are ignored:** agree, high, effort S. It is a silent wrong result, the same class as my finding 2. Match by name when names are present; otherwise check that the column count and types are plausible.
- **R6/R7/R13, unit, time-zone and °C checks:** agree, effort M, through the same validator as S3.
  - Ts in °C is a one-line `stop()`: Ts minimum < 200 K.
- **R10, cryptic crashes:** this is my findings 2 and 9. Effort S-M: `match.arg` plus anchor validation.
- **R11, csv2 file read with `read.csv`:** agree, effort S. Detect `;` in the header line and suggest `read.WSdata2`.
- **R12 and R20, noisy and mislabelled output:** these are my findings 11 and 13.
- **R14, hidden global `aoi`:** this is my finding 6. Agree, high, effort S. Remove it rather than message about it.
- **R15, cloud warning; R16, small AOI; R18, ETrF QC summary:** agree, effort M. See the conflicts below on how to show this.
- **R17, default `n` = 1 while every example uses 5:** agree, effort S. My timing data argue for keeping the default small until calcAnchors is fast.
- **R22, ET24h interface; R25, export:** agree, effort S-M.
  - `ET24h()` should accept a `waterLSEB` object and return ET24 and ETrF with named layers. This also fixes the empty layer names (my finding 17).
- **R26, seasonal ET; R31, progress information:** effort M/L, later.
  - Progress reporting only becomes meaningful after the performance work in my finding 5. One scene currently takes 58 min with no output.
- **R1-R4, R21, R23, R24, R27-R30, documentation and vignettes:** agree. These are not code risks. They pair naturally with the missing `\examples` needed for CRAN (my finding 19).

## 2. Conflicts between perspectives and a proposed resolution

1. **More diagnostics versus less console output.** Sam and Rosa want calibration diagnostics and a QC summary; Rosa and I want quieter output.
   - Resolution: do not print diagnostics. Return them in the result: `result$diagnostics` (anchor values, ETrF percentiles, share of ETrF > 1.1, ΔT between hot and cold, the `converged` flag, flex) and `result$warnings`.
   - `print.waterLSEB()` and `summary()` show them.
   - Real problems become `warning()`. Progress uses `message()`, gated by `verbose`.
   - Plots go behind a `plots` argument (default `interactive()`), plus a separate `plotAnchors()` helper for Rosa's R19.
2. **New arguments versus API consistency.** New criteria, Kt, QC thresholds and convergence settings would bloat METRIC.EB, which already has 23 arguments.
   - Resolution: group them into control lists with constructor functions, such as `anchorControl()` and `calcHControl(maxit, damping)`, and replace the untyped `extraParameters` vector.
   - Then fix the argument order once (image.SR first in all three METRIC wrappers) in the same release, with deprecation warnings.
3. **Strict validation (`stop`) versus tolerant use.** My view and Rosa's say stop on bad input; Sam's clear-sky ETr option would let questionable data through.
   - Resolution: stop on impossible values (units, Kelvin, missing day, NA or swapped anchors). Warn on unusual but possible values (Rs/Rso, ETrF > 1.1, low wind), and let an explicit option such as `ETr.inst = "clear-sky"` change the behaviour.
4. **Changing the science versus reproducibility.** Sam's S1, S2 and S4 change results by 4-31 %.
   - Resolution: first add regression tests that pin the current outputs. Then change the science one finding at a time, recording each change and its measured effect in NEWS, as the 0.8.5 NEWS already does.
5. **Performance versus scientific changes to anchors.** Do not tune the anchor criteria inside the current 97-pass code. Rewrite `calcAnchors()` first (behaviour-preserving, with tests), then apply Sam's criteria.

## 3. Updated top 5 for the team

1. **Stop silent wrong results from bad input.**
   - `match.arg` everywhere (my finding 2), `columns` matched by name (R5).
   - Remove the global `aoi` (my finding 6, R14) and the `getwd()` MTL lookup (my finding 7).
   - Sort and de-duplicate the weather records (my finding 8).
   - Mask SR out of range and DN = 0 fill (S4, my finding 12).
   - Range and unit validator for weather data and Ts (R6, R7, R13, S3).
   - Mostly effort S, high value.
2. **Correct and complete daily ET.**
   - Local overpass day and the 24-hour check (my finding 4, S9, R8, R9).
   - ASCE-compliant hourly and daily ETr (S2).
   - ET24h returns ET24 and ETrF with named layers (R22).
3. **Robust calibration.**
   - Damped r_ah iteration with a `converged` flag, a continuous u200 floor and the z0m minimum (my findings 3 and 13, S7, S11).
   - Anchor validation (my finding 9).
   - Diagnostics returned in `result$diagnostics`, not printed (resolution 1).
4. **Fix writeResults corruption, then do the performance rewrite.**
   - Unique output files (my finding 1).
   - calcAnchors rewritten as vectorised and single-pass, followed by scene-relative criteria (my finding 5, then S1).
   - calcH iterating on the anchors plus one `lapp()` pass, so that a scene takes minutes instead of about 1 hour.
5. **Tests and CRAN readiness alongside a practitioner vignette.**
   - Regression tests for the main paths (L7, mountain mode, ET24h, dailyET with `date`).
   - Examples for about 24 exported functions; `par()` restored and output silenceable (my findings 10, 11, 18, 19).
   - The Landsat 8/9 vignette as the entry point, including a per-field table and export (R1, R2, R25).
