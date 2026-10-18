# R — discussion of S's and G's reviews (practitioner view)

Update to my report: my full-scene run had not finished METRIC.EB after about 50 h (16.8 GB working set, with other reviewers' jobs on the same machine), and I stopped it. G measured about 58 min and a 33 GB peak on the same machine. For a 16 GB laptop, the practical answer is "crop to 30–60 km". I raised finding R31 to Medium.

## 1. Their findings that affect users

S's findings, ranked by how much they would mislead an irrigation advisor:

- **S1 Anchor windows (high) — Agree; this is the top scientific issue.** It explains my ETrF > 1.05 observations (R18). One anchor choice moving mean ET24 by 31% is larger than any irrigation tolerance. I support the scene-relative defaults and the diagnostics in `print()`.
- **S2 ASCE daily ETr +9% on clear days (high) — Agree, high.** It is invisible to users and passes 1:1 into every map, on exactly the days we use. It is cheap to fix.
- **S3 Rs/Rso check at the overpass — Agree, but raise it to high.** It is the same idea as my R6: a cloud-affected or badly calibrated pyranometer is the most common station fault in the field. The warning must say what to do ("check the pyranometer or use clear-sky ETr").
- **S4 SR outside 0–1 masked — Agree; priority as S has it.** It is a two-line fix with a 4% effect.
- **S6 ET at the hot pixel after rain — Agree, medium.** For advisors this is real: images after spring rains are common. A simple warning ("rain X mm in the last N days") is enough for v1; the FAO-56 Ke helper can come later.
- **S7 / G3 / G13 low-wind oscillation and the u200 floor — Agree; this is one issue across three findings.** Both Talca vignettes are affected (my R12). A 21% jump in LE between 0.50 and 0.53 m/s of wind is unacceptable for operational use.
- **S5 τNB and the SC method — Agree for scientific correctness, lower priority for practitioners.** ET survives through the calibration (±0.04 ETrF), and L8/L9 users get ST_B10. Recommend ST_B10 as the default in the docs.
- **S8 C_rad — Agree, medium, but only for the mountain mode.** Most of our valleys are flat. Document plain vs DEM first (R28).
- **S9 = G4 = my R9 — the same dailyET day bug found three times. Merge them.**
- **S10, S16–S20 — Agree, low; no user impact worth delaying anything else.**
- **S11 Z.om floor and Eq. 33 for vines and trees — Agree.** It matches my R29, where Perrier changes ET by 12% for vineyards. Add the floor and a crop-type table.

G's findings:

- **G1 writeResults corrupts results (high) — Agree, high.** Users with large scenes will turn this option on precisely because of G5. Until it is fixed, the docs should say "do not use writeResults".
- **G2 typos return stale global objects (high) — Strongly agree.** A new user typing `"olmedo"` or `"tasumi"` is the most likely mistake of all. Getting a *silently stale* G or Ts from the workspace is worse than any error.
- **G5 performance — Agree.** I would ask for the 30–60 km AOI recommendation in the docs first, then the faster code.
- **G6 = my R14 (global aoi).** G rates it Medium, I rated it High. I keep High: it breaks the first line of the Landsat 8 vignette and sends users to finding R16 (an AOI the size of a farm).
- **G7 getwd() MTL search — Agree.** It is the root cause of my R9 (dailyET without MTL). Requiring explicit paths makes the workflow safer for users.
- **G8 unsorted or duplicated records — Agree, medium-high for practice.** Excel exports and logger merges often have both. A *silently wrong* at.sat is high for users.
- **G9 hot and cold anchors swapped accepted — Agree.** Manual anchors are exactly what my R19 asks users to set, so they need this validation.
- **G10 / G11 par() changes, printing and plots.** These **help** new users only if the change keeps a visual check of the anchors. See conflict C1.
- **G12 Level-1 fill not masked — Agree, high for L7 and L1 users.** Fill pixels with Ts = −53 K and plausible Rn end up in the maps.
- **G14 bands by position (LAI 1.25 → 0.001) — Agree, medium.** Users who build their own stacks from GEE or QGIS will hit it.
- **G17 empty layer names, G19 examples — Agree.** They overlap with my R22 (ET24h example, layer names).

## 2. Conflicts and proposed resolutions

- **C1. Plots and printing (G: plots off and silent by default; my view: users need to see the anchors).**
  Resolution: `plots = interactive()` and `message()` gated by `verbose`, plus a **one-line QC summary that is always printed**, using `message()` so it can be suppressed:
  "anchors: cold 309.1 K / hot 317.0 K (ΔT 7.9 K); ETrF > 1.1 in 3% of pixels; u200 used 2.8 m/s; converged in 7 iterations".
  Add a `checkAnchors()` / `plot(x, what = "anchors")` method with an RGB background for the visual check.
- **C2. Automatic scene-relative anchors (S1) vs anchor checking by hand (R19).** They do not conflict. Better automatic defaults reduce harm; diagnostics and the manual `anchors =` route remain the safety net. Both should ship together, so that new defaults do not suddenly change maps without users being told: add a NEWS entry and give a before/after comparison in the vignette.
- **C3. The u200 fix.** S wants a continuous floor plus damping; G wants a clamp at the threshold. Agree on: damping, a convergence flag, `u200.min` as an argument with a value S documents, and a plain-language warning (my R12 text).
- **C4. Default changes that alter results** (S2 ASCE, S4 SR mask, S11 Z.om floor, S1 windows, `alb.coeff` default, `n` default). Advisors comparing seasons need stability. Ship them in one release, together with a "results changed" table in NEWS (mean ET24 before and after for both examples) and `$methods` recording the version and options.
- **C5. Severity of the global `aoi` (G6 Medium vs R14 High).** I propose High, because of its downstream effect on R16.
- **C6. Performance vs AOI size.** G's speed-ups help, but scientifically (R16, S1) the AOI must still be large enough. Document a recommended AOI of 30–60 km. That is also the size that runs in a few minutes today.

## 3. Updated team top 5

1. **Fail loudly on bad inputs.** `match.arg` for every method (G2). Weather-station checks: `columns` matched by name, units and ranges, Rs/Rso, overpass time, sorted and unique records (R5, R6, R7, S3, G8). Anchor validation (G9), surf.temp in Kelvin (R13), and masking of L1 fill and invalid SR (G12, S4).
2. **Correct daily reference ET.** ASCE Eq. 46 compliance (S2). The overpass's local day, honouring `date`, and a 24-hour completeness check (S9, G4, R8, R9). No getwd() MTL searching (G7). `ET24h()` returns ETrF (and Kc on a grass basis), with a documented conversion (R22, R24).
3. **Calibration you can trust and see.** Scene-relative, configurable anchor windows (S1). A minimum AOI and hot–cold ΔT guard (R16). A cloud-mask warning (R15). A damped, convergent H iteration with a documented u200 floor (S7, G3, G13). The always-on QC summary, `checkAnchors()`, and guidance on roughness for vines and orchards (R18, R19, R29, S11).
4. **No hidden state.** Remove the global `aoi` (R14, G6) and fix writeResults file naming (G1). Restore `par()`, and use `message()` / `plots = interactive()` (G10, G11).
5. **A practitioner path, in Spanish too.** Fix the README (CRAN archived) and make the Landsat 8/9 vignette the entry point. Ship the download guides with AOI-cropped downloads. Add a per-field ET/Kc table, GeoTIFF export, and "When not to trust" and "Choosing the station" boxes, with the recommended AOI size and run time (R1–R4, R25, R27, R28, R31, G5). Then add tests for the uncovered paths and examples for CRAN resubmission (G18, G19).
