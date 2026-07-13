# LVL Refraction Instructions (User-Focused)

This guide is for day-to-day use and result interpretation, including the formulas behind the exported summary.

Main script:
- lvl/scripts/lvl_refraction.py

Main consolidated summary output:
- lvl/output/lvl_velocity_summary.xlsx

## 1) Summary Table Layout (copy-friendly)

The first columns are exported in this order:

1. NUM
2. NAME
3. SP
4. SPREAD
5. X_UTM
6. Y_UTM
7. ELEV
8. V0E
9. V0C
10. V0
11. H0
12. V1
13. TI1
14. H1
15. V2
16. TI2
17. DR

Then extra fields are appended to the right, for traceability:

- LINE_LENGTH_SP_M
- LINE_LENGTH_M
- MIDDLE_STATION
- ACQ_DATETIME_DE
- Side-specific velocities (left/center/right)

### Field meaning

- NUM: Consecutive row number.
- NAME: LVL label, e.g. LVL150 or LVL150_new.
- SP: Number of source points (typically 3).
- SPREAD: Line length text in meters (from coordinate geometry where available).
- X_UTM, Y_UTM, ELEV: Middle-station coordinates and elevation.
- V0E: Off-end low velocity average from left and right off-end shots.
- V0C: Central-shot low velocity.
- V0: Average low velocity.
- H0: Thickness of low velocity layer.
- V1, TI1, H1: Intermediate layer velocity, intercept time, and thickness.
- V2, TI2: Refractor velocity and intercept time.
- DR: Refractor depth (H0 + H1).

## 2) Length Definitions

Two lengths are tracked:

1. LINE_LENGTH_SP_M:
- Derived from active processing shot-point geometry (receiver span used in processing).
- Typical values: around 94 m for 100 m geometry and around 192.5 m for 200 m geometry.

2. LINE_LENGTH_M:
- Derived from LVL coordinate workbook station/X/Y polyline distance.
- Used as SPREAD in the main summary columns.

## 3) Coordinate and Date Source Priority

Coordinates and date/time used for summary are resolved with this logic:

1. SEG2 coordinate search (middle shot):
- Try extracting X/Y/Z style fields from SEG2 headers.
- If unavailable, continue with Excel fallback.

2. LVL coordinate Excel fallback:
- Profile is matched by LVL profile token.
- Repeated station rows use the last row (last wins).
- Midpoint station coordinate is used for X_UTM, Y_UTM, ELEV.

Current export behavior in script:
- Coordinates/elevation are taken from coordinate Excel first.
- SEG2 coordinates are used only if that profile is not found in coordinate Excel.

3. Date/time priority:
- Prefer date/datum from the same Excel sheet used for coordinates.
- If missing, use SEG2 middle-shot acquisition time.
- Time is formatted in German timezone (Europe/Berlin, CET/CEST).

Note:
- SEGY coordinate fallback is not currently part of this script pipeline.

## 4) Geometry Rule of Thumb

If no explicit geometry is forced by CLI:

- If inferred spread is under 100 m, processing tends to resolve to 100 m geometry.
- If spread is around 192.5 m, processing tends to resolve to 200 m geometry.

For strict reproducibility in comparisons, pass geometry explicitly:

- python lvl_refraction.py 120 --geom 100
- python lvl_refraction.py 150 --geom 200

## 5) Core Equations

True offset:

```text
XO = sqrt(Xcorr^2 + PO^2)
```

Velocity from fitted slope m (ms/m):

```text
V = 1000 / m
```

Two-layer depth:

```text
H0 = (ti1/1000) * V0 * V1 / (2 * sqrt(V1^2 - V0^2))
```

Three-layer second depth (spreadsheet-compatible form):

```text
TI2_eff = TI2 - 2 * H0 * cos(arcsin(V0/V2)) / V0 * 1000
H1      = V1 * TI2_eff / (2000 * cos(arcsin(V1/V2)))
```

Refractor depth:

```text
DR = H0 + H1
```

Why 2000 appears:

- 1000 from ms to seconds conversion.
- 2 from intercept-time geometry factor.

## 6) Practical Validation Checklist

If values differ from historical Excel:

1. Confirm same geometry (100 or 200).
2. Confirm same PO and inline shift by shot.
3. Confirm same picked first breaks.
4. Confirm same layer windows and fit acceptance.
5. Confirm SPREAD and coordinate source are from the expected workbook.
6. Check whether NAME or NAME_new row is the current one.

## 7) End-to-End Processing Flow

1. Read SEG2 shot files and sort by FFID/file number.
2. Load receiver geometry (100 m or 200 m setup).
3. Resolve shot positions from config/header policy.
4. Build absolute time axis from SEG2 delay and sample interval.
5. Apply display filter/gain for picking view.
6. Pick first breaks interactively.
7. Apply bulk static correction to picks for analysis/export.
8. Build corrected offsets (inline, optional shift, perpendicular).
9. Fit layer windows per side (L/R/ALL).
10. Compute velocities, intercept times, and depths.
11. Compute RMS diagnostics.
12. Export profile files and consolidated summary.

## 8) Geometry and Offset Model

Receiver position for trace k is x_k (m), shot position is SP (m).

Signed inline offset:

```text
DX_k = x_k - SP
```

Absolute inline offset:

```text
X_inline_k = abs(DX_k)
```

Optional inline shift:

```text
X_corr_k = max(0, X_inline_k + dX)
```

True offset with perpendicular distance PO:

```text
XO_k = sqrt(X_corr_k^2 + PO^2)
```

Fitting is based on true offset XO.

## 9) Time Basis

Sample time axis (ms):

```text
t_n = delay_ms + n * dt_s * 1000
```

Bulk static correction for analysis/export:

```text
t_bulk = t_raw + t_static
```

## 10) Piecewise T-X Fitting

Each picked layer window is fitted as:

```text
t(X) = mX + b
```

Velocity from slope m (ms/m):

```text
V = 1000 / m
```

R-squared:

```text
R2 = 1 - sum((t_i - t_hat_i)^2) / sum((t_i - t_mean)^2)
```

## 11) RMS Quality Metric

RMS_ms = sqrt(mean((t_obs - t_calc)^2))

```text
RMS_ms = sqrt(mean((t_obs - t_calc)^2))
```

This is used in fit-review and diagnostics exports.

## 12) Output Files

Main outputs in lvl/output/<profile>/:

- <profile>_picks_clean.txt
- <profile>_picks.xlsx
- <profile>_tx_picks.png
- <profile>_fit_rms.png
- <profile>_corrected_qc.png
- lvl_velocity_summary.xlsx

## 13) Common Commands

- python lvl/scripts/lvl_refraction.py 150 --geom 200
- python lvl/scripts/lvl_refraction.py 150 --geom 100
- python lvl/scripts/lvl_refraction.py 150 --export-only --geom 200
- python lvl/scripts/lvl_refraction.py 150 --coord-excel lvl/data/LVL_coordinates.xlsx --device-type sw_maps

## 15) GUI Command Center

Launch the standalone GUI:

- python lvl/scripts/lvl_command_center.py

What it provides:

1. Central command panel with Run/Stop/Open Output.
2. File/folder selectors for:
	- SEG2 profile folder
	- field report Excel
	- coordinate Excel
3. Device type selector (`sw_maps` default, `geomax` optional).
4. Runs the same standalone processing script with selected options.
5. Sends live picker commands over JSON bridge while picking windows are open.

Program structure (simplified):

1. fb_picker.py: dedicated picker reference/experiments.
2. lvl_refraction.py: main processing and analysis engine.
3. lvl_command_center.py: single central GUI launcher and live control panel.
4. lvl_modules/: shared contracts used by GUI and backend.
	- run_config.py: run-request model and command construction.
	- control_bridge.py: JSON live-control read/write protocol helpers.
	- app_paths.py: shared project/data/output path definitions.

## 14) Reproducibility Checklist

For comparable results, keep these fixed between runs:

1. Geometry choice (100 vs 200) and shot positions.
2. PO and inline shift values.
3. Bulk static setting.
4. Pick set (same session/final picks).
5. Fit windows and accepted sides.
6. Offset basis (true offset XO).
7. Layer count and physically valid velocity order.
