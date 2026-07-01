# LVL Refraction Analysis

This document describes what the LVL refraction workflow does, in which order,
and with which equations. The goal is to make every processing and inversion
step transparent and reproducible.

Primary script:
- lvl/scripts/lvl_refraction.py

Data folders:
- lvl/data/<profile>/ for SEG2 files
- lvl/data/geometry100.txt and lvl/data/geometry200.txt for receiver geometry
- lvl/data/field_report*.xls(x) for optional default offsets and geometry hints

Output folders:
- lvl/output/<profile>/

## 1) End-to-end processing flow

1. Read SEG2 shot files and sort by FFID/file number.
2. Load receiver geometry from geometry100 or geometry200.
3. Resolve shot positions:
	- from config auto/manual
	- optionally from SEG2 SOURCE_LOCATION (depending on policy)
4. Build absolute time axis using SEG2 delay header and sample interval.
5. (Display/picking) apply optional zero-phase filter and gain.
6. Interactive first-break picking per shot.
7. Apply bulk static shift to picks (analysis/export stage).
8. Build corrected offsets per pick:
	- inline offset
	- optional inline shift for far-offset acquisition geometry
	- perpendicular correction to get true offset
9. Interactive layer-window fitting (per shot side L/R or ALL).
10. Compute velocities, intercept times, and layer depths.
11. Compute observed-vs-computed RMS.
12. Export text, Excel, and QC plots.

## 2) Geometry and shot-position model

Receiver position for trace $k$ is $x_k$ (m along profile).

Shot position is $SP$ (m along profile).

Signed inline offset:
$$
\Delta x_k = x_k - SP
$$

Absolute inline offset:
$$
X_{inline,k} = |\Delta x_k|
$$

Optional inline shift (for far-offset geometry cases):
$$
X_{corr,k} = \max(0, X_{inline,k} + dX)
$$
where $dX$ is per-shot inline shift in meters.

Perpendicular shot-to-line offset is $PO$ (m), then corrected true offset is:
$$
X_{O,k} = \sqrt{X_{corr,k}^2 + PO^2}
$$

Important: velocity fitting in the current implementation is based on
$X_{O,k}$ (true offset), not on raw inline distance.

## 3) Time axis and pick timing

Given:
- sample interval $dt$ (s)
- SEG2 delay $t_{delay}$ (ms)

Absolute sample time:
$$
t_n = t_{delay} + n \cdot dt \cdot 1000
$$

Picks are stored as absolute first-break time in ms from trigger-reference axis.

Bulk static correction is added at analysis/export stage:
$$
t_{bulk,k} = t_{raw,k} + t_{static}
$$

## 4) Filter and gain (display/picking domain)

### 4.1 Ormsby zero-phase bandpass

Four-corner trapezoid: $f1-f2-f3-f4$.

Frequency response $H(f)$:
- ramp up from 0 to 1 on $[f1,f2]$
- passband 1 on $(f2,f3]$
- ramp down from 1 to 0 on $[f3,f4]$

Applied in frequency domain with optional zero padding.

### 4.2 Butterworth zero-phase bandpass

Alternative display/picking filter using forward-backward SOS filtering.

### 4.3 Gain

- none
- norm: per-trace normalization by max absolute amplitude
- agc: windowed RMS/mean envelope normalization

These affect picking visibility, not the geometric equations themselves.

## 5) Piecewise T-X fitting

For each shot side (L/R) or combined side (ALL), user defines windows in offset.

Each layer segment is linear:
$$
t(X) = mX + b
$$

where:
- $m$ is slope in ms/m
- $b$ is intercept time in ms

Linear regression is done with polyfit on picks inside each selected window.

Velocity is:
$$
V = \frac{1000}{m}
$$

Coefficient of determination:
$$
R^2 = 1 - \frac{\sum (t_i-\hat t_i)^2}{\sum (t_i-\bar t)^2}
$$

## 6) RMS fit quality

For observed picks $t_i$ and computed picks $\hat t_i$ from fitted segments:
$$
RMS_{ms} = \sqrt{\frac{1}{N}\sum_{i=1}^{N}(t_i-\hat t_i)^2}
$$

This is shown in fit review and exported to analysis output.

## 7) Intercept-time depth formulas

### 7.1 Two-layer depth to first refractor

Given direct-wave velocity $V1$, refractor velocity $V2$, and intercept $ti1$:
$$
h1 = \frac{(ti1/1000) \cdot V1 \cdot V2}{2\sqrt{V2^2 - V1^2}}
$$

Valid only if $V2 > V1 > 0$.

### 7.2 Three-layer depth to second refractor

With third-layer velocity $V3$, second intercept $ti2$, and $h1$ above:

Critical-angle cosine term:
$$
\cos(i_{c12}) = \sqrt{1 - (V1/V2)^2}
$$

Delay correction from first layer:
$$
ti2_{eff} = ti2 - 2h1\cos(i_{c12})/V1 \cdot 1000
$$

Depth to second refractor:
$$
h2 = \frac{(ti2_{eff}/1000) \cdot V1 \cdot V3}{2\sqrt{V3^2 - V1^2}}
$$

Valid only for physically consistent velocity ordering and positive effective
intercept time.

## 8) Excel and report integration

Field report integration can auto-load:
- per-shot perpendicular offsets (default column D)
- optional inline shift column
- profile grouping via LVLXXX labels (vertical block parsing)

Geometry selection behavior:
- explicit CLI geometry wins
- if no explicit geometry: script can infer 100 from nearby short-profile note
- otherwise defaults to 200

Recommended for strict reproducibility:
- always pass geometry explicitly when comparing historical results.

## 9) Outputs and how to read them

Main outputs in lvl/output/<profile>/:

- <profile>_picks_clean.txt
  - exported picks table for external tools/notebooks

- <profile>_picks.xlsx
  - Config, per-shot tables, Combined, Layer_Picks, Layer_Averages, Analysis
  - includes velocities, intercepts, depths, RMS

- <profile>_tx_picks.png
  - T-X display with fitted segments

- <profile>_fit_rms.png
  - observed vs computed fit summary

- <profile>_corrected_qc.png
  - corrected pick QC across shots

## 10) Command examples

From repo root:

1. Explicit geometry 200 (recommended for historical comparability):
	- python lvl/scripts/lvl_refraction.py 150 --geom 200

2. Explicit geometry 100:
	- python lvl/scripts/lvl_refraction.py 150 --geom 100

3. Use field report defaults for offsets:
	- python lvl/scripts/lvl_refraction.py 150 --geom 200 --perp-excel lvl/data/field_report_07-04-2026_to_11-04-2026.xls

4. Export-only re-analysis without new picking:
	- python lvl/scripts/lvl_refraction.py 150 --export-only --geom 200

## 11) Reproducibility checklist

When comparing runs (script versions, Excel checks, external tools), ensure all
of the following are identical:

1. Geometry file (100 vs 200) and shot positions.
2. Per-shot PO and inline shift values.
3. Bulk static value.
4. Pick set used (same picks.json/session).
5. Fit windows selected for each shot side.
6. Segment basis uses corrected true offset XO.
7. Layer count and velocity-order validity.

If one of these differs, velocity and depth results can change significantly.

## 12) Worked numeric mini-example

This mini-example shows the same equations with concrete numbers.

Assume one shot side has corrected true-offset picks approximately on:
$$
t(X_O) = 25 + 0.55X_O
$$
with $t$ in ms and $X_O$ in m.

Then:
- slope $m = 0.55$ ms/m
- velocity:
$$
V = \frac{1000}{m} = \frac{1000}{0.55} = 1818.18\,\text{m/s}
$$

For depth, assume a 2-layer case with:
- $V1 = 600$ m/s
- $V2 = 1800$ m/s
- $ti1 = 40$ ms

Depth to first refractor:
$$
h1 = \frac{(0.040)\cdot 600\cdot 1800}{2\sqrt{1800^2-600^2}}
	= 12.73\,\text{m}
$$

Now assume a 3-layer extension:
- $V3 = 2600$ m/s
- $ti2 = 85$ ms

Critical-angle cosine term:
$$
\cos(i_{c12}) = \sqrt{1-(600/1800)^2} = 0.9428
$$

Effective second intercept time:
$$
ti2_{eff} = 85 - 2\cdot 12.73\cdot 0.9428/600\cdot 1000
			= 45.0\,\text{ms}
$$

Depth to second refractor:
$$
h2 = \frac{(0.045)\cdot 600\cdot 2600}{2\sqrt{2600^2-600^2}}
	= 13.87\,\text{m}
$$

This is exactly the same equation path implemented in
lvl/scripts/lvl_refraction.py.

## 13) Formula provenance and references

The formulas used here are standard seismic refraction relationships, not
project-specific inventions.

Primary concepts behind the equations:
1. Snell's law and critically refracted head-wave travel-time relationships.
2. Linear time-distance segment fitting for apparent velocity.
3. Intercept-time depth equations for layered media.
4. RMS misfit as least-squares fit quality metric.

Recommended references (textbooks and standard geophysics sources):
1. Sheriff, R. E., and Geldart, L. P. (1995). Exploration Seismology (2nd ed.). Cambridge University Press.
2. Telford, W. M., Geldart, L. P., and Sheriff, R. E. (1990). Applied Geophysics (2nd ed.). Cambridge University Press.
3. Kearey, P., Brooks, M., and Hill, I. (2002). An Introduction to Geophysical Exploration (3rd ed.). Blackwell Science.
4. Dobrin, M. B., and Savit, C. H. (1988). Introduction to Geophysical Prospecting (4th ed.). McGraw-Hill.
5. Yilmaz, O. (2001). Seismic Data Analysis. Society of Exploration Geophysicists.

Notes on mapping formula-to-implementation:
1. Velocity from slope and line fitting: see fit sections and exported per-shot slope/intercept in the script outputs.
2. Two-layer and three-layer intercept-time depth equations: implemented in depth_2layer and depth_3layer in lvl/scripts/lvl_refraction.py.
3. RMS equation: implemented in fit review and analysis export paths in lvl/scripts/lvl_refraction.py.
