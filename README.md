# LVL Project

This folder contains the current LVL seismic-refraction workflow for first-break picking, layer fitting, depth estimation, and export/report generation.

## Current entry points

- Desktop GUI: `lvl/apps/lvl_studio.py`
- CLI / legacy-full workflow shell: `lvl/apps/lvl_refraction.py`
- Shared backend for both: `lvl/src/`

## Current architecture

- `apps/`
	- thin entry points and UI shells
	- should orchestrate, not own geophysical algorithms
- `src/picker/`
	- first-break picking pipeline
	- legacy-style helpers plus Picker Engine V2
- `src/refraction/`
	- corrected travel times, layer fitting, velocity/depth analysis
- `src/io/`
	- SEG2 reading, project persistence, picks/session persistence, exports
- `src/utils/`
	- geometry helpers and math utilities shared across fronts ends
- `src/common/`
	- shared paths, constants, app metadata

## Start here

- If the problem is GUI behavior, start in `lvl/apps/lvl_studio.py`.
- If the problem is auto-pick behavior, start in `lvl/apps/lvl_studio.py` and then trace into `lvl/src/picker/`.
- If the problem is exported Excel/TXT/plots/reports, start in `lvl/src/io/exporters.py` and `lvl/src/refraction/pipeline.py`.
- If the problem is geometry/profile discovery, start in `lvl/src/utils/geometry.py` and `lvl/src/io/seg2_reader.py`.

## Critical dependency chain

- `apps/lvl_studio.py`
	- depends on `src/common.paths`, `src/common.settings`
	- depends on `src/io.project_io`, `src/io.pick_reader`, `src/io.seg2_reader`, `src/io.exporters`
	- depends on `src/picker.preprocessing`, `src/picker.features`, `src/picker.refinement`, `src/picker.picker`, `src/picker.settings`
	- depends on `src/refraction.pipeline`, `src/refraction.layer_analysis`, `src/refraction.velocity_model`
- `src/picker/picker.py`
	- depends on `src/picker.preprocessing`, `src/picker.features`, `src/picker.likelihood`, `src/picker.coherence`, `src/picker.optimizer`, `src/picker.confidence`, `src/picker.quality`
- `src/refraction/pipeline.py`
	- depends on `src/io.exporters`, `src/refraction.processing`, `src/refraction.layer_analysis`, `src/refraction.velocity_model`

## Picker Engine V2 settings exposed in the GUI

These are opened from `LVL Studio -> Picking -> Picker V2 Settings...`.

- Feature weights, range `0.0 .. 5.0`
	- `hilbert_weight`
	- `stalta_weight`
	- `aic_weight`
	- `gradient_weight`
	- `energy_weight`
	- `snr_weight`
	- `kurtosis_weight`
	- `skewness_weight`
	- Meaning: relative influence in the fused per-trace arrival likelihood.
- Processing
	- `sta_window`, range `0.1 .. 100.0 ms`
	- `lta_window`, range `0.5 .. 500.0 ms`
	- `hilbert_onset_pct`, range `0.0 .. 1.0`
- Coherence
	- `coherence_weight`, range `0.0 .. 1.0`
	- `coherence_radius`, range `1 .. 10` traces per side
	- `coherence_align`, boolean
	- `coherence_max_shift`, range `1 .. 200` samples
- Velocity gate
	- `use_velocity_gate`, boolean
	- `vmin_m_s`, range `1 .. 20000 m/s`
	- `vmax_m_s`, range `1 .. 20000 m/s`
	- `gate_pad_ms`, range `0 .. 500 ms`
- Path/confidence
	- `smoothness_penalty`, range `0.0 .. 5.0`
	- `jump_penalty`, range `0.0 .. 5.0`
	- `minimum_confidence`, range `0.0 .. 1.0`

## What is safe to replace

- UI layout/details in `apps/lvl_studio.py` can be changed without touching algorithms if public payloads stay the same.
- Picker algorithms should be changed in `src/picker/`, not duplicated back into `apps/`.
- Export layout can be changed in `src/io/exporters.py` if `src/refraction/pipeline.py` contracts are preserved.

## What is critical not to duplicate

- SEG2 reading logic: keep in `src/io/seg2_reader.py`
- pick/session JSON IO: keep in `src/io/pick_reader.py`
- full picker V2 pipeline: keep in `src/picker/`
- layer/velocity calculations: keep in `src/refraction/`

## Formula provenance and references

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

## Results export