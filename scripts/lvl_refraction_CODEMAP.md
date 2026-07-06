# lvl_refraction.py - Full Function Documentation

This document explains why each major function/class exists in `lvl_refraction.py` and where it fits in the workflow.

The script is organized as:
1. Environment and dependencies
2. Configuration
3. Geometry/report parsing helpers
4. Offset model setup UI (matplotlib)
5. Signal processing and picking analysis
6. Interactive picker class
7. Persistence helpers (JSON)
8. Export functions (TXT/Excel/plots)
9. Main pipeline + CLI

## 1) Environment and setup

- `_ensure(pkg, mod="")`
  - Why: keeps the script self-contained by importing or auto-installing missing Python packages in the active environment.

- `_tc()`
  - Why: central theme palette provider for all matplotlib views, ensuring visual consistency.

- `_ensure_interactive_backend()`
  - Why: guarantees a GUI backend for matplotlib when interactive windows are required.

## 2) Geometry and profile metadata

- `load_geometry(geom_type)`
  - Why: reads receiver coordinates from geometry text files (`geometry100.txt`, `geometry200.txt`).

- `auto_shot_positions(recv_pos)`
  - Why: derives standard shot positions from spread geometry (before spread, midpoint, after spread).

- `_excel_col_to_index(col)`
  - Why: converts user CLI Excel column labels (A, D, F, AA) into zero-based numeric indices.

- `_normalize_profile_token(value)`
  - Why: normalizes profile IDs (e.g. `LVL150`, `150`, `150_`) for robust matching across files.

- `discover_field_report_excels(data_dir)`
  - Why: finds all field-report Excel files; used to iterate across reports for profile lookup.

- `discover_profile_folders(data_dir)`
  - Why: discovers available profile folders dynamically (for unknown/new profiles).

- `infer_geometry_from_spread_length(length_m)`
  - Why: maps measured spread length to geometry class (`100` or `200`) by nearest reference.

- `infer_geometry_from_seg2_file(seg2_path)`
  - Why: derives geometry candidate directly from SEG2 receiver locations when available.

## 3) Geometry workbook parsing (X/Y/Z midpoint)

- `_find_station_xyz_columns(df)`
  - Why: detects flexible header layouts for profile/station/X/Y/Z in geometry Excel files.

- `_to_float_or_none(v)`
  - Why: robust numeric conversion for mixed Excel values.

- `load_midpoint_xyz_from_geometry_excels(profile_name, station_mid, excel_paths)`
  - Why: fetches midpoint coordinates for summary export; repeated station rows use last occurrence.

## 4) Summary export

- `export_velocity_summary_excel(...)`
  - Why: appends/updates one profile row in consolidated summary workbook with line/length/XYZ/velocities/times/depths.
  - Includes fallback filename if workbook is locked by Excel.

## 5) Field-report readers

- `_read_excel_table(path, sheet_name=None)`
  - Why: unified Excel read with engine fallbacks (`xlrd`/`openpyxl`).

- `infer_geometry_from_field_report(path, profile_name, sheet_name=None)`
  - Why: optional geometry hint from field-report notes (e.g. short profile -> 100m).

- `load_profile_offsets_from_excel(...)`
  - Why: loads per-shot PO and optional inline shift defaults from field report blocks.

- `read_perpendicular_offsets(path)`
  - Why: reads simple text PO map file format.

## 6) Offset and PO helpers

- `true_offset(inline_m, perp_m)`
  - Why: computes geometric true offset using inline and perpendicular components.

- `parse_shot_value_map(text)`
  - Why: parses CLI maps like `1:0,2:3.5,3:0`.

- `resolve_perp_by_shot(cfg, shot_ids, override_map=None)`
  - Why: builds per-shot PO map from config plus CLI overrides.

## 7) Offset model interactive UI

- `prompt_offset_model_by_shot(...)`
  - Why: interactive PO/X-shift editor with live preview before analysis.
  - Preview now uses signed inline x for clearer left/right visualization.

- `build_corrected_pick_data(...)`
  - Why: transforms raw picks into corrected per-trace rows for analysis/export.

- `AnalysisWorkflow`
  - Why: orchestrates PO setup, corrected-data build, and optional layer-window picking.
  - `run()` is the main entry for post-picking analysis stage.

## 8) SEG2 and filtering

- `read_seg2(path)`
  - Why: loads SEG2 traces plus key metadata (FFID, delay, shot location, receiver locations).

- `_ormsby_response(...)`, `ormsby(...)`, `apply_ormsby_all(...)`, `apply_ormsby_all_params(...)`
  - Why: ProMAX-like Ormsby filtering functions.

- `butterworth_bandpass(...)`, `apply_butterworth_all_params(...)`
  - Why: alternative zero-phase Butterworth filtering.

- `apply_bulk_static(picks_raw)`
  - Why: applies static correction at analysis/export time (raw JSON remains unchanged).

- `apply_gain(...)`
  - Why: display gain modes (`none`, `norm`, `agc`) used in picker and preview.

## 9) Refraction model math

- `fit_line(x, y)`
  - Why: linear travel-time fitting, returns slope/intercept/r2.

- `depth_2layer(ti_ms, V1, V2)`
  - Why: 2-layer intercept-time depth estimate.

- `depth_3layer(ti2_ms, V1, V2, V3, h1)`
  - Why: second-interface depth for 3-layer model.

## 10) Layer window picking and review UI

- `_pick_layer_windows_from_plot(x_vals, t_vals, title)`
  - Why: interactive selection of layer x-windows.

- `_fit_layers_from_windows(x_vals, t_vals, windows)`
  - Why: builds 1-3 layer fits from selected windows and derives velocity/depth fields.

- `_compute_fit_rms(x_vals, t_vals, fit_res)`
  - Why: RMS quality metric for fitted arrivals.

- `_recompute_fit_derived(fit_res)`
  - Why: recalculates derived velocity/intercept/depth fields after manual fit edits.

- `_drop_fit_layer(fit_res, layer_idx)`
  - Why: allows rejecting one layer from review (e.g. unstable V1L).

- `_review_layer_fit_interactive(x_vals, t_vals, fit_res, title)`
  - Why: quality-control dialog for accept/repick/skip/drop-layer decisions.

- `pick_layer_windows_interactive(profile_name, corrected_by_shot, existing_results=None)`
  - Why: per-shot and per-side fit orchestration with user decisions persisted.

## 11) Aggregation and prediction helpers

- `compute_layer_averages(layer_results, shot_pos_by_id)`
  - Why: computes off/center/average velocity/time/depth summaries.

- `_predict_time_from_fit(x_abs, fit_res)`
  - Why: computes modeled arrival time for a given corrected offset from fitted segments.

- `_choose_shot_fit(layer_by_side)`
  - Why: picks representative side fit for shot-level summary where needed.

- `build_analysis_from_layers(corrected_by_shot, layer_results)`
  - Why: composes compact shot analysis payload for exports.

## 12) Export plots and files

- `export_fit_plot(...)`
  - Why: observed vs computed scatter summary with global RMS.

- `export_corrected_qc_plot(...)`
  - Why: corrected picks QC visualization (2 panels: channel and geometry-x views).

- `export_picks_txt(...)`
  - Why: writes clean picks text format compatible with downstream notebook flow.

- `_chdr(...)`, `_autofit_xl(...)`
  - Why: Excel styling helpers for all workbook outputs.

- `export_excel(...)`
  - Why: detailed per-profile workbook export (config, shot sheets, analysis, layer picks/averages).

- `export_tx_plot(...)`
  - Why: T-X summary with picks and fitted overlays.
  - Current behavior: overlays all available fitted sides (`L`, `R`, `ALL`) instead of one side only.

- `export_arrivals_observed_computed_plot(...)`
  - Why: T-X style overlay of observed and modeled arrivals in geometry-x domain, with channel axis and RMS.

## 13) Picker class (interactive first-break picking)

- `FirstBreakPicker`
  - Why: encapsulates full matplotlib GUI for picking one shot.

Important method groups:
- Build/layout:
  - `_build_figure`, `_build_controls`, `_make_mode_buttons`, `_style_mode_buttons`, `_refresh_mode_buttons`, `_refresh_toggle_buttons`.
- State toggles/sliders:
  - `_on_tmin`, `_on_tmax`, `_on_gain_mode`, `_on_agc_stat`, `_on_agc_window`, `_on_wiggle_stretch`, `_on_display_mode`, `_on_filter_mode`, `_on_butter_order`, `_on_filter_sliders`.
- Filter execution:
  - `_schedule_filter_update`, `_apply_filter_update_now`, `_recompute_filter`.
- Rendering:
  - `_active_data`, `_draw_traces`, `_draw_picks`, `_redraw`.
- Picking events:
  - `_nearest_idx`, `_on_click`, `_on_motion`, `_on_release`, `_range_fill`, `_auto_pick`.
- Navigation and persistence hooks:
  - `_save_and_finish`, `_go_prev`, `_auto_then_redraw`, `_on_key`, `_save_qc_image`, `run`.
- Data access:
  - `picks_ms` property.

## 14) Persistence helpers

Path builders:
- `_picks_json_path`, `_session_picks_json_path`, `_layer_json_path`, `_layer_session_json_path`.

JSON transforms and IO:
- `_coerce_layer_results`
- `load_layer_json`, `load_layer_session_json`
- `save_layer_json`, `save_layer_session_json`, `clear_layer_session_json`
- `load_picks_json`, `load_session_picks_json`
- `save_picks_json`, `save_session_picks_json`, `clear_session_picks_json`

Why: robust resume/finalize behavior for both picks and layer fits.

## 15) Top-level orchestration

- `process_profile(...)`
  - Why: full end-to-end run for one profile:
    - profile resolution
    - geometry choice (SEG2 length / field report / defaults)
    - picking session
    - PO defaults + UI correction
    - layer analysis
    - all exports

- `main()`
  - Why: CLI parsing, report discovery, target profile selection, batch execution.

## Runtime flow summary

1. `main()` resolves targets and report files.
2. `process_profile()` reads SEG2 and determines geometry.
3. Optional picker (`FirstBreakPicker`) produces raw picks.
4. `AnalysisWorkflow` applies PO/X-shift model and layer fit stage.
5. Export functions generate TXT/Excel/plots and summary workbook.

## Notes for upcoming PyQt migration

- Keep the analysis/persistence/export layers as-is.
- Replace only matplotlib UI shells first:
  - `prompt_offset_model_by_shot`
  - `_pick_layer_windows_from_plot`
  - `_review_layer_fit_interactive`
  - `FirstBreakPicker`
- Preserve existing payload contracts (`all_picks`, `corrected_by_shot`, `layer_results`) so exports continue to work unchanged.
