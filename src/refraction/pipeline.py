"""
Full-profile processing pipeline: pick -> correct -> analyze -> export.

This is the orchestration layer shared by both front-ends - `process_profile`
is what `lvl_refraction.py`'s CLI and `lvl_studio.py`'s "Run Full
Computation" both call, so there is exactly one implementation of "run
everything for this profile", not two.

Two pieces are deliberately NOT here, and are imported lazily (inside the
functions that need them) from `apps.lvl_refraction` instead:

- `FirstBreakPicker` (v1) - the interactive matplotlib picking window used
  when `pick_mode=True`. (Not to be confused with
  `src.picker.picker.FirstBreakPicker`, the new engine-v2 class - see the
  `FirstBreakPickerV2` import below.)
- `prompt_offset_model_by_shot` / `pick_layer_windows_interactive` - the
  interactive matplotlib PO-entry and layer-window-picking prompts used by
  `AnalysisWorkflow.run()`.

Those are genuinely interactive, desktop-only UI and don't belong in a
framework-agnostic module - a future web backend would replace exactly
those two call sites with its own UI flow, while reusing everything else
here unchanged. The lazy import is what keeps that possible without a
circular import (this module is imported BY `apps/lvl_refraction.py`, so
it can't import from it at module load time).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from src.common import paths as core_paths
from src.common.settings import PROFILES, USE_SEG2_SHOT_POSITION, BULK_SHIFT_MS, THEME
from src.picker.preprocessing import apply_ormsby_all
from src.picker.picker import FirstBreakPicker as FirstBreakPickerV2
from src.picker.settings import PickerSettings as PickerSettingsV2
from src.io.seg2_reader import read_seg2, read_seg2_acquisition_time_de, read_seg2_mid_xyz
from src.io.pick_reader import (
    load_layer_json, load_layer_session_json, save_layer_json,
    save_layer_session_json, clear_layer_session_json,
    load_picks_json, load_session_picks_json, save_picks_json,
    save_session_picks_json, clear_session_picks_json,
)
from src.utils.geometry import (
    load_geometry, auto_shot_positions,
    infer_geometry_from_spread_length, infer_geometry_from_seg2_file,
    load_profile_geometry_from_excels, infer_geometry_from_field_report,
    load_profile_offsets_from_excel,
)
from src.refraction.processing import resolve_perp_by_shot
from src.refraction.layer_analysis import build_corrected_pick_data, build_analysis_from_layers
from src.io.exporters import (
    export_velocity_summary_excel, export_excel, export_picks_txt,
    export_tx_plot, export_arrivals_observed_computed_plot, export_layer_fit_rms_plot,
)


class AnalysisWorkflow:
    """
    Dedicated post-picking analysis workflow.

    Responsibilities
    ----------------
    - Collect/edit per-shot PO values in a UI.
    - Compute corrected offsets per trace using:
        XO_k = sqrt((x_k - SP)^2 + PO^2)
      where x_k comes from geometry, SP is shot position, and PO is
      perpendicular shot-to-line offset.
    - Show layer-window picking plots and return fitted layer metrics.
    """

    def __init__(self, profile_name: str, cfg: dict,
                 shots_info: list, shot_label_pos: dict,
                 all_picks: dict, recv_positions: Any,
                 perp_override: dict | None = None,
                 inline_shift_override: dict | None = None,
                 enable_layer_pick: bool = True,
                 existing_layer_results: dict | None = None):
        self.profile_name = profile_name
        self.cfg = cfg
        self.shots_info = shots_info
        self.shot_label_pos = shot_label_pos
        self.all_picks = all_picks
        self.recv_positions = recv_positions
        self.perp_override = perp_override or {}
        self.inline_shift_override = inline_shift_override or {}
        self.enable_layer_pick = bool(enable_layer_pick)
        self.existing_layer_results = existing_layer_results or {}

        self.perp_by_shot: dict = {}
        self.inline_shift_by_shot: dict = {}
        self.corrected_by_shot: dict = {}
        self.layer_results: dict = {}

    def run(self) -> dict:
        # Lazy import: interactive matplotlib prompts, desktop-CLI only.
        from lvl_refraction import prompt_offset_model_by_shot, pick_layer_windows_interactive

        shot_ids = [sid for sid, _ in self.shots_info]
        self.perp_by_shot = resolve_perp_by_shot(self.cfg, shot_ids, self.perp_override)
        self.inline_shift_by_shot = {int(sid): 0.0 for sid in shot_ids}
        for sid, val in (self.inline_shift_override or {}).items():
            self.inline_shift_by_shot[int(sid)] = float(val)

        print("  Step 1/2: Perpendicular offset (PO) setup")
        self.perp_by_shot, self.inline_shift_by_shot = prompt_offset_model_by_shot(
            self.perp_by_shot,
            inline_shift_by_shot=self.inline_shift_by_shot,
            shots_info=self.shots_info,
            all_picks=self.all_picks,
            recv_positions=self.recv_positions,
        )
        print("  Perp offsets by shot: "
              + ", ".join(f"S{sid}={self.perp_by_shot.get(sid, 0.0):.2f}m" for sid in shot_ids))
        print("  Inline shift by shot: "
              + ", ".join(f"S{sid}={self.inline_shift_by_shot.get(sid, 0.0):.2f}m" for sid in shot_ids))

        print("  Step 2/2: Corrected picks and layer windows")
        self.corrected_by_shot = build_corrected_pick_data(
            self.shots_info, self.all_picks, self.recv_positions,
            perp_by_shot=self.perp_by_shot,
            inline_shift_by_shot=self.inline_shift_by_shot,
        )

        if self.enable_layer_pick:
            self.layer_results = pick_layer_windows_interactive(
                self.profile_name, self.corrected_by_shot,
                existing_results=self.existing_layer_results,
            )
        else:
            self.layer_results = self.existing_layer_results

        return {
            "perp_by_shot": self.perp_by_shot,
            "inline_shift_by_shot": self.inline_shift_by_shot,
            "corrected_by_shot": self.corrected_by_shot,
            "layer_results": self.layer_results,
        }


def process_profile(profile_name: str, pick_mode: bool = True,
                    geom_override: int | None = None,
                    perp_by_shot_override: dict | None = None,
                    inline_shift_by_shot_override: dict | None = None,
                    perp_excel_cfg: dict | None = None,
                    enable_layer_pick: bool = True,
                    control_file: Path | None = None,
                    show_plot_controls: bool = True,
                    raw_dir: Path | None = None,
                    auto_pick: bool = False,
                    auto_pick_settings: "PickerSettingsV2 | None" = None):
    """
    Run the full profile pipeline: pick, correct, analyze, and export.

    Parameters
    ----------
    profile_name : str
        Profile key from PROFILES (also expected folder name under the
        project's raw data folder).
    pick_mode : bool, optional
        True: open interactive picker. False: export-only mode using saved picks.
    geom_override : int | None, optional
        Explicit geometry type (100 or 200) overriding auto/default selection.
    perp_by_shot_override : dict | None, optional
        CLI-provided perpendicular offsets by shot id.
    inline_shift_by_shot_override : dict | None, optional
        CLI-provided inline shifts (m) by shot id.
    perp_excel_cfg : dict | None, optional
        Excel mapping configuration for loading default PO/X-shift values.
    enable_layer_pick : bool, optional
        Enable interactive layer-window picking/review stage.
    raw_dir : Path | None, optional
        Folder containing one SEG2 subfolder per profile. Defaults to the
        active project's `raw_folder` if not given explicitly.
    auto_pick : bool, optional
        If True, skip the interactive picker entirely and pick every shot
        automatically with the Picker Engine V2 (`src.picker.picker.
        FirstBreakPicker`) - preprocessing -> features -> likelihood ->
        coherence -> path optimization -> confidence/quality, per shot
        gather. Results feed straight into the same layer-analysis and
        export stages as manual picks, so you get full velocity/depth
        output without picking anything by hand. Low-confidence traces
        are reported, not silently dropped - review `picks.session.json`
        / the printed summary before trusting a profile you haven't
        checked yet.
    auto_pick_settings : PickerSettingsV2 | None, optional
        Tuning for the v2 engine (feature weights, coherence radius,
        smoothness/jump penalties, ...). Defaults to `PickerSettingsV2()`.
    """
    from lvl_refraction import FirstBreakPicker  # v1 interactive picker (desktop-CLI only)
    cfg = PROFILES.get(profile_name)
    if cfg is None:
        cfg = {
            "geom": 200,
            "line_no": profile_name,
            "perp_m": 0.0,
            "shots": "auto",
        }
        print(f"  [INFO] Profile '{profile_name}' not in PROFILES; using dynamic defaults.")

    if raw_dir is None:
        print(
            "[ERROR] No raw_dir given. Pass raw_dir=project.raw_folder "
            "(see src.io.project_io.Project) or call main() with --project."
        )
        return
    data_dir = Path(raw_dir) / profile_name
    if not data_dir.exists():
        print(f"[ERROR] Data folder not found: {data_dir}")
        return

    def _ffid(p: Path) -> int:
        digits = "".join(c for c in p.stem if c.isdigit())
        return int(digits) if digits else 0

    seg2_candidates = list(data_dir.glob("*.seg2")) + list(data_dir.glob("*.SEG2"))
    seg2_unique = {str(p.resolve()).lower(): p for p in seg2_candidates}
    seg2_files = sorted(seg2_unique.values(), key=_ffid)
    if not seg2_files:
        print(f"[ERROR] No .seg2 files found in {data_dir}")
        return
    print(f"  Found {len(seg2_files)} SEG2 file(s)")

    inferred_geom: int | None = None
    report_paths = [Path(p) for p in (perp_excel_cfg or {}).get("paths", []) if p]
    if geom_override is None and report_paths:
        for rp in report_paths:
            try:
                inferred = infer_geometry_from_field_report(
                    path=rp,
                    profile_name=profile_name,
                    sheet_name=(perp_excel_cfg or {}).get("sheet"),
                )
            except Exception as exc:
                print(f"  [WARN] Geometry inference failed in {rp.name}: {exc}")
                continue
            if inferred in (100, 200):
                inferred_geom = int(inferred)
                print(f"  Geometry hint from field report: {rp.name} -> {inferred_geom} m")
                break

    geom_from_seg2, seg2_len_m = infer_geometry_from_seg2_file(seg2_files[0])

    geom_from_excel = None
    excel_len_m = None
    geometry_paths = [Path(p) for p in (perp_excel_cfg or {}).get("geometry_paths", []) if p]
    if geometry_paths:
        g_len, _gx, _gy, _gz, _gsrc, _smin, _smid, _smax, _acq_dt_de = load_profile_geometry_from_excels(
            profile_name=profile_name,
            excel_paths=geometry_paths,
        )
        if g_len is not None and g_len > 0.0:
            excel_len_m = float(g_len)
            geom_from_excel = infer_geometry_from_spread_length(excel_len_m)

    if geom_override is not None:
        geom_type = int(geom_override)
        geom_src = "CLI"
    elif geom_from_excel in (100, 200):
        geom_type = int(geom_from_excel)
        geom_src = f"LVL geometry excel ({excel_len_m:.2f} m)"
    elif int(cfg.get("geom", 200) or 200) in (100, 200):
        geom_type = int(cfg.get("geom", 200) or 200)
        geom_src = "profile config"
    elif geom_from_seg2 in (100, 200):
        geom_type = int(geom_from_seg2)
        geom_src = f"SEG2 spread-length ({seg2_len_m:.2f} m)"
    elif inferred_geom in (100, 200):
        geom_type = int(inferred_geom)
        geom_src = "field-report"
    else:
        geom_type = int(cfg.get("geom", 200) or 200)
        if geom_type not in (100, 200):
            geom_type = 200
        geom_src = "default"

    perp_cfg       = cfg.get("perp_m", 0.0)
    perp_m         = float(perp_cfg if not isinstance(perp_cfg, dict) else 0.0)
    recv_positions = load_geometry(geom_type)
    print(f"  Geometry selected: {geom_type} m ({geom_src})")

    # Resolve shot positions: "auto" derives them from geometry
    shots_cfg = cfg.get("shots", "auto")
    if shots_cfg == "auto" or not isinstance(shots_cfg, dict):
        shots_cfg = auto_shot_positions(recv_positions)
        print(f"  Shot positions (auto from geometry): "
              + "  ".join(f"Shot{k}={v:.3f}m" for k, v in shots_cfg.items()))

        print(f"  Geometry {geom_type} m : {len(recv_positions)} receivers, "
            f"{recv_positions[0]:.2f} - {recv_positions[-1]:.2f} m  "
            f"|  perp(default) = {perp_m:.1f} m  "
            f"|  bulk static = {BULK_SHIFT_MS:+.1f} ms")

    final_picks: dict = load_picks_json(profile_name)
    session_picks: dict = load_session_picks_json(profile_name)
    if session_picks:
        all_picks: dict = session_picks
        print("  Session picks found: resuming from picks.session.json")
    else:
        all_picks = {sid: dict(vals) for sid, vals in final_picks.items()}
    shots_meta: list = []
    qc_dir = core_paths.require_active_project().plots_dir_for(profile_name)
    finalized = False

    shot_cache: list = []
    for file_idx, seg2_path in enumerate(seg2_files):
        shot_id = file_idx + 1
        data_raw, dt_s, n_tr, n_samp, shot_pos_hdr, ffid_hdr, \
            delay_ms, recv_locs_hdr = read_seg2(seg2_path)
        n_show = min(n_tr, len(recv_positions))
        shot_pos_cfg = shots_cfg.get(shot_id)
        shot_pos_source = "default"
        shot_pos_m = 0.0

        if shot_pos_cfg is not None:
            shot_pos_m = float(shot_pos_cfg)
            shot_pos_source = "config"
        elif shot_pos_hdr is not None:
            shot_pos_m = float(shot_pos_hdr)
            shot_pos_source = "header"

        if USE_SEG2_SHOT_POSITION and shot_pos_hdr is not None:
            shot_pos_m = float(shot_pos_hdr)
            shot_pos_source = "header"
        shot_pos_nominal = (float(shot_pos_cfg)
                            if shot_pos_cfg is not None else float(shot_pos_m))
        shot_cache.append({
            "shot_id": shot_id,
            "seg2_path": seg2_path,
            "data_raw": data_raw,
            "dt_s": dt_s,
            "n_tr": n_tr,
            "n_samp": n_samp,
            "shot_pos_hdr": shot_pos_hdr,
            "ffid_hdr": ffid_hdr,
            "delay_ms": delay_ms,
            "n_show": n_show,
            "shot_pos_m": shot_pos_m,
            "shot_pos_source": shot_pos_source,
            "shot_pos_nominal": shot_pos_nominal,
        })

    idx = 0
    while 0 <= idx < len(shot_cache):
        shot = shot_cache[idx]
        shot_id = int(shot["shot_id"])
        seg2_path = shot["seg2_path"]
        print(f"\n  -- Shot {shot_id}  ({seg2_path.name}) --")

        data_raw = shot["data_raw"]
        dt_s = float(shot["dt_s"])
        n_tr = int(shot["n_tr"])
        n_samp = int(shot["n_samp"])
        shot_pos_hdr = shot["shot_pos_hdr"]
        shot_pos_source = shot.get("shot_pos_source", "default")
        ffid_hdr = shot["ffid_hdr"]
        delay_ms = float(shot["delay_ms"])
        n_show = int(shot["n_show"])
        shot_pos_m = float(shot["shot_pos_m"])

        dt_ms  = dt_s * 1000.0
        t_end  = delay_ms + (n_samp - 1) * dt_ms
        print(f"     FFID={ffid_hdr}  |  {n_tr} traces  |  "
              f"dt={dt_ms:.4f} ms  |  delay={delay_ms:.1f} ms  |  "
              f"{n_samp} smp  ({delay_ms:.1f} to {t_end:.1f} ms)")

        if shot_pos_source == "header":
            print(f"     Shot pos : {shot_pos_m:.2f} m  (SEG2 header)")
        elif shot_pos_source == "config":
            print(f"     Shot pos : {shot_pos_m:.2f} m  (CONFIG)")
        elif shot_pos_hdr is not None:
            print(f"     Shot pos : {shot_pos_m:.2f} m  (SEG2 header; out-of-range fallback)")
        else:
            print("     Shot pos : 0.0 m  [WARN: defaulting to 0]")

        if auto_pick:
            # Automatic picking with Picker Engine V2 - no interactive
            # window, no per-shot navigation; picks every trace of this
            # shot gather using cross-trace coherence + path optimization,
            # then moves straight to the next shot.
            geom_slice = recv_positions[:n_show]
            data_slice = data_raw[:n_show]
            data_filt  = apply_ormsby_all(data_slice, dt_s)

            offsets_m = [float(x) - float(shot_pos_m) for x in geom_slice]
            v2_settings = auto_pick_settings or PickerSettingsV2()
            v2_picker = FirstBreakPickerV2(settings=v2_settings)
            try:
                v2_results = v2_picker.pick_profile(data_filt, dt_s, offsets_m=offsets_m)
            except Exception as exc:
                print(f"     [WARN] Auto-pick (v2) failed for shot {shot_id}: {exc}")
                v2_results = []

            picks_for_shot: dict = {}
            low_conf: list = []
            for trace_idx, r in enumerate(v2_results):
                if r.sample is None:
                    continue
                t_abs = round(float(delay_ms) + float(r.sample) * dt_s * 1000.0, 2)
                picks_for_shot[trace_idx + 1] = t_abs
                if r.confidence is not None and r.confidence < v2_settings.minimum_confidence:
                    low_conf.append(trace_idx + 1)

            all_picks[shot_id] = picks_for_shot
            save_session_picks_json(profile_name, all_picks)
            conf_note = f", {len(low_conf)} flagged low-confidence (traces {low_conf})" if low_conf else ""
            print(f"     Shot {shot_id}: auto-picked {len(picks_for_shot)}/{len(v2_results)} "
                  f"trace(s) with Picker V2{conf_note}.")

        elif pick_mode:
            # Always use geometry file positions for the x-axis display and picking.
            # SEG2 RECEIVER_LOCATION may differ from the geometry file (e.g. a
            # 200-m spread recorded with 0-94 m header values).
            geom_slice = recv_positions[:n_show]

            data_slice = data_raw[:n_show]
            data_filt  = apply_ormsby_all(data_slice, dt_s)

            def _save_cb(picks_for_shot: dict, _sid: int = shot_id) -> None:
                """Called on nav/finalize to persist temporary session progress."""
                all_picks[_sid] = picks_for_shot
                save_session_picks_json(profile_name, all_picks)
                print(f"     âœ“ {len(picks_for_shot)} pick(s) saved to session.")

            picker = FirstBreakPicker(
                data_slice, data_filt, dt_s,
                geom_slice, shot_id, profile_name,
                shot_pos_m=shot_pos_m,
                delay_ms=delay_ms,
                existing_picks=all_picks.get(shot_id, {}),
                qc_dir=qc_dir,
                save_callback=_save_cb,
                header_info={
                    "ffid": ffid_hdr,
                    "n_tr": n_tr,
                    "n_samp": n_samp,
                    "shot_pos_hdr": shot_pos_hdr,
                },
                control_file=control_file,
                show_plot_controls=show_plot_controls,
            )
            result = picker.run() or {"status": "quit", "picks": all_picks.get(shot_id, {})}
            status = result.get("status", "next")
            all_picks[shot_id] = dict(result.get("picks", {}))
            save_session_picks_json(profile_name, all_picks)

            if all_picks.get(shot_id):
                print(f"     Shot {shot_id}: {len(all_picks[shot_id])} pick(s) in memory.")
            else:
                print(f"     Picking cancelled / no picks for shot {shot_id}.")

            if status == "prev":
                if idx > 0:
                    idx -= 1
                else:
                    print("     Already at first shot; staying on current shot.")
                continue
            if status == "next":
                if idx < (len(shot_cache) - 1):
                    idx += 1
                else:
                    print("     Already at last shot; staying on current shot.")
                continue
            if status == "finalize":
                finalized = True
                break
            if status == "quit":
                break

        idx += 1

    if auto_pick:
        finalized = True

    acquisition_time_de = None
    seg2_mid_xyz = (None, None, None)
    if shot_cache:
        mid_shot = shot_cache[len(shot_cache) // 2]
        mid_seg2 = Path(mid_shot.get("seg2_path"))
        acquisition_time_de = read_seg2_acquisition_time_de(mid_seg2)
        seg2_mid_xyz = read_seg2_mid_xyz(mid_seg2)
        if acquisition_time_de:
            print(f"  Acquisition time (middle shot, DE): {acquisition_time_de}")

    preview_only = bool(pick_mode and not finalized)
    if preview_only:
        print("\n  Picking session ended without Save/Close finalization.")
        print("  Continuing to correction UI in PREVIEW mode.")
        print("  Final picks.json remains unchanged until Save/Close is used.")

    for shot in shot_cache:
        shot_id = int(shot["shot_id"])
        shots_meta.append((shot_id, float(shot["shot_pos_m"]), float(shot["shot_pos_nominal"])))

    excel_perp_by_shot: dict = {}
    excel_shift_by_shot: dict = {}
    ffid_by_shot = {
        int(s["shot_id"]): int(s.get("ffid_hdr", 0) or 0)
        for s in shot_cache
    }
    for rp in report_paths:
        try:
            po_map, sh_map = load_profile_offsets_from_excel(
                path=rp,
                profile_name=profile_name,
                ffid_by_shot=ffid_by_shot,
                sheet_name=(perp_excel_cfg or {}).get("sheet"),
                ffid_col=str((perp_excel_cfg or {}).get("ffid_col", "A")),
                perp_col=str((perp_excel_cfg or {}).get("perp_col", "D")),
                profile_col=str((perp_excel_cfg or {}).get("profile_col", "F")),
                inline_shift_col=(perp_excel_cfg or {}).get("inline_shift_col"),
            )
        except Exception as exc:
            print(f"  [WARN] Could not parse PO defaults from {rp.name}: {exc}")
            continue
        if po_map:
            excel_perp_by_shot.update(po_map)
            print("  Excel defaults loaded for PO from "
                  f"{rp.name}: "
                  + ", ".join(f"S{k}={v:.2f}" for k, v in sorted(po_map.items())))
        if sh_map:
            excel_shift_by_shot.update(sh_map)
            print("  Excel defaults loaded for X-shift from "
                  f"{rp.name}: "
                  + ", ".join(f"S{k}={v:.2f}" for k, v in sorted(sh_map.items())))

    perp_effective = dict(excel_perp_by_shot)
    perp_effective.update(perp_by_shot_override or {})
    shift_effective = dict(excel_shift_by_shot)
    shift_effective.update(inline_shift_by_shot_override or {})

    if not perp_effective:
        print("\n  [WARN] No perpendicular offsets found in field reports.")
        print("         Opening geometry editor with default PO=0.0 for all shots.")
    for s in shot_cache:
        sid = int(s["shot_id"])
        perp_effective.setdefault(sid, 0.0)
        shift_effective.setdefault(sid, 0.0)

    total_picks = sum(len(v) for v in all_picks.values())
    if total_picks == 0:
        print("\n  No picks to export.")
        return

    analysis: dict = {}

    layer_final: dict = load_layer_json(profile_name)
    layer_session: dict = load_layer_session_json(profile_name)
    existing_layer_results: dict = layer_session if layer_session else layer_final
    if layer_session:
        print("  Layer session found: resuming from layer_analysis.session.json")

    # Exports / correction UI
    print(f"\n  -- Exporting / correction  ({total_picks} picks) --")
    shots_info_proc = [(sid, sp_proc) for sid, sp_proc, _ in shots_meta]
    shot_label_pos  = {sid: sp_nom for sid, _, sp_nom in shots_meta}
    analysis_ui = AnalysisWorkflow(
        profile_name=profile_name,
        cfg=cfg,
        shots_info=shots_info_proc,
        shot_label_pos=shot_label_pos,
        all_picks=all_picks,
        recv_positions=recv_positions,
        perp_override=perp_effective,
        inline_shift_override=shift_effective,
        enable_layer_pick=enable_layer_pick,
        existing_layer_results=existing_layer_results,
    )
    analysis_bundle = analysis_ui.run()
    perp_by_shot = analysis_bundle["perp_by_shot"]
    corrected_by_shot = analysis_bundle["corrected_by_shot"]
    layer_results = analysis_bundle["layer_results"]
    analysis = build_analysis_from_layers(corrected_by_shot, layer_results)
    save_layer_session_json(profile_name, layer_results)

    qc_suffix = "_preview" if preview_only else ""

    if preview_only:
        export_excel(profile_name, shots_info_proc, all_picks, recv_positions,
                     analysis, cfg,
                     corrected_by_shot=corrected_by_shot,
                     layer_results=layer_results,
                     perp_by_shot=perp_by_shot,
                     filename_suffix="_preview")
        export_tx_plot(profile_name, shots_info_proc, all_picks, recv_positions,
                       analysis, perp_m=perp_m, shot_label_pos=shot_label_pos,
                       corrected_by_shot=corrected_by_shot, layer_results=layer_results,
                       theme=THEME)
        export_arrivals_observed_computed_plot(
            profile_name, corrected_by_shot, layer_results, filename_suffix="_preview", theme=THEME
        )
        export_layer_fit_rms_plot(
            profile_name, corrected_by_shot, layer_results, filename_suffix="_preview", theme=THEME
        )
        print("  Preview files written; final picks.json not updated.")
        return

    export_picks_txt(profile_name, shots_info_proc, all_picks, recv_positions)
    export_excel(profile_name, shots_info_proc, all_picks, recv_positions,
                 analysis, cfg,
                 corrected_by_shot=corrected_by_shot,
                 layer_results=layer_results,
                 perp_by_shot=perp_by_shot)
    export_tx_plot(profile_name, shots_info_proc, all_picks, recv_positions,
                   analysis, perp_m=perp_m, shot_label_pos=shot_label_pos,
                   corrected_by_shot=corrected_by_shot, layer_results=layer_results,
                   theme=THEME)
    export_arrivals_observed_computed_plot(profile_name, corrected_by_shot, layer_results, theme=THEME)
    export_layer_fit_rms_plot(profile_name, corrected_by_shot, layer_results, theme=THEME)
    geometry_excels = [Path(p) for p in (perp_excel_cfg or {}).get("geometry_paths", []) if p]
    export_velocity_summary_excel(
        profile_name=profile_name,
        cfg=cfg,
        recv_positions=recv_positions,
        shots_info=shots_info_proc,
        layer_results=layer_results,
        analysis=analysis,
        output_dir=core_paths.require_active_project().velocity_dir,
        geometry_excel_paths=geometry_excels,
        acquisition_time_de=acquisition_time_de,
        seg2_mid_xyz=seg2_mid_xyz,
    )
    save_picks_json(profile_name, all_picks)
    save_layer_json(profile_name, layer_results)
    clear_session_picks_json(profile_name)
    clear_layer_session_json(profile_name)
    _ap = core_paths.require_active_project()
    print(f"\n  Output -> " f"{(_ap.results_dir / profile_name).relative_to(_ap.root)}")


__all__ = ["AnalysisWorkflow", "process_profile"]