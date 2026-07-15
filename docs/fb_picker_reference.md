# FB Picker Reference

This document explains:
- what the main functions in fb_picker.py do,
- the physics assumptions behind each picking/filtering method,
- which parameters are worth tuning in practice,
- a practical workflow to test and choose settings.

Scope: first-break picking only (no inversion/evaluation workflow).

## 1) High-level workflow

The picker pipeline is:
1. Read SEG2 shot gathers and basic header metadata.
2. Assign geometry and shot positions.
3. Build filtered versions of the gather (Ormsby, Butterworth bandpass, Cutpass).
4. Display traces with optional gain and polarity.
5. Place picks manually or by auto-picker.
6. Save session picks and finalize to picks.json.

Two important design points:
- Stored picks are in absolute time in ms (including SEG2 delay and optional trigger correction).
- Gain is display/picking conditioning, not a physical processing result to be exported as data.

## 2) Function map by module section

## Geometry and input

- load_geometry(geom_type)
  - Reads receiver positions from geometry100.txt or geometry200.txt.
  - Returns receiver x-positions in meters.

- auto_shot_positions(recv_pos)
  - Computes three standard shot positions from spread edges and midpoint.

- discover_profile_folders(data_dir)
  - Finds profile directories containing SEG2 files.

- discover_lvl_geometry_excels(data_dir)
  - Finds LVL Excel files likely containing station/coordinate metadata.

- infer_geometry_from_lvl_excel(data_dir, profile_name)
  - Detects whether profile likely matches 100 m or 200 m geometry from station spacing / XY distance.

- read_seg2(path)
  - Loads one SEG2 file and returns gather array plus metadata:
    - dt, trace/sample count, shot position (if present), FFID, delay, receiver locations.

## Signal processing filters and helpers

- ormsby(...)
  - Frequency-domain trapezoidal band-pass (f1-f2-f3-f4).
  - Smooth transition at low and high edges.

- butterworth_bandpass(...)
  - Zero-phase IIR band-pass using scipy SOS and filtfilt.

- butterworth_highpass(...), butterworth_lowpass(...)
  - Zero-phase high-pass and low-pass stages.

- apply_cutpass_all_params(...)
  - Sequential high-pass at f1, then low-pass at f4.
  - This is the cutpass mode.

- _zero_crossing_from_extremum_samples(...)
  - Finds zero-crossing near a local extremum in a time window.
  - Used for onset-like pick near dominant wavelet excursion.

- hilbert_envelope_pick(...)
  - Builds analytic signal envelope with Hilbert transform.
  - Returns envelope peak time and threshold-onset time.

- apply_gain(...)
  - Modes: none, norm, agc.
  - AGC uses running mean-abs or RMS window.

## UI and picker core

- FirstBreakPicker._active_data()
  - Selects current filter mode and applies gain/polarity for display and auto logic.

- FirstBreakPicker._recompute_filter()
  - Rebuilds Ormsby, Butterworth, and Cutpass cached datasets after slider/order changes.

- FirstBreakPicker._apply_filter_single()
  - Applies active filter to one trace when preparing pick-time snapping.

- FirstBreakPicker._prepare_trace_for_pick()
  - Applies process order:
    - F>G>X, G>F>X, or RAW>X
  - Then optional polarity inversion.

- FirstBreakPicker._snap_pick_time_ms(...)
  - Manual pick refinement around click hint:
    - stalta mode: keep hint (plus timing correction),
    - maxdiff_zero: peak then zero-cross,
    - hilbert_env: onset or peak from envelope.

- FirstBreakPicker._auto_pick_stalta()
  - Classical STA/LTA trigger in absolute-time search zone.

- FirstBreakPicker._auto_pick_maxdiff_zero()
  - Auto pick from extremum and nearby zero-crossing within gate.

- FirstBreakPicker._auto_pick_hilbert_env()
  - Auto pick from Hilbert envelope onset/peak within gate.

- FirstBreakPicker._auto_abs_gate_ms_for_trace(...)
  - Velocity-based plausibility window using shot-receiver offset.

## Persistence and profile processing

- load_picks_json, save_picks_json
  - Final persistent picks.

- load_session_picks_json, save_session_picks_json
  - In-progress session recovery.

- process_profile(profile_name, geom_override)
  - Runs the per-shot picker loop and manages navigation/finalization.

## 3) Physics behind methods

## 3.1 First-break concept

The first break is the earliest detectable arrival from source to receiver.
In near-surface seismic refraction, it is controlled by:
- source coupling and trigger timing,
- direct wave path at near offsets,
- critically refracted arrivals at larger offsets,
- noise (cultural, instrument, coupling, ground roll).

The picker therefore balances two goals:
- detect earliest physically plausible onset,
- stay robust against noise and wavelet variability.

## 3.2 Filter physics

- Ormsby (f1-f2-f3-f4)
  - Piecewise linear amplitude spectrum.
  - Good when you want explicit corner control and stable passband shape.

- Butterworth bandpass (f2-f3, order)
  - IIR with smooth monotonic response.
  - Higher order gives steeper roll-off but can sharpen waveform behavior.

- Cutpass (f1 low-cut, f4 high-cut, order)
  - Two-stage Butterworth:
    1) high-pass at f1 removes very low frequency drift/ground roll,
    2) low-pass at f4 removes high-frequency noise.
  - Practical interpretation: keep the useful middle band while trimming both tails.

All filter implementations here are zero-phase (forward/backward), so they minimize phase shift in pick timing.

## 3.3 Picking physics

- STA/LTA
  - Detects sudden energy increase relative to background.
  - Works best with stable noise statistics and clear onset rise.

- Maxdiff plus zero-cross
  - Uses dominant excursion then finds local zero crossing.
  - Often approximates a wavelet onset proxy for impulsive arrivals.

- Hilbert envelope
  - Envelope tracks instantaneous amplitude of analytic signal.
  - Peak target: robust for high SNR but can be late.
  - Onset target: earlier arrival estimate via threshold crossing above noise floor.

## 3.4 Velocity gate

Arrival time should roughly satisfy:
- t approximately offset / velocity

The gate uses Vmax and Vmin with padding to reject non-physical picks.
This is especially useful for far-offset false picks and noisy windows.

## 4) Parameter guide and tuning ranges

Use these as starting ranges, then calibrate per profile.

## Core filter parameters

- f1 (Hz): low-cut corner
  - Typical test range: 1 to 8
  - Raise if low-frequency ground roll dominates.

- f2 (Hz), f3 (Hz): Butter passband
  - Typical f2: 3 to 15
  - Typical f3: 80 to 180

- f4 (Hz): high-cut corner
  - Typical test range: 120 to 240
  - Lower if high-frequency noise causes jitter.

- butter order
  - Typical: 2, 4, 6
  - Start at 4.
  - Increase only if transitions too soft.

Recommended starting presets:
- Quiet data: butter, f2=5, f3=150, order=4
- Strong low-frequency noise: cutpass, f1=4, f4=170, order=4
- Broad noisy spectrum: ormsby, 2-5-130-180

## Gain and display conditioning

- gain mode
  - norm: stable general default.
  - agc: helpful when amplitudes vary strongly by trace.

- AGC window (ms)
  - Typical: 80 to 250
  - Too short: over-equalizes and can distort onset prominence.
  - Too long: weak correction.

- AGC stat
  - rms: usually smoother for seismic waveforms.
  - mean: sometimes better for impulsive/non-Gaussian noise.

## Picker method parameters

- AUTO_PICK_MODE
  - hilbert_env: robust default in many datasets.
  - maxdiff_zero: good onset alignment on impulsive wavelets.
  - stalta: simple trigger baseline.

- HILBERT_TARGET
  - onset for first-arrival timing.
  - peak for stronger but often later picks.

- HILBERT_ONSET_PCT
  - Typical: 0.08 to 0.20
  - Lower values pick earlier and can be noise-sensitive.

- MANUAL_SNAP_WIN_MS
  - Typical: 4 to 10
  - Larger allows stronger auto correction around manual click.

- ZERO_X_SEARCH_DIRECTION
  - backward usually better for onset before dominant peak.

## Plausibility and timing correction

- AUTO_USE_VELOCITY_GATE
  - Keep enabled for production autopicking.

- AUTO_VMIN_M_S, AUTO_VMAX_M_S
  - Set from near-surface expectation for your site.
  - Example broad start: 120 to 3500.

- AUTO_GATE_PAD_MS
  - Typical: 10 to 30
  - Increase if gate clips true picks at transitions.

- TRIGGER_STATIC_MS
  - Use to correct fixed radio/trigger delay.
  - Calibrate from known offsets or repeated shots.

- TRIGGER_MS_PER_M
  - Use only if delay changes with offset/cable geometry.

## 5) Practical tuning protocol

1. Start with mode hilbert_env, target onset, filter butter.
2. Set a reasonable band (f2/f3), order 4.
3. Pick a few near, mid, far traces manually.
4. Run auto picker and inspect residual pattern vs offset.
5. If picks are late and smooth, lower onset threshold slightly.
6. If picks are noisy/jittery, lower f4 or raise f1.
7. If far offsets are wrong, tighten velocity gate and pad.
8. Apply TRIGGER_STATIC_MS only after waveform behavior is stable.

## 6) Common failure patterns and fixes

- Picks on wavelet tail:
  - Use hilbert onset, lower onset_pct slightly, or use maxdiff_zero backward.

- Too few auto picks:
  - Widen velocity gate, increase pad, relax filter aggressiveness.

- Picks mostly near offsets only:
  - Check Vmin too high, window too narrow, or far-offset SNR too low.

- Strong line-to-line variability:
  - Prefer per-profile parameter presets and keep session picks.

## 7) Notes on interpretation

- First-break time precision is limited by dt, noise, and trigger uncertainty.
- Zero-phase filtering preserves timing better than causal filtering for this use.
- Envelope peak is not the same as physical first motion; onset is generally closer to first break.

## 8) Suggested future extension

If you want, a small config block with named presets (for example soft_ground, urban_noise, clean_refraction) can be added to switch parameter sets quickly per profile.

## 9) Practical examples (copy, run, adjust)

These examples are designed so a new user can run the picker quickly and tune only what matters.

Important:
- The CLI call selects profile and optional geometry.
- Most tuning values are currently set in fb_picker.py config constants.
- After changing constants, rerun the script.

### Example A: Balanced default for most lines

Use when data quality is moderate and arrivals are visible on most traces.

Suggested constants in fb_picker.py:
- DEFAULT_FILTER_MODE = "butter"
- BP_F2 = 5.0
- BP_F3 = 150.0
- BUTTER_ORDER = 4
- AUTO_PICK_MODE = "hilbert_env"
- HILBERT_TARGET = "onset"
- HILBERT_ONSET_PCT = 0.10
- GAIN_MODE = "norm"
- AUTO_USE_VELOCITY_GATE = True
- AUTO_VMIN_M_S = 120.0
- AUTO_VMAX_M_S = 3500.0
- AUTO_GATE_PAD_MS = 20.0

Run:

```powershell
python fb_picker.py 110
```

What to check:
1. Near, mid, and far offsets all receive picks.
2. Picks sit near first onset, not on late wavelet tail.
3. If picks are too late, reduce HILBERT_ONSET_PCT to 0.08.

### Example B: Strong low-frequency noise (ground roll)

Use when low-frequency energy masks early arrivals.

Suggested constants:
- DEFAULT_FILTER_MODE = "cutpass"
- BP_F1 = 4.0
- BP_F4 = 170.0
- BUTTER_ORDER = 4
- AUTO_PICK_MODE = "hilbert_env"
- HILBERT_TARGET = "onset"
- GAIN_MODE = "agc"
- AGC_WINDOW_MS = 140.0
- AGC_STAT = "rms"

Run:

```powershell
python fb_picker.py 110
```

What to check:
1. Low-frequency waviness is reduced.
2. First breaks become easier to follow trace-to-trace.
3. If picks become unstable, lower BP_F4 (for example 150).

### Example C: Impulsive wavelet, onset alignment focus

Use when peak is strong but onset timing is needed for refraction picks.

Suggested constants:
- DEFAULT_FILTER_MODE = "butter"
- AUTO_PICK_MODE = "maxdiff_zero"
- ZERO_X_SEARCH_DIRECTION = "backward"
- MANUAL_SNAP_WIN_MS = 6.0
- PICK_PROCESS_ORDER = "F>G>X"

Run:

```powershell
python fb_picker.py 110
```

What to check:
1. Auto picks align to onset side of wavelet.
2. Manual click near wavelet body snaps toward earlier crossing.
3. If still late, increase manual backward bias by raising MANUAL_SNAP_WIN_MS slightly.

### Example D: Auto-picks missing far offsets

Use when picks appear only near source.

Suggested constants:
- AUTO_USE_VELOCITY_GATE = True
- AUTO_VMIN_M_S = 90.0
- AUTO_VMAX_M_S = 4000.0
- AUTO_GATE_PAD_MS = 30.0
- AUTO_PICK_MODE = "hilbert_env"
- HILBERT_TARGET = "onset"

Run:

```powershell
python fb_picker.py 110
```

What to check:
1. Far traces now receive candidate picks.
2. If false picks increase, tighten AUTO_VMIN_M_S upward gradually.
3. Keep pad large enough to avoid clipping real arrivals.

### Example E: Trigger/radio delay correction

Use when all picks are consistently shifted in time.

Suggested constants:
- TRIGGER_STATIC_MS = 6.0
- TRIGGER_MS_PER_M = 0.0

If delay grows with offset, test:
- TRIGGER_STATIC_MS = 2.0
- TRIGGER_MS_PER_M = 0.005

Run:

```powershell
python fb_picker.py 110
```

What to check:
1. Pick line aligns better with expected travel-time trend.
2. Near and far offsets both improve, not just one side.
3. Do not tune delay until filter and picker mode are already stable.

### Example F: Process all profiles with one setup

Use when applying one tested preset to all lines before fine adjustments.

Run:

```powershell
python fb_picker.py --all
```

If geometry must be forced:

```powershell
python fb_picker.py 150 200
```

### Fast operator checklist

1. Choose filter mode first (butter, cutpass, or ormsby).
2. Set picker mode second (hilbert_env or maxdiff_zero).
3. Verify gate settings on near and far offsets.
4. Calibrate trigger correction only after steps 1 to 3 are stable.
5. Save session often, finalize only when the line is consistent end-to-end.
