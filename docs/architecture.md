# Architecture

```
GUI (lvl_studio.py, PyQt)     CLI (lvl_refraction.py)     (future: web backend)
        \                              |                              /
         \                             |                             /
          '---------------  src/  (framework-agnostic core)  -------'
```

Both front-ends call into `src/` for anything that isn't pure UI - no
picking, geometry, or export logic should live only in one app. This is
what makes a future web version possible without duplicating the
science: the web backend would be a third caller of the same `src/`
modules, same as the desktop GUI and the CLI are today.

## Layers

- **`src/common/`** - `paths.py` (project system, see `project_format.md`),
  `settings.py` (shared tunables, app identity).
- **`src/io/`** - `project_io.py` (create/open/save projects),
  `seg2_reader.py`, `pick_reader.py`, `exporters.py` (Excel/plot output).
- **`src/utils/geometry.py`** - receiver geometry templates, profile
  discovery, and matching survey coordinate tables (Excel/csv/txt/dat) to
  a profile by name.
- **`src/picker/`** - the "Picker Engine V2" pipeline (see
  `picker_v2_design.md` if present, or the module docstrings):
  `preprocessing -> features -> likelihood -> coherence -> optimizer ->
  confidence -> quality`, all orchestrated through the single public
  class `picker.FirstBreakPicker`. GUI/CLI code should only ever call
  `FirstBreakPicker`, never the individual pipeline modules.
- **`src/refraction/`** - layer analysis, travel-time processing, velocity/
  depth model computation, and `pipeline.py` (`process_profile`,
  `AnalysisWorkflow`) - the "run everything for this profile" orchestrator
  shared by both front-ends. Both `lvl_studio.py`'s "Run Full Computation"
  and `lvl_refraction.py`'s CLI call the exact same functions here; neither
  has its own copy. `lvl_studio.py` has no import-time or runtime
  dependency on `lvl_refraction.py` at all anymore (verified: instantiating
  the GUI never loads the `lvl_refraction` module).

  Two genuinely interactive, desktop-only steps are lazily imported from
  `lvl_refraction.py` inside `pipeline.py` rather than duplicated: the v1
  matplotlib picking window (`FirstBreakPicker`, used when picking
  interactively rather than with `--auto-pick`/Picker V2) and the
  PO-entry / layer-window-picking prompts inside `AnalysisWorkflow.run()`.
  These are exactly the seams a future web backend would replace with its
  own UI flow - everything else in `pipeline.py` is already reusable as-is.

## Per-project output

Every project (see `project_format.md`) gets its own:

```
<project>/results/
├── velocity/
│   └── lvl_velocity_summary.xlsx
├── reports/
│   └── <profile>/
│       └── layer_analysis.json
└── plots/
    └── <profile>/
        ├── arrivals.png
        └── rms.png
```

This is per-project now, not a single shared `data/results/` tree - two
different surveys never collide, and opening a different project just
means a different `results/` folder underneath it.