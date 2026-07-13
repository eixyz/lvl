# LVL Studio (Independent Program)

This is a standalone seismic workspace application (PyQt-based).

Entry point:
- python lvl/scripts/lvl_studio.py

## Architecture

- Dedicated PyQt GUI application and display workflow.
- Uses lvl_refraction core functions as backend for processing, picking modes, and analysis math.
- Does not use fb_picker.py.

## Core capabilities

1. Open any SEG2 folder with arbitrary shot count.
2. Handles arbitrary station/trace count per shot (not limited to 48).
3. Display-first layout with resizable side docks:
   - Controls dock (compact)
   - Shot Browser dock
   - Shot Summary dock
4. Undock and maximize the display to a separate window, then dock it again.
5. Built-in zoom/pan/home via matplotlib navigation toolbar.
6. Display modes:
   - wiggle
   - vd (variable density)
   - both
7. Native picker:
   - Left-click: set pick on selected trace
   - Right-click: remove pick
   - Save picks to JSON
8. Processing controls mirrored from lvl_refraction:
   - SEG2 delay usage, TRIGGER_STATIC_MS and extra delay
   - Polarity normal/invert
   - Gain modes none/norm/agc with AGC window/stat
   - Filter modes none/butter/ormsby/cutpass with frequency controls
   - Auto-pick modes stalta/maxdiff_zero/hilbert_env
9. Shot-level summary and header inspection:
   - Summary table shows shot-point metrics and selected SEG2 header fields
   - No per-receiver row list in the summary panel
10. Full computation backend bridge:
   - "Run Full Computation" writes studio picks into lvl_refraction session format
   - Executes full lvl_refraction computations and exports from GUI
   - Requirement: opened SEG2 folder must be under `lvl/data/<profile>`
11. Interactive analysis bridge:
   - Runs perpendicular-offset setup and interactive layer-window picking from GUI
   - Shows average velocity summary after analysis
12. Performance controls for large projects:
   - Max traces in view (decimation for rendering)
   - Clip percentile
   - Colormap selection

## Typical workflow

1. Start app: python lvl/scripts/lvl_studio.py
2. Click "Open SEG2 Folder"
3. Select project profile folder
4. Browse shots from left dock
5. Undock display from View > Undock Display for full-screen interpretation
6. Re-dock with View > Dock Display
7. Pick arrivals manually or with Auto Pick mode
8. (Optional) Run interactive analysis for perpendicular offsets/layer picks
9. (Optional) Run full computation with export

## Dependencies

Required packages:
- PyQt6 or PyQt5
- matplotlib
- numpy
- obspy

If missing, install into your active environment, e.g.:
- pip install PyQt6 matplotlib numpy obspy
- or: pip install PyQt5 matplotlib numpy obspy
