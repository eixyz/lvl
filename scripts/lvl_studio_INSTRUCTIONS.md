# LVL Studio (Independent Program)

This is a standalone seismic workspace application (PyQt-based).

Entry point:
- python lvl/scripts/lvl_studio.py

## What makes it independent

- Does not import or call lvl_refraction.py
- Does not import or call fb_picker.py
- Uses its own GUI, plotting, and data-loading pipeline

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
8. Shot-level summary and header inspection:
   - Summary table shows shot-point metrics and selected SEG2 header fields
   - No per-receiver row list in the summary panel
9. Full computation backend bridge:
   - "Run Full Computation" writes studio picks into lvl_refraction session format
   - Executes full lvl_refraction computations and exports from GUI
   - Requirement: opened SEG2 folder must be under `lvl/data/<profile>`
6. Performance controls for large projects:
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
7. (Optional) Run full computation with the button/menu command

## Dependencies

Required packages:
- PyQt6 or PyQt5
- matplotlib
- numpy
- obspy

If missing, install into your active environment, e.g.:
- pip install PyQt6 matplotlib numpy obspy
- or: pip install PyQt5 matplotlib numpy obspy
