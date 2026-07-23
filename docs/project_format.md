# Project folder format

A "project" is one field survey. It can live anywhere on disk (by default
under `projects/<name>/` in the app install folder, but the "New Project"
dialog lets you pick any location). There is no more flat `data/` layout -
every project owns its own `processing/` and `results/` folders, and
remembers where its raw/geometry/metadata inputs live externally.

```
<project root>/
    project.json

    processing/
        cache/
        picks/<profile>/picks.json
        sessions/<profile>/picks.session.json
                          /layer_analysis.session.json

    results/
        plots/<profile>/...
        reports/<profile>/layer_analysis.json
        velocity/lvl_velocity_summary.xlsx
```

Implemented in `src/common/paths.py` (`ProjectPaths`) and
`src/io/project_io.py` (`Project`, `create_project`, `open_project`,
`save_project`, `list_projects`).

## project.json

```json
{
  "name": "LVL_Erfurt",
  "raw_folder": "D:\\Daten\\seismic\\lvl\\data\\input\\raw",
  "geometry_folder": "D:\\Daten\\seismic\\lvl\\data\\input\\geometry",
  "metadata_folder": "D:\\Daten\\seismic\\lvl\\data\\input\\metadata",
  "manual_geometry_files": [],
  "profiles": ["110", "120", "125"],
  "notes": "",
  "created": "2026-07-20T10:19:52+00:00",
  "modified": "2026-07-21T08:04:36+00:00"
}
```

- **raw_folder** - the folder that contains *one subfolder per profile*,
  each holding that profile's SEG2 files (e.g. `raw_folder/125/Rec_*.seg2`).
  It must be the *parent* of the profile folders, not a profile folder
  itself - "Project Settings..." now warns you if you pick a folder that
  has SEG2 files directly in it, since that usually means you picked one
  profile's folder by mistake.
- **geometry_folder** - where survey coordinate tables live
  (`LVL*.xlsx/.xls/.csv/.txt/.dat`, matched to a profile by a "Profile" /
  "LVL Number" column - see `src/utils/geometry.py`). Not required to be
  set; `manual_geometry_files` below covers files that live elsewhere.
- **metadata_folder** - where field-report workbooks
  (`field_report*.xls(x)`) live, for perpendicular-offset auto-detection.
- **manual_geometry_files** - individual geometry/coordinate files you
  explicitly registered via **File > Add Geometry File(s)...**, for cases
  where a file isn't under `geometry_folder` or doesn't follow the
  `LVL*` naming convention. Always searched, and takes priority over the
  auto-discovered folder scan.
- **profiles** - profiles that have been opened/picked at least once in
  this project (auto-updated; informational, not required for anything
  to work).

None of `raw_folder` / `geometry_folder` / `metadata_folder` are copied
into the project - they're referenced by path only, since raw survey data
is often large and lives on a field laptop or network share.

## Application-level resources (not part of any project)

```
resources/
    geometry_templates/
        geometry100.txt      # 100 m spread receiver geometry
        geometry200.txt      # 200 m spread receiver geometry
    assets/
        lvl_logo.png         # About dialog / header
        lvl_icon.png         # window/taskbar icon
```

These are fixed, bundled with the app itself (see `src/common/paths.py`:
`RESOURCES_DIR`, `GEOM_TEMPLATES_DIR`, `ASSETS_DIR`) - swap
`resources/assets/lvl_logo.png` / `lvl_icon.png` for the real company logo
at any time; missing files are handled gracefully (app just runs without
a custom icon).