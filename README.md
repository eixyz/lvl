# LVL Project

This folder contains the Low Velocity Layer (LVL) seismic refraction workflow used to pick first breaks, fit layer velocities, estimate depths, and export standardized reports.

Main goals:
- Process SEG2 shot gathers for LVL profiles.
- Apply geometry and offset corrections with interactive quality control.
- Export reproducible outputs (plots, profile workbooks, consolidated summary).

Primary script:
- lvl/scripts/lvl_refraction.py

User instructions (non-programmer friendly, full formulas and field definitions):
- lvl/scripts/lvl_refraction_INSTRUCTIONS.md

Developer/function reference:
- lvl/scripts/lvl_refraction_CODEMAP.md

Data location:
- lvl/data/

Outputs:
- lvl/output/
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

## 14) Results export