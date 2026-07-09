# 20241106 Low-Speed Rotating Calibration

This folder keeps the `20241106` numbering-first rotating-calibration route
separate from the gap-prior / gap-library experiments.

## Current Main Route

The current clear main route is:

```matlab
Step00_BuildLowSpeedSensorConfig_20241106
Step01_BuildLowSpeedReferenceFingerprint_20241106
Step02_NumberLowSpeedSensors_20241106
Step03_SelectHighSpeedRegion_20241106
Step04P_BuildHighSpeedPeakCache_20241106
Step04_NumberHighSpeedSensorsInRegion_20241106
Step05_BuildLowSpeedTemplateLibrary_20241106
Step06_RunIdentificationByBlade_20241106
Step07_AuditNumberingAndCorrespondence_20241106
```

Only these main-route scripts stay in the root folder.

You can run the whole route with either:

```matlab
Run_Main_OPRCenterStd_20241106
```

or

```matlab
Run_NewFlow_20241106
```

Both entries now call the same current route.

They no longer auto-run visualization helpers from the root folder.

## Step Meaning

- `Step00`: build local low-speed `Sensor_Config_20241106.mat`.
- `Step01`: build low-speed six-peak reference fingerprint.
- `Step02`: number the low-speed sensor subset.
- `Step03`: select the high-speed time region.
- `Step04P`: build reusable high-speed six-peak caches. This expensive
  whole-case peak fitting is done once and reused when scanning start times.
- `Step04`: number the high-speed sensor subset in the selected region.
- `Step05`: build the low-speed template library.
- `Step06`: directly run blade-wise identification. This script now extracts the needed high-speed waveform slices inline from raw data, so it is the main high-speed processing and identification step.
- `Step07`: audit numbering and correspondence.

## Important Update

The older standalone high-speed waveform-map route is no longer the main path.

- Old `Step06_BuildHighSpeedWaveformMap_20241106.m`
- Old `Step06V_VisualizeHighSpeedWaveformMap_20241106.m`

have been moved into:

```text
old_programs/
```

The current `Step06_RunIdentificationByBlade_20241106.m` does not require
those old `Step06` artifacts. It directly reads:

- Step04 high-speed numbering
- Step05 low-speed templates
- raw high-speed sensor data

and then builds the identification input internally.

## Tuning Rule

The preferred tuning rule is now:

- each active script has its own parameter block at the top
- do not jump back to `NewFlow_Config_20241106.m` just to tune one script

In particular:

- tune low-speed template settings in `Step05_BuildLowSpeedTemplateLibrary_20241106.m`
- tune identification settings in `Step06_RunIdentificationByBlade_20241106.m`

For scanning high-speed start times, the recommended path is now:

- run `Step04P_BuildHighSpeedPeakCache_20241106` once for the sensor set
- change `S06.startTimeSec` in `Step06_RunIdentificationByBlade_20241106.m`
- run `Step06_RunIdentificationByBlade_20241106`

`Step06` will first look for existing numbering for that time. If it is
missing, it rebuilds the current-region numbering from the cached fitted
peaks instead of refitting every pulse segment.

## Folder Layout

Root folder keeps only:

- `Run_*`
- `NewFlow_Config_20241106.m`
- `README.md`
- current main-route scripts `Step00` to `Step07`

Diagnostic / visualization / scan helpers are moved into:

```text
diagnostics/
```

Older replaced routes remain under:

```text
old_programs/
```

## Diagnostics

The following helper scripts are now under `diagnostics/`:

- low/high visualization helpers such as `Step02V`, `Step04V`, `Step05V`, `Step06M`, `Step06V`, `Step06SG`
- scan helpers such as `Step06C`, `Step06D`, `Step06L`, `Step06S`
- deep audit helpers such as `Step06E`, `Step06F`, `Step06H`, `Step07A` to `Step07E`

These scripts can still be run directly. They now auto-locate the main route
folder and load `NewFlow_Config_20241106` from there.

Figures are shown on screen. Existing saved outputs still go under:

```text
output/new_flow/figures/
```

## Output Folders

The saved artifact folder names remain:

- `05A_low_speed_template_library`
- `06_identification`
- `07_audit`

This is intentional, so existing outputs stay compatible while the script names
are now cleaner.

## Archived Old Programs

Older comparison / transition / replaced main-route programs are kept under:

```text
old_programs/
```

The main root folder should only contain the current active route.
