# 20241106 Low-Speed Gap-Prior Decoupling

This folder now keeps the `20241106` gap-prior route in the same split style
as `20250527` / `20251222`:

- `gap_tilt` is the only paper-facing main flow
- `fixed` and `gap_only` are kept only in a separate comparison script
- low-speed templates and high-speed event bases are taken from
  `20241106_low_speed_rotating_calibration` whenever possible

## Current Main Route

```matlab
Step00_CheckCase_20241106
Step00_Run_ProjectionFlow_20241106
```

`ProjectionFlow_Config_20241106.m` is the persistent parameter entry point.
Use `Step00_CheckCase_20241106` before running when changing the time,
blade, sensors, or window settings.

The recommended runner dispatches the active bank route:

```matlab
Step05_Build_Response_Surface_20241106
Step05I_Learn_OffsetTilt_Shared_Response_Surface_20241106
Step06A_BuildLowSpeedTemplateBank_20241106
Step06I_Calibrate_AllBladeGapLibrary_20241106
Step06_BuildGapAwareDynamicMap_20241106
Step07J_NestedStaticWarp_VPFullWave_20241106
```

## Comparison Route

```matlab
Step07K_Compare_NestedStaticWarp_20241106
```

It reads the same `DynamicMap`, `Template`, `CorrectedGapLibrary`, and direct
reference as `Step07J`, then compares:

- `fixed`
- `gap_only`
- `gap_tilt`

## What Each Step Does

`Step06I_Calibrate_AllBladeGapLibrary_20241106.m`

- builds one all-blade gap calibration bank for capacitive sensors S5/S7
- consumes the all-blade low-speed template bank from `Step06A`
- writes `outputs/GapCalibrationBank_20241106_B1toB6_S57.mat`

`Step06_BuildGapAwareDynamicMap_20241106.m`

- keeps the mature continuous direct-style `DynamicMap.Window(w).Sensor(s)`
  waveform structure
- attaches the matching low-speed template, corrected gap library, and direct
  baseline reference metadata
- serves as the unified event/window input for the new `Step07J` and
  `Step07K`

`Step07J_NestedStaticWarp_VPFullWave_20241106.m`

- main route only
- default main model is `gap_tilt`
- defaults to the rotating-calibration low-speed baseline template, with `Step06I_Generated_LowSpeed_Template_...` kept only as a fallback
- keeps `eoCandidateMode = direct_plus_vp`
- keeps previous-window soft EO candidate logic and phase-safe expansion
- bundle source order is strict: existing direct bundle first, then local
  `compact_case` reconstruction. DynamicMap fallback is debug-only and must
  be explicitly enabled with `STEP07J_ALLOW_DYNAMIC_MAP_FALLBACK=1`.
- writes `WindowResult`, `Trend`, `Summary`, `BestWindow`

`Step07K_Compare_NestedStaticWarp_20241106.m`

- comparison only
- reruns `Step07J` in `comparison` mode when needed
- exports method summary and window-level comparison tables
- does not change the paper-facing main route

## Direct Baseline

`20241106_low_speed_rotating_calibration/Step06_RunIdentificationByBlade_20241106.m`
is now treated as the direct-template baseline source, not the place to keep
growing the new gap-aware main flow.

## Old Programs

Older direct-only or transitional scripts stay under:

```text
old_programs/
legacy/
_archive_candidates_20260701/
```

The root folder should keep only the active gap-prior main route and the
separate comparison entry.
