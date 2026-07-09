# 20251222 low-speed template + gap-library decoupling branch

This folder is an independent 20251222 branch adapted from the 20250527
gap-prior work. It does not overwrite the existing
`20251222_low_speed_rotating_calibration` template-only workflow.

The default dynamic data source is:

```text
../20251222_low_speed_rotating_calibration
target blade: B1
analysis sensors: S123
```

## Recommended Route

Run these scripts in order:

```matlab
Step05_Build_Response_Surface_20251222
Step05H_Build_OffsetCorrected_Shared_Response_Surface_20251222
Step05I_Learn_OffsetTilt_Shared_Response_Surface_20251222
Step06I_Calibrate_OffsetTiltShared_GapLibrary_20251222
Step07J_NestedStaticWarp_VPFullWave_20251222
Step08J_CrossValidate_StaticWarp_20251222
```

The main idea is:

```text
high-speed non-vibration template
  = low-speed rotating template T_low
  + static gap-library increment
```

The low-speed rotating template remains the baseline because it carries the
real mounted rotating state. The gap library supplies only the static waveform
increment from the low-speed equivalent clearance to the high-speed equivalent
clearance.

## Static Library Logic

`Step05` builds the raw static response surface. By default it keeps the
single-reference-blade route rather than averaging all six blades.

`Step05H` and `Step05I` then use the blade-dependent initial clearance offsets
from:

```text
../../reference_blade_gap_analysis/results/analysis_04_blade_offset_summary_matched.csv
```

The intended shared static model is:

```text
g_eff = g_j + delta_g_b + mu_b x
```

Here `delta_g_b` is taken from the independent blade-gap analysis, and `mu_b`
is learned as a secondary tilt correction.

## High-Speed Models

`Step07J` compares nested high-speed static-warp models:

```text
fixed          : no high-speed static increment
gap_only       : dg_s
gap_tilt       : dg_s + dmu_s x
gap_tilt_shift : dg_s + dmu_s (x - dtau_s)
```

The paper-facing main model is now `gap_only`. The `gap_tilt` and
`gap_tilt_shift` branches are kept as ablations, because the two experimental
datasets do not support claiming a robust independent improvement from the
extra tilt/shift degrees of freedom.

VP only screens EO candidates. The final score is the full waveform objective.
No global EO consistency constraint is used.

`Step08J` cross-validates these models by refitting on held-out folds. It is
used to decide whether `gap_tilt` is genuinely useful or merely absorbs
training residuals.

## Notes

This branch has been run through the preferred 20251222 route:

```text
Step05 -> Step05H -> Step05I -> Step06I -> Step07J -> Step08J
```

Main generated files are under `outputs/`. The current default result set is
`B1_S123`.

Step07J corrected 18-window fitting shows that static gap correction strongly
reduces waveform RMSE relative to the direct low-template baseline:

```text
direct_low_template_main mean RMSE : 21.98 mV
fixed mean RMSE                    : 21.29 mV
gap_only main mean RMSE            : 14.48 mV
gap_tilt ablation mean RMSE        : 14.40 mV
gap_tilt_shift ablation mean RMSE  : 14.39 mV
```

The corrected 80-lap diagnostic gives the same practical conclusion:

```text
gap_only EO14-only RMSE            : 14.53 mV
gap_tilt EO14-only RMSE            : 14.52 mV
gap_tilt_shift EO14-only RMSE      : 14.50 mV
```

Step08J held-out cross-validation is more conservative. `gap_only` and
`gap_tilt` are essentially tied:

```text
gap_only mean validation RMSE       : 15.63 mV
gap_tilt mean validation RMSE       : 15.63 mV
gap_tilt validation gain            : 0.0065 mV
gap_tilt_shift validation RMSE      : 15.83 mV
```

Current interpretation: use `gap_only` as the robust 20251222 transfer
correction. The extra high-speed `dmu_s` tilt term is stable and not
boundary-seeking in this dataset, but its held-out gain is too weak to use as
the main model.

## EO Diagnostic Note

Coordinate diagnostics must keep the Step07J input file explicit. The validated
direct low-speed-template reference uses the GradientXRange030 template with the
base-named DynamicMap:

```text
../20251222_low_speed_rotating_calibration/output/dynamic_maps/DynamicMap_B1_S123_SlidingWindows_20251222.mat
```

Runs made with `DynamicMap_B1_S123_SlidingWindows_GradientXRange030_20251222.mat`
are useful diagnostics, but they should not be mixed with the base-DynamicMap
direct reference. A clean Step07J rerun can be saved with a suffix, for example:

```matlab
setenv('STEP07J_RESULT_SUFFIX','AuditBaseXc')
setenv('STEP07J_DYNAMIC_MAP_FILE','')
Step07J_NestedStaticWarp_VPFullWave_20251222
```

Step08J can then read that exact Step07J run without overwriting the default
cross-validation files:

```matlab
setenv('STEP08J_STEP07J_SUFFIX','AuditBaseXc')
Step08J_CrossValidate_StaticWarp_20251222
```

More-window diagnostics still show that the free EO search is weak for short windows:

```text
direct, 80 laps, 3-lap windows : EO14 = 65/78, wrong EO10 = 13
direct, 80 laps, 5-lap windows : EO14 = 68/76, wrong EO10 = 8
gap-prior, 80 laps, 3-lap windows:
  direct  : EO14 = 65/78, wrong EO10 = 13
  fixed   : EO14 = 68/78, wrong EO10/EO3/EO13 = 6/3/1
  gap     : EO14 = 73/78, wrong EO3/EO7 = 3/2
  tilt    : EO14 = 72/78, wrong EO3/EO7 = 4/2
  shift   : EO14 = 73/78, wrong EO3/EO7 = 3/2
```

The 80-lap corrected comparison is saved as:

```text
Step07J_MethodComparison_main80_gradient_xc_20251222_B1_S123.csv
Step07J_MethodBars_main80_gradient_xc_20251222_B1_S123.png
Step07J_WindowTrends_main80_gradient_xc_20251222_B1_S123.png
```

This indicates that the main failure source is not simply the static gap
library. The unconstrained EO objective has close local alternatives in short
windows. `gap` correction removes several direct EO10 mistakes, but one window
still flips to EO7. For paper-facing results, use global EO consistency
(`EO=14`) or an EO-continuity/lock step before comparing amplitude and waveform
RMSE.

The older `Step06F/Step06G/Step07C/Step07D/Step07I` scripts were also ported as
comparison branches, but the current preferred diagnostic path is
`Step05I -> Step06I -> Step07J -> Step08J`.
