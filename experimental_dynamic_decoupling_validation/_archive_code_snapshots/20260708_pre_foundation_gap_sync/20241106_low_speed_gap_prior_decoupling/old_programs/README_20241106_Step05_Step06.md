# 20241106 Step05 main workflow

20241106 now follows the same response-surface workflow as the 20250527 and
20251222 folders.

## Run order

```matlab
Step05_Build_Response_Surface_20241106
Step05I_Learn_OffsetTilt_Shared_Response_Surface_20241106
```

## What Step05 does

`Step05_Build_Response_Surface_20241106.m` reads the staged static gap
waveform library:

```text
experimental_dynamic_decoupling_validation/static_gap_waveform_library_20260524_2000Hz
```

It fits the static gap response surface:

```text
F(g,x) = B0(x) + B1(x)/g + B2(x)*log(g/g0)
```

Main outputs:

```text
outputs/Step05_Response_Surface_20241106.mat
outputs/Step05_Response_Surface_Fit_20241106.csv
outputs/Step05_Calibration_Motion_20241106.csv
outputs/Step05_Static_Gap_Continuous_Baseline_20241106.csv
```

The figure opens directly when the script runs, so the fitted waveform,
surface, fit error and sensitivity can be checked visually.

## What Step05I does

`Step05I_Learn_OffsetTilt_Shared_Response_Surface_20241106.m` loads the Step05
response surface and the blade offset table:

```text
reference_blade_gap_analysis/results/analysis_04_blade_offset_summary_matched.csv
```

It keeps the blade gap offset fixed and learns one tilt coefficient for each
blade:

```text
g_eff,bj(x) = g_j + delta_g_b + mu_b*x
```

Main outputs:

```text
outputs/Step05I_OffsetTilt_Shared_Response_Surface_20241106.mat
outputs/Step05I_OffsetTilt_Shared_Response_Surface_Fit_20241106.csv
outputs/Step05I_OffsetTilt_Shared_Response_Surface_Tilt_20241106.csv
outputs/Step05I_OffsetTilt_Shared_Response_Surface_BladeSummary_20241106.csv
```

## Relation to Step04E

Step04E is used to find the large resonance regions from OPR-aligned strain
STFT and BTT vibration context. Step05 does not select BTT time windows. It
prepares the static waveform model needed by the later dynamic identification
steps.

`Step05_Select_BTT_Windows_From_Resonance_20241106.m` is kept only as an
optional diagnostic script. It should not be treated as the formal Step05 main
workflow, because Step02 BTT displacement is an instantaneous blade-pass
displacement sample, not a direct vibration-amplitude response for one strain
gauge blade.
